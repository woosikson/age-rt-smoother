#=
age-rt-smoother / code/01_synthetic.jl
STEP 1 - synthetic age-structured SEIR epidemic with a prescribed R_j(t)  (paper Fig 2, Table S5)

- Generates the synthetic epidemic used to validate the estimator, together with the
  prescribed (true) age-specific instantaneous reproduction numbers.

- Model: age-structured SEIR with six bands (0-5, 6-11, 12-17, 18-44, 45-64, 65+)

    dS_i/dt = -chi_i(t) (S_i/N_i) sum_j c_ij I_j
    dE_i/dt =  chi_i(t) (S_i/N_i) sum_j c_ij I_j - sigma E_i
    dI_i/dt =  sigma E_i - gamma I_i
    dR_i/dt =  gamma I_i

  sigma = gamma = 1/1.8 per day (mean latent and infectious periods of 1.8 days). The latent
  stage shifts timing only and leaves the reproduction number unchanged, because the
  next-generation matrix of this SEIR model equals that of the corresponding SIR model.

- True instantaneous reproduction number = column sum of the time-dependent next-generation
  matrix,  K_ij(t) = chi_i(t) c_ij / gamma * S_i(t)/N_i ,  R_j(t) = sum_i K_ij(t).

- Prescribed transmission probability per contact:
      chi_i(t) = chi_b,i + chi_m,i * sin(2*pi/T_i * (t - phi_i)),   chi_m = alpha * chi_b
  The parameters (Table S5) were tuned so that both the peak height per 100,000 and the peak
  timing follow  6-11 > 0-5 > 12-17 > 18-44 > 45-64 > 65+ :
  chi_b sets the level and ordering, phi staggers the timing, and a single 300-day period
  gives one hump over the simulated window.
  Note that 18-44 is a contact hub and therefore peaks within a day of 12-17; this near-tie
  is structural and cannot be removed by tuning.

- Initial condition: I(0) = E(0) = [597, 453, 188, 677, 302, 223] (integer seeds),
  R(0) = 0, S(0) = N - E(0) - I(0), simulated for T = 150 days.

- Age bands follow the paper. Population N_i is aggregated from data/Pops_Dec2022.csv
  (Ministry of the Interior and Safety, resident registration, December 2022, national row);
  the contact matrix c_ij is the survey-based matrix of Son et al. (2025).

- No figures are drawn here (see code/figures.R). Output, in tidy long form:
    results/pre_confirm.csv : day, age_group, incidence   (daily new recoveries, i.e. cases)
    results/pre_rt.csv      : day, age_group, Rt          (true instantaneous R_j(t))
  plus results/pre_confirm_noisy.csv, an overdispersed replicate used for the weight study.

- Run: julia --project=. code/01_synthetic.jl
=#

using DifferentialEquations
using CSV
using DataFrames
using Printf
using Distributions, Random         # for the overdispersed (NB) observation replicate

const AGE = ["0-5", "6-11", "12-17", "18-44", "45-64", "65+"]

# --- population N_i: Pops_Dec2022.csv (EUC-KR, national row) aggregated to the six bands ---
popfile = joinpath(normpath(joinpath(@__DIR__,"..")), "data", "Pops_Dec2022.csv")
popdf = CSV.read(popfile, DataFrame; header=1)
pops = [parse(Int, replace(String(popdf[1, c]), "," => "")) for c in 4:ncol(popdf)]  # age 0..100+
lb = [1, 7, 13, 19, 46, 66]                     # 1-indexed: 0-5, 6-11, 12-17, 18-44, 45-64, 65+
ub = [6, 12, 18, 45, 65, length(pops)]
N = Float64[sum(pops[lb[i]:ub[i]]) for i in 1:6]
@assert sum(N) == sum(pops)
println("N (paper bands): ", join(["$(AGE[i])=$(Int(N[i]))" for i in 1:6], ", "))

# --- contact matrix c_ij (survey-based, Son et al. 2025); column = infector ---
C = [1.76 0.21 0.03 0.20 0.05 0.05;
     0.33 3.75 0.30 0.31 0.14 0.13;
     0.05 0.31 3.65 0.22 0.31 0.12;
     2.07 2.11 1.38 1.78 1.26 0.84;
     0.47 0.88 1.87 1.18 1.97 1.47;
     0.30 0.47 0.42 0.45 0.83 2.50]

# ─── SEIR parameters ───
σ = 1 / 1.8
γ = 1 / 1.8

# --- prescribed chi_i(t): seasonal sine, level set by chi_b and timing by phi ---
# With gamma = 1/1.8 rather than 1/5, the same chi would scale R_j by 1.8/5 and no epidemic
# would occur, so chi is scaled by 5/1.8 to restore the intended level of R_j.
# Tuning target: both the peak height per 100,000 and the peak timing must follow
#   6-11 > 0-5 > 12-17 > 18-44 > 45-64 > 65+ .
# chi_b controls level and ordering, phi controls the timing stagger.
# Caveat: 18-44 is a contact hub (c = 2.07 and 2.11 with 0-5 and 6-11), so it rides the
# crest of the young groups and peaks within a day of the largely self-contacting 12-17 group.
T = 150   # length of the simulated series in days, t = 1..150
# chi_4 raised and chi_6 lowered so that the three adult curves are visually separable;
# this also breaks the timing tie, making both orderings strict.
χ_base = [0.200, 0.185, 0.120, 0.135, 0.075, 0.035]   # 0-5, 6-11, 12-17, 18-44, 45-64, 65+
χ_amp  = 0.45                                          # relative seasonal amplitude (chi_m = alpha * chi_b)
χ_max  = χ_amp .* χ_base
period = fill(300.0, 6)
φ      = [-20.0, -40.0, -10.0, 30.0, 45.0, 60.0]       # phase stagger (6-11 earliest ... 65+ latest)
χ = zeros(6, T)
for i in 1:6, t in 1:T
    χ[i, t] = χ_base[i] + χ_max[i] * sin(2π / period[i] * (t - φ[i]))
end

# --- initial condition ---
# integer seeds, so that the initial condition can be stated exactly in the paper:
#   round(5·[119.327,90.5102,37.6122,135.469,60.4898,44.6327]) = [597,453,188,677,302,223]
I0 = round.(5 .* [119.327, 90.5102, 37.6122, 135.469, 60.4898, 44.6327])
E0 = copy(I0)               # E(0) = I(0)
R0 = zeros(6)
S0 = N .- E0 .- I0

# --- one-day propagation; chi_i is held constant within each day (piecewise constant) ---
# state u = [S(1:6); E(7:12); I(13:18); R(19:24)]
function seir_ag!(du, u, p, t)
    S = @view u[1:6]; E = @view u[7:12]; I = @view u[13:18]
    foi = p.χt .* (S ./ p.N) .* (p.C * I)          # χ_i (S_i/N_i) Σ_j c_ij I_j
    @views @. du[1:6]   = -foi
    @views @. du[7:12]  =  foi - p.σ * E
    @views @. du[13:18] =  p.σ * E - p.γ * I
    @views @. du[19:24] =  p.γ * I
    return nothing
end

function propagate_day(u0, χt)
    prob = ODEProblem(seir_ag!, u0, (0.0, 1.0), (χt=χt, C=C, N=N, σ=σ, γ=γ))
    sol = solve(prob, Tsit5(); save_everystep=false, abstol=1e-9, reltol=1e-9)
    return sol.u[end]
end

# state history: S/E/I/R as 6 x T, with R cumulative
function simulate(S0, E0, I0, R0, χ, T)
    S = zeros(6, T); E = zeros(6, T); Icomp = zeros(6, T); Rcum = zeros(6, T)
    S[:, 1] = S0; E[:, 1] = E0; Icomp[:, 1] = I0; Rcum[:, 1] = R0
    u = vcat(S0, E0, I0, R0)
    for t in 1:(T-1)
        u = propagate_day(u, χ[:, t])              # chi held constant during day t
        S[:, t+1] = u[1:6]; E[:, t+1] = u[7:12]; Icomp[:, t+1] = u[13:18]; Rcum[:, t+1] = u[19:24]
    end
    return S, E, Icomp, Rcum
end
S, E, Icomp, Rcum = simulate(S0, E0, I0, R0, χ, T)

# incidence (cases) = daily new recoveries, dR over day t = Rcum[t+1] - Rcum[t]
incid = Rcum[:, 2:T] .- Rcum[:, 1:(T-1)]           # 6 × (T-1), day 1..T-1

# true instantaneous 𝓡_j(t) = Σ_i (χ_i(t)/γ)(S_i(t)/N_i) c_ij   (day 1..T-1)
Rt = zeros(6, T-1)
for t in 1:(T-1), j in 1:6
    s = 0.0
    for i in 1:6
        s += (χ[i, t] / γ) * (S[i, t] / N[i]) * C[i, j]
    end
    Rt[j, t] = s
end

# ─── tidy long CSV ───
dn = joinpath(normpath(joinpath(@__DIR__,"..")), "results"); isdir(dn) || mkpath(dn)
# incidence is an absolute count (dR), which the observation model consumes directly;
# the per-1,000 and per-100,000 columns are for display only.
open(joinpath(dn, "pre_confirm.csv"), "w") do io
    println(io, "day,age_group,incidence,incidence_per1000,incidence_per100k")
    for t in 1:(T-1), i in 1:6
        @printf(io, "%d,%s,%.6f,%.6f,%.6f\n", t, AGE[i], incid[i, t],
                incid[i, t] * 1e3 / N[i], incid[i, t] * 1e5 / N[i])
    end
end
open(joinpath(dn, "pre_rt.csv"), "w") do io
    println(io, "day,age_group,Rt")
    for t in 1:(T-1), j in 1:6
        @printf(io, "%d,%s,%.6f\n", t, AGE[j], Rt[j, t])
    end
end

# --- overdispersed replicate for the observation-model study: Y ~ NB(mean = dR, size = r_gen).
#   Only the observation is perturbed; the true R_j is unchanged. Var = mu + mu^2/r_gen.
const R_GEN = 20                                    # true overdispersion (about 22% relative noise)
Random.seed!(20260808)                              # fixed seed: a single reproducible realization
incid_noisy = similar(incid)
for t in 1:(T-1), i in 1:6
    μ = max(incid[i, t], 1e-8)
    incid_noisy[i, t] = rand(NegativeBinomial(R_GEN, R_GEN / (R_GEN + μ)))   # mean μ, size r_gen
end
open(joinpath(dn, "pre_confirm_noisy.csv"), "w") do io
    println(io, "day,age_group,incidence,incidence_per1000,incidence_per100k")
    for t in 1:(T-1), i in 1:6
        @printf(io, "%d,%s,%.6f,%.6f,%.6f\n", t, AGE[i], incid_noisy[i, t],
                incid_noisy[i, t] * 1e3 / N[i], incid_noisy[i, t] * 1e5 / N[i])
    end
end

@printf("Done. peak incidence per age: %s | noisy(NB r=%d) peak: %s\n",
        join([@sprintf("%s=%.0f", AGE[i], maximum(incid[i, :])) for i in 1:6], ", "), R_GEN,
        join([@sprintf("%s=%.0f", AGE[i], maximum(incid_noisy[i, :])) for i in 1:6], ", "))
@printf("true 𝓡_j(t=1) range: %.2f–%.2f  | overall Rt range %.2f–%.2f\n",
        minimum(Rt[:, 1]), maximum(Rt[:, 1]), minimum(Rt), maximum(Rt))
println("CSV -> results/{pre_confirm,pre_rt}.csv  |  next: code/02_particle_smoother.jl")
