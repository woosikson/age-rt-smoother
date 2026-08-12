#=
age-rt-smoother / code/02_particle_smoother.jl
STEP 2 - age-structured SEIR particle smoother  (paper Fig 3, S3, S4; Tables S1-S4)

- Estimates the age-specific instantaneous reproduction number R_j(t) and the model incidence
  dR_i(t) from the synthetic cases written by code/01_synthetic.jl (results/pre_confirm.csv),
  and computes the diagnostics reported in the supplement.

- Method: per-stratum multinomial resampling (age groups resampled independently), an NNLS
  initial guess for chi, a log-random-walk on chi with scale sigma_rw, and a fixed-lag window L
  (whole-history resampling is available for comparison).

- The initial state is a random variable: E(0), I(0) ~ iid Poisson(mean = true E(0), I(0)) per
  particle, R(0) = 0, S(0) = N - E(0) - I(0). Treating it as a point mass would understate the
  uncertainty of the first days.

- Observation model: Y_i[t] = dR_i (daily new recoveries, absolute counts) ~ Poisson or
  NB(r) with mean equal to the model dR_i of each particle.

- R_j(t) is the column sum of the time-dependent next-generation matrix; the latent stage
  cancels, so the SEIR and SIR next-generation matrices coincide:
    K_ij(t) = chi_i(t) c_ij / gamma * S_i(t)/N_i ,   R_j(t) = sum_i K_ij(t).

- ODE solver: Tsit5 with one fixed step per day.

- Run: julia --project=. --threads=auto code/02_particle_smoother.jl --mode sweep --weight-r 30
       (see the CLI block below for the other modes)
=#

using DifferentialEquations
using CSV, DataFrames
using Distributions, StatsBase, LinearAlgebra, NonNegLeastSquares
using Random, Printf

# ── CLI options (the corresponding environment variables also work) ─────────
#   --mode      sweep | weight_sweep | weight_misspec | weight_curves_r25 |
#               save_noisy | baseline_fit | pfps | test        (default: test)
#   --weight-r  R      observation model NB(r); 0 = Poisson
#   --noise-r   R      overdispersion of the synthetic observation noise; 0 = clean
#   --np-list   comma-separated particle counts for the sweep, e.g. 10000,100000,1000000
#   --np / --window / --sigma-rw : single run in test mode
let opt = Dict("--mode" => "RT_PS_MODE", "--weight-r" => "RT_NBR", "--noise-r" => "RT_NOISY_RGEN",
               "--np-list" => "NPS", "--np" => "NP", "--window" => "L", "--sigma-rw" => "SRW",
               "--sigma-list" => "SRWS", "--baseline-threshold" => "BASETHR")
    i = 1
    while i <= length(ARGS)
        haskey(opt, ARGS[i]) || error("unknown option: $(ARGS[i])")
        i < length(ARGS) || error("missing value for $(ARGS[i])")
        ENV[opt[ARGS[i]]] = ARGS[i + 1]
        i += 2
    end
end

const AGE = ["0-5", "6-11", "12-17", "18-44", "45-64", "65+"]
const C = [1.76 0.21 0.03 0.20 0.05 0.05;
           0.33 3.75 0.30 0.31 0.14 0.13;
           0.05 0.31 3.65 0.22 0.31 0.12;
           2.07 2.11 1.38 1.78 1.26 0.84;
           0.47 0.88 1.87 1.18 1.97 1.47;
           0.30 0.47 0.42 0.45 0.83 2.50]
const σ = 1 / 1.8
const γ = 1 / 1.8

# --- population N_i ---
popdf = CSV.read(joinpath(normpath(joinpath(@__DIR__,"..")), "data", "Pops_Dec2022.csv"), DataFrame; header=1)
pops = [parse(Int, replace(String(popdf[1, c]), "," => "")) for c in 4:ncol(popdf)]
lb = [1, 7, 13, 19, 46, 66]; ub = [6, 12, 18, 45, 65, length(pops)]
const N = Float64[sum(pops[lb[i]:ub[i]]) for i in 1:6]

# --- synthetic observations Y (= dR, integer counts) and the true R_j(t) ---
conf = CSV.read(joinpath(normpath(joinpath(@__DIR__,"..")), "results", "pre_confirm.csv"), DataFrame)
rtdf = CSV.read(joinpath(normpath(joinpath(@__DIR__,"..")), "results", "pre_rt.csv"), DataFrame)
const Tobs = maximum(conf.day)
const Y = round.(Int, reshape(conf.incidence, 6, Tobs))    # (6, Tobs) absolute counts
const trueRj = reshape(rtdf.Rt, 6, Tobs)                   # (6, Tobs)
const trueDR = reshape(conf.incidence, 6, Tobs)            # true ΔR (float)

# Initial condition: these are the MEANS of the Poisson prior (equal to the true values used by
#   code/01_synthetic.jl). E(0) and I(0) are drawn per particle inside run_ps, so the initial state
const I0 = [597.0, 453, 188, 677, 302, 223]      # Poisson mean (= true I(0))
const E0 = copy(I0)                               # Poisson mean (= true E(0))
const R0 = zeros(6)
const S0 = N .- E0 .- I0                          # S(0) at the mean, used for the NNLS matrix only

# --- batched SEIR right-hand side (u: Np x 24 = [S(1:6) E(7:12) I(13:18) R(19:24)]; p.chi: Np x 6) ---
function seir_rhs!(du, u, p, t)
    S = @view u[:, 1:6]; E = @view u[:, 7:12]; I = @view u[:, 13:18]
    mul!(p.CI, I, p.Ct)                              # CI[p,i] = Σ_j c_ij I[p,j]
    @. p.foi = p.χ * (S / p.Nr) * p.CI               # foi = χ_i (S_i/N_i) Σ_j c_ij I_j
    @views @. du[:, 1:6]   = -p.foi
    @views @. du[:, 7:12]  =  p.foi - σ * E
    @views @. du[:, 13:18] =  σ * E - γ * I
    @views @. du[:, 19:24] =  γ * I
    return nothing
end

"""Advance one day: state u0 (Np x 24) and chi (Np x 6) -> next day. Buffers are reused."""
function step_day!(u0, χ, Ct, Nr, CI, foi)
    p = (χ=χ, Ct=Ct, Nr=Nr, CI=CI, foi=foi)
    prob = ODEProblem(seir_rhs!, u0, (0.0, 1.0), p)
    sol = solve(prob, Tsit5(); dt=1.0, adaptive=false, save_everystep=false)
    return sol.u[end]
end

# --- NNLS initial chi: A[j,i] = c_ij/gamma * S0_i/N_i,  chi_p = nnls(A, Rtarget_p) ---
function init_chi(Np, σ_rw; R0target=1.05)
    A = permutedims((C .* (S0 ./ N)) ./ γ)           # A[j,i] = C[i,j]·S0_i/N_i/γ
    χ0 = Matrix{Float64}(undef, Np, 6)
    for p in 1:Np
        Rtarget = R0target .* exp.(randn(6) .* σ_rw)
        χ0[p, :] = nonneg_lsq(A, Rtarget; alg=:nnls)
    end
    return χ0
end

# ─── particle smoother (per-stratum resample, Poisson weight, fixed-lag L / whole-history) ───
# adaptive = false: resample every step. adaptive = true: only when ESS_i < ess_frac * Np.
# Returns R_j and dR (median, 2.5%, 97.5%), ESS(6,T), unique ancestors(6,T), resample flags, runtime.
function run_ps(Np, L, σ_rw; adaptive=false, ess_frac=0.5, seed=1, Yobs=Y, nb_r=0, robust=false)
    Random.seed!(seed)
    T = Tobs; whole = (L >= T)
    t0 = time()
    # trajectories: [particle, time (1 = day 0, 2 = day 1, ..., T+1 = day T), age group]
    Xh = Array{Float64}(undef, Np, T + 1, 6)         # χ
    Sh = Array{Float64}(undef, Np, T + 1, 6)         # S
    Rh = Array{Float64}(undef, Np, T + 1, 6)         # R (cumulative)
    Ih = zeros(Float64, Np, T + 1, 6)                # model incidence dR, stored along the trajectory
    Ecur = Matrix{Float64}(undef, Np, 6)
    Icur = Matrix{Float64}(undef, Np, 6)
    # initialisation
    χ0 = init_chi(Np, σ_rw)
    @views Xh[:, 1, :] .= χ0
    # The initial state is a random variable: E(0), I(0) ~ iid Poisson(mean) per particle and group.
    for i in 1:6
        @views Ecur[:, i] .= rand(Poisson(E0[i]), Np)
        @views Icur[:, i] .= rand(Poisson(I0[i]), Np)
        @views Sh[:, 1, i] .= N[i] .- Ecur[:, i] .- Icur[:, i] .- R0[i]
        @views Rh[:, 1, i] .= R0[i]
    end
    # logs and diagnostics
    W = fill(1.0 / Np, Np, 6)                        # accumulated normalised weights, used when adaptive
    logW = zeros(Np, 6)
    ess = zeros(6, T); uniq = zeros(Int, 6, T); rflag = falses(6, T)
    # buffers
    Ct = permutedims(C); Nr = reshape(N, 1, 6)
    CI = Matrix{Float64}(undef, Np, 6); foi = Matrix{Float64}(undef, Np, 6)
    u = Matrix{Float64}(undef, Np, 24)

    for t in 1:T
        ti = t + 1                                   # time index in history
        # 1) χ log-random-walk:  χ_i(t) = χ_i(t-1)·exp(N(0,σ_rw))
        χt = @view Xh[:, ti, :]
        @inbounds for a in 1:6, p in 1:Np
            χt[p, a] = Xh[p, ti-1, a] * exp(σ_rw * randn())
        end
        # 2) propagate SEIR one day
        @views u[:, 1:6]   .= Sh[:, ti-1, :]
        @views u[:, 7:12]  .= Ecur
        @views u[:, 13:18] .= Icur
        @views u[:, 19:24] .= Rh[:, ti-1, :]
        unext = step_day!(u, χt, Ct, Nr, CI, foi)
        @views Sh[:, ti, :] .= unext[:, 1:6]
        @views Ecur          .= unext[:, 7:12]
        @views Icur          .= unext[:, 13:18]
        @views Rh[:, ti, :] .= unext[:, 19:24]
        # 3) model dR_i = R(ti) - R(ti-1), then weight and resample within each stratum
        for i in 1:6
            dR = @views Rh[:, ti, i] .- Rh[:, ti-1, i]
            @views Ih[:, ti, i] .= dR         # keep the incidence aligned with the resampled path
            yi = Yobs[i, t]
            # nb_r = 0 gives Poisson, nb_r > 0 gives NB(mean = model dR); diverged particles get weight 0
            if robust
                # log-space weights (max subtraction) for noisy data: a Poisson likelihood can underflow
                # for every particle at an outlier, leaving sum(w) = 0. Normalised weights are unchanged.
                lw = map(dR) do d
                    isfinite(d) || return -Inf
                    dd = max(d, 1e-8)
                    nb_r > 0 ? logpdf(NegativeBinomial(nb_r, nb_r / (nb_r + dd)), yi) : logpdf(Poisson(dd), yi)
                end
                m = maximum(lw)
                w = isfinite(m) ? exp.(lw .- m) : fill(1.0, Np)   # all -Inf -> fall back to uniform
            else                                            # clean data: plain linear weights
                w = map(dR) do d
                    isfinite(d) || return 0.0
                    dd = max(d, 1e-8)
                    nb_r > 0 ? pdf(NegativeBinomial(nb_r, nb_r / (nb_r + dd)), yi) : pdf(Poisson(dd), yi)
                end
            end
            if adaptive
                @. logW[:, i] += log(w + 1e-300)
                m = maximum(@view logW[:, i]); wn = exp.(@view(logW[:, i]) .- m)
                sw = sum(wn); wn = sw > 0 ? wn ./ sw : fill(1/Np, Np)
                @views W[:, i] .= wn
            else
                sw = sum(w); wn = sw > 0 ? w ./ sw : fill(1/Np, Np)
                @views W[:, i] .= wn
            end
            essi = 1.0 / sum(abs2, @view W[:, i]); ess[i, t] = essi
            doresample = adaptive ? (essi < ess_frac * Np) : true
            if doresample
                rflag[i, t] = true
                inds = wsample(1:Np, @view(W[:, i]), Np; replace=true)
                win = whole ? (1:ti) : (max(ti - L, 1):ti)
                @views Xh[:, win, i] .= Xh[inds, win, i]
                @views Sh[:, win, i] .= Sh[inds, win, i]
                @views Rh[:, win, i] .= Rh[inds, win, i]
                @views Ih[:, win, i] .= Ih[inds, win, i]   # resample the incidence path with the same indices
                @views Ecur[:, i]    .= Ecur[inds, i]
                @views Icur[:, i]    .= Icur[inds, i]
                if adaptive; @views logW[:, i] .= 0.0; end
            end
        end
    end
    # degeneracy: number of distinct ancestors (distinct chi values) per day in the FINAL trajectory,
    for t in 1:T, i in 1:6
        uniq[i, t] = length(unique(@view Xh[:, t+1, i]))
    end

    # --- summaries (per-stratum weights W; uniform when resampling every step) ---
    Rj = zeros(3, 6, T); DR = zeros(3, 6, T)         # [med,lo,hi] × group × day
    Rjp = Matrix{Float64}(undef, Np, 6)
    for t in 1:T
        ti = t + 1
        # 𝓡_j^{(p)}(t) = Σ_i (χ_i/γ)(S_i/N_i) c_ij
        fill!(Rjp, 0.0)
        @inbounds for j in 1:6, i in 1:6
            cij = C[i, j]
            @views @. Rjp[:, j] += (Xh[:, ti, i] / γ) * (Sh[:, ti, i] / N[i]) * cij
        end
        for j in 1:6
            wj = adaptive ? Weights(@view W[:, j]) : Weights(fill(1.0, Np))  # R_j uses group j weight
            q = quantile(@view(Rjp[:, j]), wj, [0.5, 0.025, 0.975])
            Rj[:, j, t] = q
            qd = quantile(@view(Ih[:, ti, j]), wj, [0.5, 0.025, 0.975])  # cases CrI from the coherent path
            DR[:, j, t] = qd
        end
    end
    runtime = time() - t0
    return (Rj=Rj, DR=DR, ess=ess, uniq=uniq, rflag=rflag, runtime=runtime)
end

# --- metrics against the truth: MAPE, RMSE, bias, coverage, interval width ---
function metrics(est, truth; warmup=5)
    # est: (3,6,T) med/lo/hi ; truth: (6,T)
    T = size(truth, 2); rng = (warmup+1):T
    med = est[1, :, :]; lo = est[2, :, :]; hi = est[3, :, :]
    mape = zeros(6); rmse = zeros(6); bias = zeros(6); cov = zeros(6); wid = zeros(6)
    for i in 1:6
        tr = truth[i, rng]; m = med[i, rng]
        keep = tr .> 0
        mape[i] = mean(abs.(m[keep] .- tr[keep]) ./ tr[keep]) * 100
        rmse[i] = sqrt(mean((m .- tr) .^ 2))
        bias[i] = mean(m .- tr)
        cov[i]  = mean((tr .>= lo[i, rng]) .& (tr .<= hi[i, rng])) * 100
        wid[i]  = mean(hi[i, rng] .- lo[i, rng])
    end
    return (mape=mape, rmse=rmse, bias=bias, cov=cov, wid=wid,
            mape_all=mean(mape), rmse_all=mean(rmse), cov_all=mean(cov))
end


# --- output helpers ---
Llab(L) = L >= Tobs ? "wh" : string(L)
dn = joinpath(normpath(joinpath(@__DIR__,"..")), "results"); isdir(dn) || mkpath(dn)

function save_curves!(rows, Rj, DR, Np, Ll, srw, rs)
    for t in 1:Tobs, i in 1:6
        push!(rows, (Np=Np, L=Ll, sigma_rw=srw, resample=rs, day=t, age_group=AGE[i],
                     true_Rj=trueRj[i,t], Rj_med=Rj[1,i,t], Rj_lo=Rj[2,i,t], Rj_hi=Rj[3,i,t],
                     true_dR=trueDR[i,t], dR_med=DR[1,i,t], dR_lo=DR[2,i,t], dR_hi=DR[3,i,t]))
    end
end

# ═════════════════════════════════════════════════════════════════════
#  Full sweep (RT_PS_MODE=sweep) : main 36(every-step) + adaptive baseline
# ═════════════════════════════════════════════════════════════════════
if get(ENV, "RT_PS_MODE", "test") == "sweep"
    Nps  = [parse(Int, x) for x in split(get(ENV, "NPS", "10000,100000,1000000"), ",")]
    Ls   = [16, 24, 32, Tobs + 10]          # the last entry means whole-history resampling
    SRWs = [0.1, 0.15, 0.2, 0.25, 0.3]
    basethr = parse(Int, get(ENV, "BASETHR", "100000"))
    baseNp = [n for n in Nps if n >= basethr]  # baseline Np (1e5·1e6)
    bL = 16; bsrw = 0.15                       # adopted window and random-walk scale
    seed = 1
    # noisy mode (--noise-r > 0): NB(r_gen) observation noise is added to the true incidence and the
    #   observation model switches to NB(30) in log space. --noise-r 0 keeps the clean data.
    RGEN = parse(Int, get(ENV, "RT_NOISY_RGEN", "0"))
    if RGEN > 0
        Random.seed!(20260808 + RGEN)          # reproducible noise, one seed per r_gen
        Yuse = Matrix{Int}(undef, 6, Tobs)
        for t in 1:Tobs, i in 1:6
            Yuse[i, t] = rand(NegativeBinomial(RGEN, RGEN / (RGEN + max(trueDR[i, t], 1e-8))))
        end
        NBR = 30; ROB = true; sfx = "_noisy_r$(RGEN)"
        println(">>> NOISY sweep: r_gen=$RGEN, weight=NB($NBR), log-space robust, suffix=$sfx")
    else
        Yuse = Y                                # clean data
        nbr_c = parse(Int, get(ENV, "RT_NBR", "0"))
        if nbr_c > 0                            # clean data with an NB weight: the main configuration
            NBR = nbr_c; ROB = true; sfx = ""   # empty suffix: writes the main output files
            println(">>> CLEAN sweep: weight=NB($NBR), log-space robust (clean main)")
        else                                    # clean data with a Poisson weight: recovery check
            NBR = 0; ROB = false; sfx = ""
        end
    end

    met = DataFrame()                          # metrics table
    curves = NamedTuple[]                       # baseline fit curves
    essrows = NamedTuple[]                      # baseline ESS(t)
    degrows = NamedTuple[]                      # degeneracy uniq(t) per L
    adrows  = NamedTuple[]                      # adaptive monitor
    addmetrow!(cfg, mR, mD, r) = push!(met, (
        Np=cfg.Np, L=Llab(cfg.L), sigma_rw=cfg.srw, resample=cfg.rs,
        mape_Rj=mR.mape_all, rmse_Rj=mR.rmse_all, bias_Rj=mean(abs, mR.bias),
        cov_Rj=mR.cov_all, width_Rj=mean(mR.wid),
        mape_dR=mD.mape_all, cov_dR=mD.cov_all,
        ess_min=minimum(r.ess), ess_med=median(r.ess), uniq_min=minimum(r.uniq),
        runtime_s=r.runtime,
        mape_Rj_0_5=mR.mape[1], mape_Rj_6_11=mR.mape[2], mape_Rj_12_17=mR.mape[3],
        mape_Rj_18_44=mR.mape[4], mape_Rj_45_64=mR.mape[5], mape_Rj_65=mR.mape[6]); cols=:union)

    # incremental writes, so a crash does not lose the sweep
    flush_csvs() = begin
        CSV.write(joinpath(dn, "ps_metrics$(sfx).csv"), met)
        !isempty(curves)  && CSV.write(joinpath(dn, "ps_baseline_curves$(sfx).csv"), DataFrame(curves))
        !isempty(essrows) && CSV.write(joinpath(dn, "ps_ess$(sfx).csv"), DataFrame(essrows))
        !isempty(degrows) && CSV.write(joinpath(dn, "ps_degeneracy$(sfx).csv"), DataFrame(degrows))
        !isempty(adrows)  && CSV.write(joinpath(dn, "ps_adaptive$(sfx).csv"), DataFrame(adrows))
    end

    println("== MAIN SWEEP (every-step): Np×L×σ_rw = $(length(Nps))×$(length(Ls))×$(length(SRWs)) ==")
    for Np in Nps, L in Ls, srw in SRWs
        @printf("[main] Np=%d L=%s σ_rw=%.1f ... ", Np, Llab(L), srw); flush(stdout)
        r = run_ps(Np, L, srw; adaptive=false, seed=seed, Yobs=Yuse, nb_r=NBR, robust=ROB)
        mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR)
        addmetrow!((Np=Np, L=L, srw=srw, rs="every"), mR, mD, r)
        @printf("MAPE_Rj=%.2f%% cov=%.0f%% t=%.0fs\n", mR.mape_all, mR.cov_all, r.runtime)
        # baseline(L=16,σ_rw=0.3, Np∈base) : fit curves + ESS
        if L == bL && srw == bsrw && Np in baseNp
            save_curves!(curves, r.Rj, r.DR, Np, Llab(L), srw, "every")
            for t in 1:Tobs, i in 1:6
                push!(essrows, (Np=Np, day=t, age_group=AGE[i], ess=r.ess[i,t]))
            end
        end
        # degeneracy at the baseline Np for every window length
        if srw == bsrw && Np in baseNp
            for t in 1:Tobs, i in 1:6
                push!(degrows, (Np=Np, L=Llab(L), day=t, age_group=AGE[i], uniq=r.uniq[i,t]))
            end
        end
        flush_csvs()                       # write after every run so a crash loses nothing
    end

    println("== ADAPTIVE baseline (Np∈base, L=16, σ_rw=0.3, ESS<N/2) ==")
    for Np in baseNp
        @printf("[adapt] Np=%d ... ", Np); flush(stdout)
        r = run_ps(Np, bL, bsrw; adaptive=true, ess_frac=0.5, seed=seed, Yobs=Yuse, nb_r=NBR, robust=ROB)
        mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR)
        addmetrow!((Np=Np, L=bL, srw=bsrw, rs="adaptive"), mR, mD, r)
        @printf("MAPE_Rj=%.2f%% cov=%.0f%% resamples (mean over groups)=%.0f t=%.0fs\n",
                mR.mape_all, mR.cov_all, mean(sum(r.rflag, dims=2)), r.runtime)
        for t in 1:Tobs, i in 1:6
            push!(adrows, (Np=Np, day=t, age_group=AGE[i], ess=r.ess[i,t], resampled=r.rflag[i,t]))
        end
        flush_csvs()
    end
    println("\nDone. CSV → results/ps_{metrics,baseline_curves,ess,degeneracy,adaptive}.csv")
elseif get(ENV, "RT_PS_MODE", "test") == "test"
    Np = parse(Int, get(ENV, "NP", "10000"))
    L  = parse(Int, get(ENV, "L", "16"))
    σ_rw = parse(Float64, get(ENV, "SRW", "0.3"))
    adapt = get(ENV, "ADAPT", "0") == "1"
    println("== PS test: Np=$Np L=$L σ_rw=$σ_rw adaptive=$adapt every/adapt, Tobs=$Tobs ==")
    r = run_ps(Np, L, σ_rw; adaptive=adapt, seed=1)
    mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR)
    @printf("runtime = %.1f s\n", r.runtime)
    @printf("R_j  MAPE per group: %s | overall %.2f%%\n", join([@sprintf("%s=%.1f", AGE[i], mR.mape[i]) for i in 1:6], " "), mR.mape_all)
    @printf("R_j  95%% coverage per group: %s | overall %.0f%%\n", join([@sprintf("%s=%.0f", AGE[i], mR.cov[i]) for i in 1:6], " "), mR.cov_all)
    @printf("dR   MAPE overall %.2f%% | coverage %.0f%%\n", mD.mape_all, mD.cov_all)
    @printf("ESS(t) min/median over all groups: %.0f / %.0f\n", minimum(r.ess), median(r.ess))
    @printf("unique ancestors on day 1, per group: %s\n", join([string(r.uniq[i,1]) for i in 1:6], " "))
elseif get(ENV, "RT_PS_MODE", "test") == "sigma_extra"
    # extra sigma_rw points at L = 16, appended to ps_metrics.csv and de-duplicated
    seed = 1
    SRWs = [parse(Float64, x) for x in split(get(ENV, "SRWS", "0.15,0.25"), ",")]
    Nps  = [parse(Int, x) for x in split(get(ENV, "NPS", "100000,1000000"), ",")]
    mf = joinpath(dn, "ps_metrics.csv"); existing = CSV.read(mf, DataFrame)
    newrows = NamedTuple[]
    for Np in Nps, srw in SRWs
        @printf("[sigma_extra] Np=%d L=16 σ_rw=%.2f ... ", Np, srw); flush(stdout)
        r = run_ps(Np, 16, srw; adaptive=false, seed=seed)
        mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR)
        push!(newrows, (Np=Np, L="16", sigma_rw=srw, resample="every",
            mape_Rj=mR.mape_all, rmse_Rj=mR.rmse_all, bias_Rj=mean(abs, mR.bias),
            cov_Rj=mR.cov_all, width_Rj=mean(mR.wid), mape_dR=mD.mape_all, cov_dR=mD.cov_all,
            ess_min=minimum(r.ess), ess_med=median(r.ess), uniq_min=minimum(r.uniq), runtime_s=r.runtime,
            mape_Rj_0_5=mR.mape[1], mape_Rj_6_11=mR.mape[2], mape_Rj_12_17=mR.mape[3],
            mape_Rj_18_44=mR.mape[4], mape_Rj_45_64=mR.mape[5], mape_Rj_65=mR.mape[6]))
        @printf("MAPE=%.3f%% cov=%.1f%%\n", mR.mape_all, mR.cov_all)
    end
    combined = vcat(existing, DataFrame(newrows); cols=:union)
    unique!(combined, [:Np, :L, :sigma_rw, :resample])
    sort!(combined, [:Np, :L, :sigma_rw])
    CSV.write(mf, combined); println("ps_metrics.csv updated (extra sigma points appended, duplicates removed)")
elseif get(ENV, "RT_PS_MODE", "test") == "rebaseline"
    # recompute the baseline at a different sigma_rw: regenerate every baseline output and re-check
    #   L, Np and the resampling schedule. ps_metrics is appended and de-duplicated; the other files
    srw = parse(Float64, get(ENV, "SRW", "0.15")); seed = 1
    Ls = [16, 24, 32, Tobs + 10]; baseNp = [100000, 1000000]
    econfigs = vcat([(10000, 16)], [(Np, L) for Np in baseNp for L in Ls])  # 1e4 at L16 plus every L
    curves = NamedTuple[]; essrows = NamedTuple[]; degrows = NamedTuple[]; adrows = NamedTuple[]; newmet = NamedTuple[]
    addrow!(Np, L, rs, mR, mD, r) = push!(newmet, (Np=Np, L=Llab(L), sigma_rw=srw, resample=rs,
        mape_Rj=mR.mape_all, rmse_Rj=mR.rmse_all, bias_Rj=mean(abs, mR.bias),
        cov_Rj=mR.cov_all, width_Rj=mean(mR.wid), mape_dR=mD.mape_all, cov_dR=mD.cov_all,
        ess_min=minimum(r.ess), ess_med=median(r.ess), uniq_min=minimum(r.uniq), runtime_s=r.runtime,
        mape_Rj_0_5=mR.mape[1], mape_Rj_6_11=mR.mape[2], mape_Rj_12_17=mR.mape[3],
        mape_Rj_18_44=mR.mape[4], mape_Rj_45_64=mR.mape[5], mape_Rj_65=mR.mape[6]))
    println("== REBASELINE σ_rw=$srw : every $(length(econfigs)) + adaptive $(length(baseNp)) ==")
    for (Np, L) in econfigs
        @printf("[every] Np=%d L=%s ... ", Np, Llab(L)); flush(stdout)
        r = run_ps(Np, L, srw; adaptive=false, seed=seed)
        mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR); addrow!(Np, L, "every", mR, mD, r)
        @printf("MAPE=%.2f%% cov=%.0f%% t=%.0fs\n", mR.mape_all, mR.cov_all, r.runtime)
        if L == 16 && Np in baseNp
            save_curves!(curves, r.Rj, r.DR, Np, Llab(L), srw, "every")
            for t in 1:Tobs, i in 1:6; push!(essrows, (Np=Np, day=t, age_group=AGE[i], ess=r.ess[i, t])); end
        end
        if Np in baseNp
            for t in 1:Tobs, i in 1:6; push!(degrows, (Np=Np, L=Llab(L), day=t, age_group=AGE[i], uniq=r.uniq[i, t])); end
        end
    end
    for Np in baseNp
        @printf("[adapt] Np=%d ... ", Np); flush(stdout)
        r = run_ps(Np, 16, srw; adaptive=true, ess_frac=0.5, seed=seed)
        mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR); addrow!(Np, 16, "adaptive", mR, mD, r)
        @printf("MAPE=%.2f%% cov=%.0f%% resample=%.0f t=%.0fs\n", mR.mape_all, mR.cov_all, mean(sum(r.rflag, dims=2)), r.runtime)
        for t in 1:Tobs, i in 1:6; push!(adrows, (Np=Np, day=t, age_group=AGE[i], ess=r.ess[i, t], resampled=r.rflag[i, t])); end
    end
    mf = joinpath(dn, "ps_metrics.csv"); existing = CSV.read(mf, DataFrame)
    combined = vcat(existing, DataFrame(newmet); cols=:union)
    unique!(combined, [:Np, :L, :sigma_rw, :resample]); sort!(combined, [:Np, :L, :sigma_rw])
    CSV.write(mf, combined)
    CSV.write(joinpath(dn, "ps_baseline_curves.csv"), DataFrame(curves))
    CSV.write(joinpath(dn, "ps_ess.csv"), DataFrame(essrows))
    CSV.write(joinpath(dn, "ps_degeneracy.csv"), DataFrame(degrows))
    CSV.write(joinpath(dn, "ps_adaptive.csv"), DataFrame(adrows))
    println("REBASELINE done @σ=$srw : ps_metrics append(dedup) + curves/ess/degeneracy/adaptive overwrite")
elseif get(ENV, "RT_PS_MODE", "test") == "weight_sweep"
    # Observation-model sensitivity on noisy synthetic data, scored against the true R_j:
    #   Poisson vs NB(r) at the baseline Np = 1e5, L = 16, every step.
    seed = 1
    noisy = CSV.read(joinpath(dn, "pre_confirm_noisy.csv"), DataFrame)
    Ynoisy = round.(Int, reshape(noisy.incidence, 6, Tobs))
    wts = [("Poisson", 0), ("NB10", 10), ("NB20", 20), ("NB30", 30), ("NB40", 40)]
    SRWs = [0.1, 0.15, 0.2]
    rough(Rj) = mean(mean(abs.(diff(@view Rj[1, j, :]))) for j in 1:6)   # mean absolute day-to-day change
    wrows = NamedTuple[]; wcurves = NamedTuple[]
    println("== WEIGHT SWEEP (noisy synthetic data, r_gen=20): weight x sigma_rw = $(length(wts)) x $(length(SRWs)) ==")
    for (wname, nbr) in wts, srw in SRWs
        @printf("[weight] %-7s σ=%.2f ... ", wname, srw); flush(stdout)
        r = run_ps(100000, 16, srw; adaptive=false, seed=seed, Yobs=Ynoisy, nb_r=nbr, robust=true)
        mR = metrics(r.Rj, trueRj)
        push!(wrows, (weight=wname, nb_r=nbr, sigma_rw=srw, mape_Rj=mR.mape_all, cov_Rj=mR.cov_all,
            width_Rj=mean(mR.wid), roughness=rough(r.Rj),
            mape_Rj_0_5=mR.mape[1], mape_Rj_6_11=mR.mape[2], mape_Rj_12_17=mR.mape[3],
            mape_Rj_18_44=mR.mape[4], mape_Rj_45_64=mR.mape[5], mape_Rj_65=mR.mape[6]))
        @printf("MAPE=%.2f%% cov=%.0f%% width=%.3f rough=%.4f\n", mR.mape_all, mR.cov_all, mean(mR.wid), rough(r.Rj))
        if srw == 0.15 && (nbr == 0 || nbr == 20)       # representative curves: Poisson and NB20
            for t in 1:Tobs, j in 1:6
                push!(wcurves, (weight=wname, day=t, age_group=AGE[j],
                    true_Rj=trueRj[j, t], Rj_med=r.Rj[1, j, t], Rj_lo=r.Rj[2, j, t], Rj_hi=r.Rj[3, j, t],
                    true_dR=trueDR[j, t], dR_med=r.DR[1, j, t], dR_lo=r.DR[2, j, t], dR_hi=r.DR[3, j, t],
                    obs=Ynoisy[j, t]))
            end
        end
    end
    CSV.write(joinpath(dn, "ps_weight_sweep.csv"), DataFrame(wrows))
    CSV.write(joinpath(dn, "ps_weight_curves.csv"), DataFrame(wcurves))
    println("weight_sweep done → ps_weight_{sweep,curves}.csv")

elseif get(ENV, "RT_PS_MODE", "test") == "weight_misspec"
    # Non-circular justification of the NB dispersion: the true overdispersion r_gen is set OUTSIDE
    #   the estimation grid (15/25/50 against Poisson and NB 10/20/30/40), so the chosen r cannot be
    #   'the right answer by construction'. In the real data the true overdispersion is unknown.
    seed = 1
    Rgens = [15, 25, 50]                        # true overdispersion, outside the estimation grid
    wts = [("Poisson", 0), ("NB10", 10), ("NB20", 20), ("NB30", 30), ("NB40", 40)]
    srw = 0.15
    rough(Rj) = mean(mean(abs.(diff(@view Rj[1, j, :]))) for j in 1:6)
    mrows = NamedTuple[]
    println("== WEIGHT MISSPECIFICATION: true r_gen x estimated r = $(length(Rgens)) x $(length(wts)), sigma_rw=$srw ==")
    for rg in Rgens
        Random.seed!(20260808 + rg)             # one realisation per r_gen, shared by every r_est
        Ynoisy = Matrix{Int}(undef, 6, Tobs)
        for t in 1:Tobs, i in 1:6
            μ = max(trueDR[i, t], 1e-8)
            Ynoisy[i, t] = rand(NegativeBinomial(rg, rg / (rg + μ)))
        end
        for (wname, nbr) in wts
            @printf("[misspec] r_gen=%2d %-7s ... ", rg, wname); flush(stdout)
            r = run_ps(100000, 16, srw; adaptive=false, seed=seed, Yobs=Ynoisy, nb_r=nbr, robust=true)
            mR = metrics(r.Rj, trueRj)
            push!(mrows, (r_gen=rg, weight=wname, nb_r=nbr, mape_Rj=mR.mape_all,
                cov_Rj=mR.cov_all, width_Rj=mean(mR.wid), roughness=rough(r.Rj)))
            @printf("MAPE=%.2f%% cov=%.0f%% width=%.3f\n", mR.mape_all, mR.cov_all, mean(mR.wid))
        end
    end
    CSV.write(joinpath(dn, "ps_weight_misspec.csv"), DataFrame(mrows))
    println("weight_misspec done → ps_weight_misspec.csv")

elseif get(ENV, "RT_PS_MODE", "test") == "baseline_fit"
    # Baseline fit curves only, for one noise level: skips the full sweep and runs the adopted
    #   configuration once (Np = 1e5, L = 16, sigma_rw = 0.15, NB 30).
    RGEN = parse(Int, get(ENV, "RT_NOISY_RGEN", "15"))
    seed = 1
    if RGEN == 0                                    # clean data with NB30
        Yn = Y; tag = "clean_nb"
        @printf("== BASELINE FIT (clean, NB30, Np=1e5·L16·σ0.15) ==\n"); flush(stdout)
    else
        Random.seed!(20260808 + RGEN)
        Yn = Matrix{Int}(undef, 6, Tobs)
        for t in 1:Tobs, i in 1:6
            Yn[i, t] = rand(NegativeBinomial(RGEN, RGEN / (RGEN + max(trueDR[i, t], 1e-8))))
        end
        tag = "noisy_r$RGEN"
        @printf("== BASELINE FIT (noisy r_gen=%d, NB30, Np=1e5·L16·σ0.15) ==\n", RGEN); flush(stdout)
    end
    r = run_ps(100000, 16, 0.15; adaptive=false, seed=seed, Yobs=Yn, nb_r=30, robust=true)
    mR = metrics(r.Rj, trueRj); mD = metrics(r.DR, trueDR)
    @printf("MAPE_Rj=%.2f%% cov=%.0f%% MAPE_dR=%.2f%% t=%.0fs\n", mR.mape_all, mR.cov_all, mD.mape_all, r.runtime)
    curves = NamedTuple[]
    save_curves!(curves, r.Rj, r.DR, 100000, "16", 0.15, "every")
    # --noise-r 0 overwrites the main ps_baseline_curves.csv (clean, NB 30)
    outname = RGEN == 0 ? "ps_baseline_curves.csv" : "ps_baseline_curves_noisy_r$(RGEN).csv"
    CSV.write(joinpath(dn, outname), DataFrame(curves))
    println("baseline_fit done → $outname")

elseif get(ENV, "RT_PS_MODE", "test") == "save_noisy"
    # Write the noisy incidence to CSV so that code/03_renewal.jl reads exactly the same realisation.
    #   The generating code is identical to the sweep (seed 20260808 + r_gen, same loop order).
    for rg in [15, 25, 50]
        Random.seed!(20260808 + rg)
        Yn = Matrix{Int}(undef, 6, Tobs)
        for t in 1:Tobs, i in 1:6
            Yn[i, t] = rand(NegativeBinomial(rg, rg / (rg + max(trueDR[i, t], 1e-8))))
        end
        rows = NamedTuple[]
        for t in 1:Tobs, i in 1:6
            push!(rows, (day=t, age_group=AGE[i], incidence=Yn[i, t]))
        end
        CSV.write(joinpath(dn, "pre_confirm_noisy_r$(rg).csv"), DataFrame(rows))
        @printf("saved pre_confirm_noisy_r%d.csv (peak=%d)\n", rg, maximum(Yn))
    end
    println("save_noisy done")

elseif get(ENV, "RT_PS_MODE", "test") == "weight_curves_r25"
    # Observation-model comparison curves at r_gen = 25: Poisson vs NB(30), same realisation as above.
    #   Np = 1e5, L = 16, sigma_rw = 0.15.
    seed = 1
    noisy = CSV.read(joinpath(dn, "pre_confirm_noisy_r25.csv"), DataFrame)
    Ynoisy = round.(Int, reshape(noisy.incidence, 6, Tobs))
    wcurves = NamedTuple[]
    for (wname, nbr) in [("Poisson", 0), ("NB30", 30)]
        @printf("[wcurve25] %-7s ... ", wname); flush(stdout)
        r = run_ps(100000, 16, 0.15; adaptive=false, seed=seed, Yobs=Ynoisy, nb_r=nbr, robust=true)
        mR = metrics(r.Rj, trueRj)
        @printf("MAPE=%.2f%% cov=%.0f%%\n", mR.mape_all, mR.cov_all)
        for t in 1:Tobs, j in 1:6
            push!(wcurves, (weight=wname, day=t, age_group=AGE[j],
                true_Rj=trueRj[j, t], Rj_med=r.Rj[1, j, t], Rj_lo=r.Rj[2, j, t], Rj_hi=r.Rj[3, j, t],
                true_dR=trueDR[j, t], dR_med=r.DR[1, j, t], dR_lo=r.DR[2, j, t], dR_hi=r.DR[3, j, t],
                obs=Ynoisy[j, t]))
        end
    end
    CSV.write(joinpath(dn, "ps_weight_curves_r25.csv"), DataFrame(wcurves))
    println("weight_curves_r25 done → ps_weight_curves_r25.csv")

elseif get(ENV, "RT_PS_MODE", "test") == "pfps"
    # Particle filter (L = 0) vs particle smoother (L = 16) on clean synthetic data, where the true
    #   R_j is known. Structural differences are clearest without observation noise.
    #   Measures: cross-correlation lag, 95% coverage and interval score (Winkler), ESS and ancestors.
    #   (5) degeneracy: ESS(t)·unique-ancestor(t). PF=L0(filter)·PS=L16(fixed-lag smoother).
    seed = 1; Np = 100000; srw = 0.15; NBR = 30; warm = 5; α = 0.05
    # tau > 0 means the estimate lags the truth by tau days (argmax over the post-warmup window)
    function xlag(m, tr; maxlag = 14)
        best = 0; bestc = -Inf
        for τ in -maxlag:maxlag
            a = τ >= 0 ? tr[1:end-τ] : tr[1-τ:end]
            b = τ >= 0 ? m[1+τ:end]  : m[1:end+τ]
            length(a) < 15 && continue
            c = cor(a, b)
            if c > bestc; bestc = c; best = τ; end
        end
        return best
    end
    # Two regimes: clean, and noisy with r_gen = 25 (the realisation written by --mode save_noisy).
    #   Clean separates sharpness at equal coverage; noisy separates calibration.
    noisy = CSV.read(joinpath(dn, "pre_confirm_noisy_r25.csv"), DataFrame)
    Ynoisy = round.(Int, reshape(noisy.incidence, 6, Tobs))
    srows = NamedTuple[]; crows = NamedTuple[]; erows = NamedTuple[]
    for (reg, Yuse) in [("clean", Y), ("noisy", Ynoisy)], (mname, L) in [("PF", 0), ("PS", 16)]
        @printf("[pfps] %-5s %s (L=%d) ... ", reg, mname, L); flush(stdout)
        r = run_ps(Np, L, srw; adaptive = false, seed = seed, Yobs = Yuse, nb_r = NBR, robust = true)
        mR = metrics(r.Rj, trueRj)
        rng = (warm+1):Tobs
        for j in 1:6
            med = r.Rj[1, j, :]; lo = r.Rj[2, j, :]; hi = r.Rj[3, j, :]; tr = trueRj[j, :]
            lag = xlag(med[rng], tr[rng]; maxlag = 14)
            is = 0.0                                   # Winkler interval score: width plus a miscalibration penalty
            for t in rng
                y = tr[t]; l = lo[t]; h = hi[t]
                is += (h - l) + (y < l ? (2/α)*(l - y) : 0.0) + (y > h ? (2/α)*(y - h) : 0.0)
            end
            is /= length(rng)
            push!(srows, (regime = reg, method = mname, L = L, age_group = AGE[j], mape_Rj = mR.mape[j],
                cov_Rj = mR.cov[j], width_Rj = mR.wid[j], iscore = is, lag_days = lag))
        end
        @printf("MAPE=%.2f%% cov=%.0f%% ESS med=%.0f uniq med=%.0f t=%.0fs\n",
                mR.mape_all, mR.cov_all, median(r.ess), median(r.uniq), r.runtime)
        for t in 1:Tobs, j in 1:6
            push!(crows, (regime = reg, method = mname, day = t, age_group = AGE[j],
                true_Rj = trueRj[j, t], Rj_med = r.Rj[1, j, t], Rj_lo = r.Rj[2, j, t], Rj_hi = r.Rj[3, j, t]))
            push!(erows, (regime = reg, method = mname, day = t, age_group = AGE[j], ess = r.ess[j, t], uniq = r.uniq[j, t]))
        end
    end
    CSV.write(joinpath(dn, "ps_pfps_summary.csv"), DataFrame(srows))
    CSV.write(joinpath(dn, "ps_pfps_curves.csv"), DataFrame(crows))
    CSV.write(joinpath(dn, "ps_pfps_ess.csv"), DataFrame(erows))
    println("pfps done → ps_pfps_{summary,curves,ess}.csv")
end
