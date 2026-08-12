# age-rt-smoother / code/04_nhis.jl — STEP 4: NHIS, 11 seasons, plus sensitivity analyses
# ─────────────────────────────────────────────────────────────────────────────
# Estimates the age-specific instantaneous reproduction number for each influenza season in the
#   Korean NHIS claims data, with the age-structured SEIR particle smoother and, for comparison,
#   the renewal equation. Configuration: per-stratum resampling, NB(r = 30) likelihood evaluated
#   in log space, NNLS initial guess, window L = 16, sigma_rw = 0.15, Np = 1e5; the renewal
#   estimator uses the same Erlang generation interval truncated at K = 16.
# Preprocessing: pool provinces and sexes into a national daily series, fill missing days with 0,
#   apply a 7-day centred moving average, and cut seasons from week 36 to week 35 of the next year.
# Initial condition: E(0), I(0) ~ iid Poisson(mean of the 7 days preceding the season / gamma)
#   drawn per particle, R(0) = 0, S(0) = N - E(0) - I(0). Raw and moving-average runs use the same
#   rule, each on its own source. Since there is no ground truth here, the fit is judged by how
#   well the model incidence reproduces the observed cases.
# ─────────────────────────────────────────────────────────────────────────────
using CSV, DataFrames, Dates, Distributions, StatsBase, Statistics, LinearAlgebra
using DifferentialEquations, NonNegLeastSquares, Random, Printf

# ── CLI options (the corresponding environment variables also work) ──────────
#   with a value: --seasons 2018,2022 · --contact prem2017|prem2021 · --weight-r 30 (0 = Poisson)
#                 --window 16 (0 = particle filter) · --tag NAME (suffix of the output file)
#   flags       : --raw (no moving average) · --healthcare-seeking · --immunity · --immunity-hseek
let val = Dict("--seasons" => "NHIS_SEASONS", "--contact" => "RT_CM", "--weight-r" => "RT_NBR",
               "--window" => "RT_L", "--tag" => "RT_TAG"),
    flg = Dict("--raw" => "RT_NOSMOOTH", "--healthcare-seeking" => "RT_HSEEK",
               "--immunity" => "RT_IMMUNITY", "--immunity-hseek" => "RT_IMHS")
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if haskey(flg, a)
            ENV[flg[a]] = "1"; i += 1
        elseif haskey(val, a)
            i < length(ARGS) || error("missing value for $a")
            ENV[val[a]] = ARGS[i + 1]; i += 2
        else
            error("unknown option: $a")
        end
    end
end

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DDIR = joinpath(ROOT, "data")        # inputs: NHIS workbook export, population, contact matrices
const RDIR = joinpath(ROOT, "results")     # outputs: estimated series
isdir(RDIR) || mkpath(RDIR)
const AGE  = ["0-5", "6-11", "12-17", "18-44", "45-64", "65+"]
const σ = 1 / 1.8
const γ = 1 / 1.8
C = [1.76 0.21 0.03 0.20 0.05 0.05;                # contact matrix; replaced by --contact
     0.33 3.75 0.30 0.31 0.14 0.13;
     0.05 0.31 3.65 0.22 0.31 0.12;
     2.07 2.11 1.38 1.78 1.26 0.84;
     0.47 0.88 1.87 1.18 1.97 1.47;
     0.30 0.47 0.42 0.45 0.83 2.50]

# population N_i: Pops_Dec2022.csv aggregated to the six bands
let popdf = CSV.read(joinpath(DDIR, "Pops_Dec2022.csv"), DataFrame; header = 1)
    pops = [parse(Int, replace(String(popdf[1, c]), "," => "")) for c in 4:ncol(popdf)]
    lb = [1, 7, 13, 19, 46, 66]; ub = [6, 12, 18, 45, 65, length(pops)]
    global const N = Float64[sum(pops[lb[i]:ub[i]]) for i in 1:6]
end

# Erlang generation interval (sigma = gamma = 1/1.8), discretised as in code/03_renewal.jl
const Kgi = 16
const WGI = let Ga = Gamma(2.0, 1.8), F(x) = x <= 0 ? 0.0 : cdf(Ga, x)
    w = zeros(Kgi); w[1] = F(1.5) - F(0.5)
    for d in 2:Kgi; w[d] = F(d + 0.5) - F(d - 0.5); end
    w ./ sum(w)
end

# ── KDCA epidemiological year and week ──
function kdca_yw(dt::Date)
    cy = year(dt); fdcy = Date(cy, 1, 1)
    dowf = dayofweek(fdcy) % 7                       # Sun=0..Sat=6 (Python (weekday()+1)%7)
    kw = ((dt - fdcy).value + dowf) ÷ 7 + 1; ky = cy
    fdny = Date(cy + 1, 1, 1); down = dayofweek(fdny) % 7
    if (fdny - dt).value <= down; ky = cy + 1; kw = 1; end
    return ky, kw
end

# ── NHIS export → national daily counts (6 x dates) → missing days 0 → 7-day centred MA ──
const NHIS_CSV = joinpath(RDIR, "NHIS_Flu.csv")
isfile(NHIS_CSV) || error("""
    results/NHIS_Flu.csv not found.
    The NHIS claims data are not distributed with this repository. Download the workbook from
    https://www.data.go.kr/data/15089429/fileData.do , rename it to data/NHIS.xlsx and run
        julia --project=. code/00_prepare_nhis.jl
    which writes this file (see data/README.md).""")
println("loading NHIS_Flu.csv ...")
raw = CSV.read(NHIS_CSV, DataFrame)
agg_daily = combine(groupby(raw, [:date, :agg]), :N => sum => :N)   # pool provinces and sexes
agg_daily.date = Date.(agg_daily.date)
dmin, dmax = minimum(agg_daily.date), maximum(agg_daily.date)
alldates = collect(dmin:Day(1):dmax); nd = length(alldates)
didx = Dict(d => i for (i, d) in enumerate(alldates))
counts = zeros(Float64, 6, nd)                       # [age, date]
for r in eachrow(agg_daily); counts[r.agg, didx[r.date]] = r.N; end
# 7-day centred moving average (partial window clamped at the ends)
smooth = similar(counts)
for a in 1:6, t in 1:nd
    lo = max(1, t - 3); hi = min(nd, t + 3)
    smooth[a, t] = mean(@view counts[a, lo:hi])
end
kyw = [kdca_yw(d) for d in alldates]                 # (ky,kw) per date

# --raw keeps the unsmoothed daily counts instead of the moving average.
nosmooth = get(ENV, "RT_NOSMOOTH", "") == "1"

# ── season matrix: week 36 of start_year to week 35 of the next year ──
function season_data(start_year)
    idx = [t for t in 1:nd if (kyw[t][1] == start_year && kyw[t][2] >= 36) ||
                              (kyw[t][1] == start_year + 1 && kyw[t][2] <= 35)]
    sort!(idx, by = t -> alldates[t])
    Yf = (nosmooth ? counts : smooth)[:, idx]        # 6 x T, raw or 7-day MA cases
    return Yf, alldates[idx], [kyw[t][2] for t in idx]
end

# ── initial Poisson mean: mean cases over the 7 days preceding the season, divided by gamma ──
#   raw runs average the raw counts, moving-average runs average the smoothed values.
#   The window is causal (past only) and avoids the day-of-week effect of a single first day.
#   E(0), I(0) ~ iid Poisson(mEI); dividing by gamma matches the data scale because Y is
#   approximately gamma * I in quasi-steady state, and E = I follows from sigma = gamma.
function init_mean(dates1)
    src = nosmooth ? counts : smooth
    t1 = didx[dates1]; pw = max(1, t1 - 7):(t1 - 1)
    seed7 = vec(mean(view(src, :, pw), dims = 2))
    return max.(seed7 ./ γ, 1.0), seed7
end

# ═════════════════════════════════════════════════════════════════════
#  Particle smoother (age-structured SEIR, per-stratum resampling)
# ═════════════════════════════════════════════════════════════════════
function seir_rhs!(du, u, p, t)
    S = @view u[:, 1:6]; E = @view u[:, 7:12]; I = @view u[:, 13:18]
    mul!(p.CI, I, p.Ct); @. p.foi = p.χ * (S / p.Nr) * p.CI
    @views @. du[:, 1:6]   = -p.foi
    @views @. du[:, 7:12]  =  p.foi - σ * E
    @views @. du[:, 13:18] =  σ * E - γ * I
    @views @. du[:, 19:24] =  γ * I
    return nothing
end
step_day!(u0, χ, Ct, Nr, CI, foi) =
    solve(ODEProblem(seir_rhs!, u0, (0.0, 1.0), (χ=χ, Ct=Ct, Nr=Nr, CI=CI, foi=foi)),
          Tsit5(); dt=1.0, adaptive=false, save_everystep=false).u[end]

function init_chi(Np, σ_rw, S0; R0target=1.05)
    A = permutedims((C .* (S0 ./ N)) ./ γ)
    χ0 = Matrix{Float64}(undef, Np, 6)
    for p in 1:Np; χ0[p, :] = nonneg_lsq(A, R0target .* exp.(randn(6) .* σ_rw); alg=:nnls); end
    return χ0
end

# Y (6 x T integer counts), mEI = initial Poisson mean (length 6) → Rj and dR as 3 x 6 x T
#   (median, 2.5%, 97.5%). The initial state is a random variable: E(0), I(0) ~ iid Poisson(mEI)
#   per particle. The NNLS matrix uses the mean S(0); the per-particle spread is O(1e-4) relative.
function run_ps(Y, mEI; Np=100000, L=16, σ_rw=0.15, seed=1, nb_r=30, R0=zeros(6))
    Random.seed!(seed); T = size(Y, 2); whole = (L >= T)
    Xh = Array{Float64}(undef, Np, T + 1, 6); Sh = similar(Xh); Rh = similar(Xh)
    Ih = zeros(Float64, Np, T + 1, 6)                 # model incidence dR, stored to keep cases CrI coherent
    Ecur = Matrix{Float64}(undef, Np, 6); Icur = similar(Ecur)
    S0mean = N .- 2 .* mEI .- R0
    χ0 = init_chi(Np, σ_rw, S0mean); @views Xh[:, 1, :] .= χ0
    for i in 1:6
        @views Ecur[:, i] .= rand(Poisson(mEI[i]), Np)
        @views Icur[:, i] .= rand(Poisson(mEI[i]), Np)
        @views Sh[:, 1, i] .= N[i] .- Ecur[:, i] .- Icur[:, i] .- R0[i]
        @views Rh[:, 1, i] .= R0[i]
    end
    Ct = permutedims(C); Nr = reshape(N, 1, 6)
    CI = Matrix{Float64}(undef, Np, 6); foi = similar(CI); u = Matrix{Float64}(undef, Np, 24)
    for t in 1:T
        ti = t + 1; χt = @view Xh[:, ti, :]
        @inbounds for a in 1:6, p in 1:Np; χt[p, a] = Xh[p, ti-1, a] * exp(σ_rw * randn()); end
        @views u[:, 1:6] .= Sh[:, ti-1, :]; @views u[:, 7:12] .= Ecur
        @views u[:, 13:18] .= Icur; @views u[:, 19:24] .= Rh[:, ti-1, :]
        unext = step_day!(u, χt, Ct, Nr, CI, foi)
        @views Sh[:, ti, :] .= unext[:, 1:6]; @views Ecur .= unext[:, 7:12]
        @views Icur .= unext[:, 13:18]; @views Rh[:, ti, :] .= unext[:, 19:24]
        for i in 1:6                                  # per-stratum NB(r) weight (log-space robust) + resample
            dR = @views Rh[:, ti, i] .- Rh[:, ti-1, i]
            @views Ih[:, ti, i] .= dR                 # store aligned incidence along the trajectory
            yi = Y[i, t]                              # NB heavy tail avoids Poisson degeneracy at high counts
            lw = map(dR) do d
                isfinite(d) || return -Inf
                dd = max(d, 1e-8)
                nb_r > 0 ? logpdf(NegativeBinomial(nb_r, nb_r / (nb_r + dd)), yi) : logpdf(Poisson(dd), yi)
            end
            m = maximum(lw)
            w = isfinite(m) ? exp.(lw .- m) : fill(1.0, Np)   # log-space max subtraction guards underflow
            sw = sum(w); wn = sw > 0 ? w ./ sw : fill(1/Np, Np)
            inds = wsample(1:Np, wn, Np; replace=true)
            win = whole ? (1:ti) : (max(ti - L, 1):ti)
            @views Xh[:, win, i] .= Xh[inds, win, i]; @views Sh[:, win, i] .= Sh[inds, win, i]
            @views Rh[:, win, i] .= Rh[inds, win, i]
            @views Ih[:, win, i] .= Ih[inds, win, i]   # resample the incidence trajectory with the same indices
            @views Ecur[:, i] .= Ecur[inds, i]; @views Icur[:, i] .= Icur[inds, i]
        end
    end
    Rj = zeros(3, 6, T); DR = zeros(3, 6, T); RR = zeros(3, 6, T); SS = zeros(3, 6, T)
    Rjp = Matrix{Float64}(undef, Np, 6)
    for t in 1:T
        ti = t + 1; fill!(Rjp, 0.0)
        @inbounds for j in 1:6, i in 1:6
            @views @. Rjp[:, j] += (Xh[:, ti, i] / γ) * (Sh[:, ti, i] / N[i]) * C[i, j]
        end
        for j in 1:6
            Rj[:, j, t] = quantile(@view(Rjp[:, j]), [0.5, 0.025, 0.975])
            DR[:, j, t] = quantile(@view(Ih[:, ti, j]), [0.5, 0.025, 0.975])  # cases CrI from the coherent path
            RR[:, j, t] = quantile(@view(Rh[:, ti, j]), [0.5, 0.025, 0.975])   # recovered level R(t)
            SS[:, j, t] = quantile(@view(Sh[:, ti, j]), [0.5, 0.025, 0.975])   # susceptible level S(t)
        end
    end
    return Rj, DR, RR, SS
end

# ── renewal estimator (same convention as code/03_renewal.jl: past-only denominator) ──
function gorji(mat)                                  # mat: 6×T (smoothed cases)
    T = size(mat, 2); Rj = zeros(6, T)
    Λ(k, t) = (umax = t <= Kgi ? t - 1 : Kgi; s = 0.0;
               for u in 1:umax; s += mat[k, t-u] * WGI[u]; end; s)
    for t in 1:T
        D = [sum(C[i, k] * Λ(k, t) for k in 1:6) for i in 1:6]
        for j in 1:6
            Rj[j, t] = sum(D[i] == 0 ? NaN : (C[i, j] * mat[i, t]) / D[i] for i in 1:6)
        end
    end
    return Rj
end

# ═════════════════════════════════════════════════════════════════════
#  main: loop over seasons
# ═════════════════════════════════════════════════════════════════════
const ALL_SEASONS = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2022]
seasons = haskey(ENV, "NHIS_SEASONS") ?
    parse.(Int, split(ENV["NHIS_SEASONS"], ",")) : ALL_SEASONS

# Alternative contact matrix: --contact prem2017|prem2021 loads data/contact_<name>.csv.
#   In this mode only the R_j median is written (nhis_rj_<name>_<season>.csv).
altcm = get(ENV, "RT_CM", "")
if altcm != ""
    global C = Matrix{Float64}(CSV.read(joinpath(DDIR, "contact_$(altcm).csv"), DataFrame; header=false))
    @printf(">>> alternative contact matrix: %s (aggregated and transposed)\n", altcm)
end

# Healthcare-seeking correction: --healthcare-seeking multiplies cases by 100/p_seek per age group.
#   p_seek for influenza-like illness (Zhang 2020, Table 3): 91.7 (0-15), 84.0 (16-54), 100 (55+),
#   converted to the six bands with single-year population weights.
hseek = get(ENV, "RT_HSEEK", "") == "1"
const HS_FACTOR = [1.0905, 1.0905, 1.1206, 1.1905, 1.0886, 1.0000]   # 0-5·6-11·12-17·18-44·45-64·65+
if hseek; @printf(">>> healthcare-seeking correction: factor = %s\n", string(HS_FACTOR)); end

# Pre-existing immunity (one-year duration): --immunity sets R(0)_i to the previous season's
#   cumulative cases times the healthcare-seeking factor, and reduces S(0) accordingly.
#   Only the R_j median is written; R(0) is echoed to stdout.
immun = get(ENV, "RT_IMMUNITY", "") == "1"
if immun; println(">>> pre-existing immunity: R0 = HS_FACTOR * (previous season cases), S0 reduced"); end

# --immunity-hseek runs, for each season, the uncorrected baseline and the combined correction
#   and writes both to one file: recovered R(t), susceptible S(t), cases and R_j median.
imhs = get(ENV, "RT_IMHS", "") == "1"
if imhs; println(">>> immunity+healthcare-seeking overlay (baseline vs corrected)"); end

# Single-season method comparisons: --weight-r 0 (Poisson), --window 0 (particle filter), --tag NAME.
nbr_env = parse(Int, get(ENV, "RT_NBR", "30"))
L_env   = parse(Int, get(ENV, "RT_L", "16"))
ctag    = get(ENV, "RT_TAG", "")
if ctag != ""; @printf(">>> custom run: tag=%s nb_r=%d L=%d nosmooth=%s\n", ctag, nbr_env, L_env, string(nosmooth)); end

for sy in seasons
    if imhs                                           # ── baseline and corrected run, side by side ──
        Yf0, dates, wks = season_data(sy); T0 = size(Yf0, 2)
        mEIb, _ = init_mean(dates[1])                                          # baseline: R0 = 0, no correction
        @printf("[%d-%d IMHS] baseline … ", sy, sy+1); tb = time()
        Rjb, _, RRb, SSb = run_ps(round.(Int, Yf0), mEIb)
        Yf1 = HS_FACTOR .* Yf0                                                 # corrected: hseek cases + immunity R0
        mEIc = HS_FACTOR .* mEIb                                               # the seed is corrected too
        Yprev, = season_data(sy - 1); R0c = HS_FACTOR .* vec(sum(Yprev, dims = 2))
        R0c = min.(R0c, N .- 2 .* mEIc)
        Rjc, _, RRc, SSc = run_ps(round.(Int, Yf1), mEIc; R0 = R0c)
        @printf("corrected R0/age: %s | %.0fs\n",
                join([@sprintf("%s=%.0f", AGE[i], R0c[i]) for i in 1:6], " "), time() - tb)
        open(joinpath(RDIR, "nhis_overlay_imhs_$(sy).csv"), "w") do io
            println(io, "day,date,wk,age_group,cond,R_med,S_med,cases,Rj_med")
            for t in 1:T0, i in 1:6
                @printf(io, "%d,%s,%d,%s,baseline,%.4f,%.4f,%.4f,%.4f\n",
                        t, dates[t], wks[t], AGE[i], RRb[1, i, t], SSb[1, i, t], Yf0[i, t], Rjb[1, i, t])
                @printf(io, "%d,%s,%d,%s,corrected,%.4f,%.4f,%.4f,%.4f\n",
                        t, dates[t], wks[t], AGE[i], RRc[1, i, t], SSc[1, i, t], Yf1[i, t], Rjc[1, i, t])
            end
        end
        println("  → results/nhis_overlay_imhs_$(sy).csv")
        continue
    end
    Yf, dates, wks = season_data(sy)
    if hseek; Yf = HS_FACTOR .* Yf; end              # per-age scaling of observed cases
    T = size(Yf, 2)
    Y = round.(Int, Yf)                              # integer counts for the observation model
    mEI, seed7 = init_mean(dates[1])                 # mean of the preceding 7 days / gamma
    if hseek; mEI = HS_FACTOR .* mEI; seed7 = HS_FACTOR .* seed7; end   # same correction as the cases
    @printf(">>> [%d-%d] seed (mean cases over the preceding 7 days): %s | Poisson mean E0=I0 (=/gamma): %s\n", sy, sy+1,
            join([@sprintf("%s=%.1f", AGE[i], seed7[i]) for i in 1:6], " "),
            join([@sprintf("%s=%.1f", AGE[i], mEI[i]) for i in 1:6], " "))
    R0 = zeros(6)
    if immun                                          # R(0) = previous season cumulative cases x HS_FACTOR
        Yprev, = season_data(sy - 1)
        R0 = HS_FACTOR .* vec(sum(Yprev, dims = 2))
        R0 = min.(R0, N .- 2 .* mEI)                  # keep S(0) non-negative
        @printf(">>> [%d-%d] R(t=0)/age: %s | R0/N(%%): %s\n", sy, sy+1,
                join([@sprintf("%s=%.0f", AGE[i], R0[i]) for i in 1:6], " "),
                join([@sprintf("%s=%.2f", AGE[i], 100*R0[i]/N[i]) for i in 1:6], " "))
    end
    @printf("[%d-%d] T=%d days, observed peak per age: %s ... ", sy, sy+1, T,
            join([@sprintf("%s=%.0f", AGE[i], maximum(Yf[i, :])) for i in 1:6], " "))
    t0 = time()
    Rj, DR, RR, SS = run_ps(Y, mEI; R0 = R0, nb_r = nbr_env, L = L_env)
    @printf("done %.0fs\n", time() - t0)
    if altcm != "" || hseek || immun || nosmooth || ctag != ""   # compact output: cases and Rj median +- CrI
        tag = ctag != "" ? ctag : altcm != "" ? altcm : hseek ? "hseek" : immun ? "immunity" : "raw"
        open(joinpath(RDIR, "nhis_rj_$(tag)_$(sy).csv"), "w") do io
            println(io, "day,date,wk,age_group,cases,Rj_med,Rj_lo,Rj_hi")
            for t in 1:T, i in 1:6
                @printf(io, "%d,%s,%d,%s,%.4f,%.4f,%.4f,%.4f\n",
                        t, dates[t], wks[t], AGE[i], Yf[i, t], Rj[1, i, t], Rj[2, i, t], Rj[3, i, t])
            end
        end
        println("  → results/nhis_rj_$(tag)_$(sy).csv")
    else                                             # default: full output including the renewal estimate
        Rjr = gorji(Yf)
        open(joinpath(RDIR, "nhis_rt_$(sy).csv"), "w") do io
            println(io, "day,date,wk,age_group,obs,dR_med,dR_lo,dR_hi,Rj_med,Rj_lo,Rj_hi,Rj_renewal")
            for t in 1:T, i in 1:6
                @printf(io, "%d,%s,%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n",
                        t, dates[t], wks[t], AGE[i], Yf[i, t],
                        DR[1, i, t], DR[2, i, t], DR[3, i, t], Rj[1, i, t], Rj[2, i, t], Rj[3, i, t],
                        isnan(Rjr[i, t]) ? "NA" : @sprintf("%.4f", Rjr[i, t]))
            end
        end
        println("  → results/nhis_rt_$(sy).csv")
    end
end
