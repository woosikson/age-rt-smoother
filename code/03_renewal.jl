# age-rt-smoother / code/03_renewal.jl — STEP 3: renewal-equation estimate vs particle smoother (paper Tables S6, S7)
# ─────────────────────────────────────────────────────────────────────────────
# Applies the age-structured renewal equation of Gorji et al. to the synthetic incidence written by
# code/01_synthetic.jl (results/pre_confirm.csv, daily new recoveries dR_i) and obtains a point
# estimate of the instantaneous R_j(t). No credible intervals: the renewal estimator as used here
# produces a point estimate only.
# The serial/generation interval is the Erlang implied by the SEIR rates, Gamma(shape 2, scale 1.8),
#   discretised by rounding, dropping day 0 and renormalising, truncated at K = 16 days.
# The whole series (day 1..T) is kept; no initial points are trimmed and no shift is applied.
#
# Estimator (beta = 1, M = contact matrix):
#   R_{ij}(t) = M[i,j]·mat[i,t] / Σ_k M[i,k]·Λ_k(t)      (i=recipient, j=infector)
#   Lambda_k(t) = sum_u mat[k, t-u] * w[u]
#   R_j(t) = sum_i R_ij(t)      (summing over recipients i gives the R of infector group j)
#
# Output: results/renewal_rt.csv          day, age_group, Rj_est, true_Rj
#         results/compare_metrics.csv     method x convention x age_group: mape, rmse, bias, lag
#                                         (the particle smoother is re-scored the same way)
# ─────────────────────────────────────────────────────────────────────────────

using CSV, DataFrames, Distributions, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DDIR = joinpath(ROOT, "results")     # this step only reads and writes generated files
const AGE  = ["0-5", "6-11", "12-17", "18-44", "45-64", "65+"]

# contact matrix C (same as code/01_synthetic.jl; Son et al. 2025) ────────────
const C = [1.76 0.21 0.03 0.20 0.05 0.05;
           0.33 3.75 0.30 0.31 0.14 0.13;
           0.05 0.31 3.65 0.22 0.31 0.12;
           2.07 2.11 1.38 1.78 1.26 0.84;
           0.47 0.88 1.87 1.18 1.97 1.47;
           0.30 0.47 0.42 0.45 0.83 2.50]

# ── 1. Erlang(σ=γ=1/1.8) discretized generation interval (gi_discretization.png ①) ──
#    Gamma(shape 2, scale 1.8); rounding, day 0 ([0,0.5)) removed, d = 1..K renormalised, K = 16.
const Kgi = 16
function erlang_gi()
    Ga = Gamma(2.0, 1.8)                       # mean 3.6 days
    F(x) = x <= 0 ? 0.0 : cdf(Ga, x)
    w = zeros(Kgi)
    w[1] = F(1.5) - F(0.5)
    for d in 2:Kgi; w[d] = F(d + 0.5) - F(d - 0.5); end
    w ./= sum(w)                              # redistribute the mass removed with day 0
    return w
end
const WGI = erlang_gi()                        # WGI[d] = probability that the GI equals d days
@printf("Erlang GI (Gamma(2,1.8)) discrete mean = %.3f days (K=%d)\n", sum((1:Kgi) .* WGI), Kgi)

# ── 2. inputs: pre_confirm.csv (incidence dR) → mat[6xT], pre_rt.csv → trueRj[6xT] ──
function load_matrix(csvfile, valcol)
    df = CSV.read(joinpath(DDIR, csvfile), DataFrame)
    days = sort(unique(df.day)); T = length(days)
    M = zeros(6, T)
    for r in eachrow(df)
        a = findfirst(==(r.age_group), AGE)
        t = findfirst(==(r.day), days)
        M[a, t] = r[valcol]
    end
    return M, T
end
mat, T   = load_matrix("pre_confirm.csv", :incidence)
trueRj,_ = load_matrix("pre_rt.csv", :Rt)
@printf("input: incidence mat 6x%d, true R_j 6x%d\n", T, T)

# ── 3. renewal estimate of R_j(t) (beta = 1) ──
# The denominator Lambda_k(t) uses past days only (day t excluded), consistently in both branches,
# and a partial window is left unnormalised, which is what produces the early spike.
#   t <= K : Lambda_k(t) = sum_{u=1}^{t-1} mat[k, t-u] * WGI[u]   (t-1 past days, partial weight)
#   t >  K : Lambda_k(t) = sum_{u=1}^{K}   mat[k, t-u] * WGI[u]   (K past days, full weight)
#   t =  1 : empty window, Lambda = 0, so the estimate is NaN and is kept as such
function gorji_lambda(M, k, t)
    umax = t <= Kgi ? (t - 1) : Kgi
    s = 0.0
    for u in 1:umax; s += M[k, t - u] * WGI[u]; end
    return s
end

function compute_renewal(M, T)                       # Gorji 𝓡_j(t) (M = incidence matrix 6×T)
    Rj = zeros(6, T)
    for t in 1:T
        Λ = [gorji_lambda(M, k, t) for k in 1:6]            # k=source group force-of-infection
        D = [sum(C[i, k] * Λ[k] for k in 1:6) for i in 1:6]  # D_i(t) = sum_k M[i,k] Lambda_k(t)
        for j in 1:6                                      # j = infector (column): R_j = sum_i R_ij
            Rj[j, t] = sum(D[i] == 0 ? NaN : (C[i, j] * M[i, t]) / D[i] for i in 1:6)
        end
    end
    return Rj
end
Rj_est = compute_renewal(mat, T)

open(joinpath(DDIR, "renewal_rt.csv"), "w") do io
    println(io, "day,age_group,Rj_est,true_Rj")
    for t in 1:T, j in 1:6
        @printf(io, "%d,%s,%.6f,%.6f\n", t, AGE[j], Rj_est[j, t], trueRj[j, t])
    end
end
println("wrote results/renewal_rt.csv")

# ── 4. particle-smoother baseline median, re-scored on the same convention ──
psdf = CSV.read(joinpath(DDIR, "ps_baseline_curves.csv"), DataFrame)
psb  = psdf[(psdf.Np .== 100000) .& (string.(psdf.L) .== "16") .& (psdf.sigma_rw .== 0.15) .& (psdf.resample .== "every"), :]
Rj_ps = zeros(6, T)
for r in eachrow(psb)
    a = findfirst(==(r.age_group), AGE); Rj_ps[a, r.day] = r.Rj_med
end

# ── 5. metrics: MAPE, RMSE, bias and lag under two conventions (all days / warmup 5) ──
# lag = argmin_{d in -20..20} sum_t (est[t] - truth[t-d])^2   (positive d = estimate lags truth)
function sse_shift(est, tr, δ)          # pairs containing NaN are skipped
    T = length(est); s = 0.0
    for t in 1:T
        (1 <= t - δ <= T) || continue
        (isnan(est[t]) || isnan(tr[t-δ])) && continue
        s += (est[t] - tr[t-δ])^2
    end
    return s
end
function best_lag(est, tr; maxlag=20)
    δs = -maxlag:maxlag
    return δs[argmin([sse_shift(est, tr, δ) for δ in δs])]
end
# sub-day resolution: fit a parabola to the three SSE values around the integer argmin
function best_lag_cont(est, tr; maxlag=20)
    δ0 = best_lag(est, tr; maxlag=maxlag)
    (δ0 <= -maxlag || δ0 >= maxlag) && return Float64(δ0)
    y1 = sse_shift(est, tr, δ0 - 1); y0 = sse_shift(est, tr, δ0); y2 = sse_shift(est, tr, δ0 + 1)
    d = y1 - 2y0 + y2
    return d == 0 ? Float64(δ0) : δ0 + 0.5 * (y1 - y2) / d
end
function metric_row(est6, warmup)
    rng = (warmup + 1):T
    mape = Float64[]; rmse = Float64[]; bias = Float64[]; lag = Int[]; lagc = Float64[]
    for j in 1:6
        e = est6[j, rng]; tr = trueRj[j, rng]
        keep = (tr .!= 0) .& .!isnan.(e)                 # truth non-zero and estimate defined
        push!(mape, mean(abs.(e[keep] .- tr[keep]) ./ tr[keep]) * 100)
        push!(rmse, sqrt(mean((e[keep] .- tr[keep]) .^ 2)))
        push!(bias, mean(e[keep] .- tr[keep]))
        push!(lag, best_lag(e, tr))
        push!(lagc, best_lag_cont(e, tr))
    end
    return mape, rmse, bias, lag, lagc
end

open(joinpath(DDIR, "compare_metrics.csv"), "w") do io
    println(io, "method,convention,age_group,mape,rmse,bias,lag,lag_cont")
    for (mname, est6) in (("renewal", Rj_est), ("PS", Rj_ps))
        for (cname, wu) in (("full", 0), ("warmup5", 5))
            mape, rmse, bias, lag, lagc = metric_row(est6, wu)
            for j in 1:6
                @printf(io, "%s,%s,%s,%.4f,%.5f,%.5f,%d,%.3f\n",
                        mname, cname, AGE[j], mape[j], rmse[j], bias[j], lag[j], lagc[j])
            end
            @printf(io, "%s,%s,all,%.4f,%.5f,%.5f,%.2f,%.3f\n",
                    mname, cname, mean(mape), mean(rmse), mean(bias), mean(lag), mean(lagc))
        end
    end
end
println("wrote results/compare_metrics.csv")

# ── 5b. initial estimation bias ──
# (i) magnitude: MAPE over the first W days   (ii) persistence: settling day
#   relerr_j(t) = |R̂_j(t) − R_j(t)| / R_j(t)
#   init-MAPE(W) = mean_{t=1..W} relerr_j(t) × 100
#   settle_eps = 1 + last{ t : relerr_j(t) >= eps }  (first day after which the error stays below eps)
relerr(est6, j) = abs.(est6[j, :] .- trueRj[j, :]) ./ trueRj[j, :]   # NaN is preserved
nanmean(v) = (w = filter(!isnan, v); isempty(w) ? NaN : mean(w))
init_mape(est6, j, W) = nanmean(relerr(est6, j)[1:W]) * 100
function settle_day(est6, j; ε=0.05)          # NaN counts as below the threshold
    re = relerr(est6, j); lb = findlast(x -> !isnan(x) && x >= ε, re)
    return lb === nothing ? 1 : lb + 1
end

open(joinpath(DDIR, "compare_initbias.csv"), "w") do io
    println(io, "method,age_group,mape_d1_5,mape_d1_10,mape_d1_16,mape_d6_end,settle5,settle10,maxrelerr_d1_16")
    for (mname, est6) in (("renewal", Rj_est), ("PS", Rj_ps))
        m15=Float64[]; m10=Float64[]; m16=Float64[]; m6=Float64[]; s5=Int[]; s10=Int[]; mx=Float64[]
        for j in 1:6
            push!(m15, init_mape(est6, j, 5));  push!(m10, init_mape(est6, j, 10))
            push!(m16, init_mape(est6, j, Kgi)); push!(m6, nanmean(relerr(est6, j)[6:T]) * 100)
            push!(s5, settle_day(est6, j; ε=0.05)); push!(s10, settle_day(est6, j; ε=0.10))
            push!(mx, maximum(filter(!isnan, relerr(est6, j)[1:Kgi])) * 100)
            @printf(io, "%s,%s,%.3f,%.3f,%.3f,%.3f,%d,%d,%.3f\n",
                    mname, AGE[j], m15[j], m10[j], m16[j], m6[j], s5[j], s10[j], mx[j])
        end
        @printf(io, "%s,all,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.3f\n",
                mname, mean(m15), mean(m10), mean(m16), mean(m6), mean(s5), mean(s10), mean(mx))
    end
end
println("wrote results/compare_initbias.csv")

println("\n============ initial estimation bias (early-window MAPE and settling day) ============")
cib = CSV.read(joinpath(DDIR, "compare_initbias.csv"), DataFrame)
@printf("%-8s %-7s | %8s %8s %8s %9s | %7s %7s %9s\n",
        "method","grp","MAPE1-5","1-10","1-16","MAPE6-end","settle5","settle10","max%1-16")
for mname in ("renewal", "PS")
    for r in eachrow(cib[cib.method .== mname, :])
        @printf("%-8s %-7s | %8.2f %8.2f %8.2f %9.2f | %7s %7s %9.1f\n",
                mname, r.age_group, r.mape_d1_5, r.mape_d1_10, r.mape_d1_16,
                r.mape_d6_end, string(r.settle5), string(r.settle10), r.maxrelerr_d1_16)
    end
end

# ── 6. console summary ──
println("\n================ renewal(Gorji) vs PS — MAPE·RMSE·bias·lag ================")
cdf = CSV.read(joinpath(DDIR, "compare_metrics.csv"), DataFrame)
for cname in ("full", "warmup5")
    println("\n[$cname]  (all = mean over the six age groups)")
    @printf("%-8s %-7s | %8s %9s %9s %6s %9s\n", "method", "grp", "MAPE%", "RMSE", "bias", "lag", "lag_cont")
    for mname in ("renewal", "PS")
        sub = cdf[(cdf.method .== mname) .& (cdf.convention .== cname), :]
        for r in eachrow(sub)
            @printf("%-8s %-7s | %8.3f %9.4f %9.4f %6s %9.2f\n",
                    mname, r.age_group, r.mape, r.rmse, r.bias, string(r.lag), r.lag_cont)
        end
    end
end

# ── 7. overdispersed synthetic data (r_gen in {15,25,50}): renewal vs particle smoother ──
#   The same realisation as the smoother (results/pre_confirm_noisy_r{rg}.csv) is used, so the
#   comparison is like for like. The renewal estimator has no observation model and therefore
#   follows the noise; the smoother with an NB(30) likelihood does not.
function load_ps_noisy(rg)
    f = joinpath(DDIR, "ps_baseline_curves_noisy_r$(rg).csv")
    isfile(f) || return nothing
    df = CSV.read(f, DataFrame); R = fill(NaN, 6, T)
    for r in eachrow(df); a = findfirst(==(r.age_group), AGE); R[a, r.day] = r.Rj_med; end
    return R
end
open(joinpath(DDIR, "compare_metrics_noisy.csv"), "w") do io
    println(io, "r_gen,method,convention,age_group,mape,rmse,bias,lag,lag_cont")
    for rg in [15, 25, 50]
        matn, _ = load_matrix("pre_confirm_noisy_r$(rg).csv", :incidence)
        Rj_rn = compute_renewal(matn, T)
        open(joinpath(DDIR, "renewal_rt_noisy_r$(rg).csv"), "w") do io2
            println(io2, "day,age_group,Rj_est,true_Rj")
            for t in 1:T, j in 1:6; @printf(io2, "%d,%s,%.6f,%.6f\n", t, AGE[j], Rj_rn[j, t], trueRj[j, t]); end
        end
        Rj_psn = load_ps_noisy(rg)
        ests = Rj_psn === nothing ? (("renewal", Rj_rn),) : (("renewal", Rj_rn), ("PS", Rj_psn))
        for (mname, est6) in ests, (cname, wu) in (("full", 0), ("warmup5", 5))
            mape, rmse, bias, lag, lagc = metric_row(est6, wu)
            for j in 1:6
                @printf(io, "%d,%s,%s,%s,%.4f,%.5f,%.5f,%d,%.3f\n",
                        rg, mname, cname, AGE[j], mape[j], rmse[j], bias[j], lag[j], lagc[j])
            end
            @printf(io, "%d,%s,%s,all,%.4f,%.5f,%.5f,%.2f,%.3f\n",
                    rg, mname, cname, mean(mape), mean(rmse), mean(bias), mean(lag), mean(lagc))
        end
    end
end
println("wrote results/{renewal_rt_noisy_r*, compare_metrics_noisy}.csv")

cdn = CSV.read(joinpath(DDIR, "compare_metrics_noisy.csv"), DataFrame)
println("\n========= noisy: renewal vs particle smoother (NB 30), warmup 5, mean over age groups =========")
@printf("%-6s %-8s | %8s %9s %9s %9s\n", "r_gen", "method", "MAPE%", "RMSE", "bias", "lag_cont")
for rg in [15, 25, 50], mname in ("renewal", "PS")
    sub = cdn[(cdn.r_gen .== rg) .& (cdn.method .== mname) .& (cdn.convention .== "warmup5") .& (cdn.age_group .== "all"), :]
    isempty(sub) && continue
    r = sub[1, :]
    @printf("%-6d %-8s | %8.3f %9.4f %9.4f %9.2f\n", rg, mname, r.mape, r.rmse, r.bias, r.lag_cont)
end
