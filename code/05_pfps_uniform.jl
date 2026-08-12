# age-rt-smoother / code/05_pfps_uniform.jl — particle filter vs particle smoother (paper Fig S2)
# ─────────────────────────────────────────────────────────────────────────────
# Uniform-mixing (non-age-structured) SEIR with a prescribed R_t (step and sinusoidal), estimated
#   by a particle filter and by the fixed-window particle smoother. Rates are set to the values of
#   the main analysis, sigma = gamma = 1/1.8. Only the adopted configuration is run — NB(r = 30)
#   observation model, window L = 16, Np = 1e5, sigma_rw = 0.3, E(0) = I(0) = 10 — and reported as
#   median with 95% credible interval. Purpose: make the intrinsic lag of filtering visible against
#   a truth that changes abruptly, which the smooth age-structured truth does not reveal.
# Output: results/pfps_uniform_{curves,meta}.csv    Run: julia --project=. code/05_pfps_uniform.jl
# ─────────────────────────────────────────────────────────────────────────────
using DifferentialEquations, Distributions, Random, StatsBase, Statistics, Printf
using Base.Threads

function seir!(du, u, p, t)                     # ODE used to generate the synthetic epidemic
    β = p.β(t); S, E, I, R = u
    du[1] = -β * S * I / p.N
    du[2] =  β * S * I / p.N - p.σ * E
    du[3] =  p.σ * E - p.γ * I
    du[4] =  p.γ * I
    return nothing
end
function generate_fractions(β_func, σ, γ, Nref, u0, T; tstops=Float64[])
    p = (β=β_func, σ=σ, γ=γ, N=Nref)
    sol = solve(ODEProblem(seir!, u0, (0.0, Float64(T)), p), Tsit5();
                saveat=0:T, tstops=tstops, abstol=1e-9, reltol=1e-9)
    R_vec = [u[4] for u in sol.u]; S_vec = [u[1] for u in sol.u]
    return diff(R_vec) ./ Nref, S_vec[1:end-1] ./ Nref       # ΔR/N, S/N
end

function seir_batch!(du, u, p, t)               # batched propagation of Np particles
    σ, γ, N = p.σ, p.γ, p.N; betas = p.betas
    @inbounds @threads for i in axes(u, 1)
        S = u[i, 1]; E = u[i, 2]; I = u[i, 3]; β = betas[i]; sIN = S * I / N
        du[i, 1] = -β * sIN; du[i, 2] = β * sIN - σ * E
        du[i, 3] = σ * E - γ * I; du[i, 4] = γ * I
    end
    return nothing
end
function propagate_particles(S, E, I, betas, σ, γ, N)
    Np = length(S); u0 = Matrix{Float64}(undef, Np, 4)
    @views u0[:, 1] .= S; @views u0[:, 2] .= E; @views u0[:, 3] .= I; @views u0[:, 4] .= 0.0
    sol = solve(ODEProblem(seir_batch!, u0, (0.0, 1.0), (betas=betas, σ=σ, γ=γ, N=N)), Tsit5();
                save_everystep=false, abstol=1e-8, reltol=1e-8)
    u_end = sol.u[end]
    return u_end[:, 1], u_end[:, 2], u_end[:, 3], u_end[:, 4]
end

function importance_weight(weight, μ, y)        # observation model p(Y_t | 𝓡_t^(i))
    if weight === :poisson
        return pdf.(Poisson.(μ), y)
    else
        r = Float64(weight[2]); p = r ./ (r .+ max.(μ, 1e-10))
        return pdf.(NegativeBinomial.(r, p), y)
    end
end

# PF (:PF, win = t:t) and fixed-window PS (:PS_win, win = max(t-L,1):t) differ only in the resample slice.
function run_method(method::Symbol, weight, Y, σ, γ, N, Np, L, sigma_rw; seed=1, quantiles=true)
    Random.seed!(seed); leng = length(Y)
    S = zeros(Np, leng); E = zeros(Np, leng); I = zeros(Np, leng); R = zeros(Np, leng); beta = zeros(Np, leng)
    E[:, 1] = rand(Poisson(10), Np); I[:, 1] = rand(Poisson(10), Np)   # E0=I0=10
    S[:, 1] = N .- E[:, 1] .- I[:, 1]; R[:, 1] .= 0.0
    beta[:, 1] = 1.05 .* exp.(rand(Normal(0, sigma_rw), Np)) .* γ ./ (S[:, 1] ./ N)
    for t in 2:leng
        beta[:, t] = beta[:, t-1] .* exp.(rand(Normal(0, sigma_rw), Np))
        Sn, En, In, Rn = propagate_particles(S[:, t-1], E[:, t-1], I[:, t-1], beta[:, t], σ, γ, N)
        S[:, t] = Sn; E[:, t] = En; I[:, t] = In; R[:, t] = Rn
        w = importance_weight(weight, R[:, t], Y[t]); sw = sum(w)
        w = sw == 0.0 ? fill(1 / Np, Np) : w ./ sw
        inds = wsample(1:Np, w, Np, replace=true)
        win = method == :PF ? (t:t) : (max(t - L, 1):t)
        S[:, win] = S[inds, win]; E[:, t] = E[inds, t]; I[:, t] = I[inds, t]; R[:, t] = R[inds, t]
        beta[:, win] = beta[inds, win]
    end
    Rt_mean  = [mean(@views beta[:, t] ./ γ .* (S[:, t] ./ N)) for t in 1:leng]
    obs_mean = [mean(view(R, :, t)) for t in 1:leng]
    Rt_p = beta ./ γ .* (S ./ N); q(A, p) = [quantile(view(A, :, t), p) for t in 1:leng]
    return (Rt_mean=Rt_mean, Rt_med=q(Rt_p, 0.5), Rt_lo=q(Rt_p, 0.025), Rt_hi=q(Rt_p, 0.975),
            obs_mean=obs_mean, obs_med=q(R, 0.5), obs_lo=q(R, 0.025), obs_hi=q(R, 0.975))
end
mape(truth, est) = (m = truth .> 0; mean(abs.((truth[m] .- est[m]) ./ truth[m])) * 100)

# ── setup: sigma = gamma = 1/1.8 as in the main analysis (GI = Erlang(2, 1.8), mean 3.6 d) ──
σ = 1/1.8; γ = 1/1.8; T = 100; N = 1e7
u0 = [N * (1 - 2e-6), N * 1e-6, N * 1e-6, 0.0]           # E0=I0=10
sigma_rw = 0.3; Np = 100_000; L = 16; best_weight = (:nb, 30.0)

ratio_step = vcat(1.2 * ones(25), 1.8 * ones(25), 0.6 * ones(25), 1.2 * ones(25))
β_step = t -> (ratio_step[clamp(ceil(Int, t), 1, length(ratio_step))]) * γ
ratio_sin = 1.2 .+ 0.6 .* sin.(2π * (5/4) / 100 .* (1:100))
β_sin = t -> (ratio_sin[clamp(ceil(Int, t), 1, length(ratio_sin))]) * γ

println("generating synthetic (σ=γ=1/1.8) ...")
frac_step, s_step = generate_fractions(β_step, σ, γ, N, u0, T; tstops=collect(1.0:1.0:T-1))
frac_sin,  s_sin  = generate_fractions(β_sin,  σ, γ, N, u0, T)
Y_step = round.(Int, frac_step .* N); trueRt_step = ratio_step .* s_step
Y_sin  = round.(Int, frac_sin  .* N); trueRt_sin  = ratio_sin  .* s_sin
@printf("  peak Y: step=%d sin=%d\n", maximum(Y_step), maximum(Y_sin))

dn = joinpath(normpath(joinpath(@__DIR__,"..")), "results"); isdir(dn) || mkpath(dn)
open(joinpath(dn, "pfps_uniform_curves.csv"), "w") do io
    println(io, "scenario,day,quantity,series,mid,lo,hi")
    meta = IOBuffer(); println(meta, "scenario,Np,weight,L,ps_rt_mape,pf_rt_mape,ps_obs_mape,pf_obs_mape")
    for (scen, Y, trueR) in [("step", Y_step, trueRt_step), ("sin", Y_sin, trueRt_sin)]
        @printf("[%s] PS(L=%d) & PF, NB30, Np=%d ...\n", scen, L, Np)
        rPS = run_method(:PS_win, best_weight, Y, σ, γ, N, Np, L, sigma_rw; quantiles=true)
        rPF = run_method(:PF,     best_weight, Y, σ, γ, N, Np, 0, sigma_rw; quantiles=true)
        for t in 1:T
            @printf(io, "%s,%d,Rt,true,%.6f,%.6f,%.6f\n", scen, t, trueR[t], trueR[t], trueR[t])
            @printf(io, "%s,%d,Rt,PS,%.6f,%.6f,%.6f\n",   scen, t, rPS.Rt_med[t], rPS.Rt_lo[t], rPS.Rt_hi[t])
            @printf(io, "%s,%d,Rt,PF,%.6f,%.6f,%.6f\n",   scen, t, rPF.Rt_med[t], rPF.Rt_lo[t], rPF.Rt_hi[t])
            @printf(io, "%s,%d,obs,true,%d,%d,%d\n",       scen, t, Y[t], Y[t], Y[t])
            @printf(io, "%s,%d,obs,PS,%.6f,%.6f,%.6f\n",   scen, t, rPS.obs_med[t], rPS.obs_lo[t], rPS.obs_hi[t])
            @printf(io, "%s,%d,obs,PF,%.6f,%.6f,%.6f\n",   scen, t, rPF.obs_med[t], rPF.obs_lo[t], rPF.obs_hi[t])
        end
        ps_rt = mape(trueR[2:end], rPS.Rt_mean[2:end]); pf_rt = mape(trueR[2:end], rPF.Rt_mean[2:end])
        ps_ob = mape(Y[2:end], rPS.obs_mean[2:end]);    pf_ob = mape(Y[2:end], rPF.obs_mean[2:end])
        @printf("  Rt MAPE: PS=%.2f%% PF=%.2f%% | obs MAPE: PS=%.2f%% PF=%.2f%%\n", ps_rt, pf_rt, ps_ob, pf_ob)
        @printf(meta, "%s,%d,nb30,%d,%.4f,%.4f,%.4f,%.4f\n", scen, Np, L, ps_rt, pf_rt, ps_ob, pf_ob)
        GC.gc()
    end
    open(joinpath(dn, "pfps_uniform_meta.csv"), "w") do mio; write(mio, take!(meta)); end
end
println("done → results/pfps_uniform_{curves,meta}.csv")
