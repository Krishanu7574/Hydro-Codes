using CairoMakie, LinearAlgebra

# ====== PDE & helpers ======
const λ  = 7.0               # source strength
const Lx = 1.0                # domain length

@inline prev(i,N)  = (i==1 ? N : i-1)
@inline nexti(i,N) = (i==N ? 1 : i+1)

@inline function minmod(a::Float64, b::Float64)
    (a*b <= 0.0) && return 0.0
    return sign(a) * min(abs(a), abs(b))
end

"""
reconstruct_minmod(u, dx) -> uL, uR
Slope-limited linear reconstruction per cell.
uL[i] = left state at face i+1/2 (from cell i)
uR[i] = right state at face i-1/2 (from cell i)
"""
function reconstruct_minmod(u::Vector{Float64}, dx::Float64)
    N  = length(u)
    uL = similar(u)
    uR = similar(u)
    @inbounds for i in 1:N
        im, ip = prev(i,N), nexti(i,N)
        σ = minmod((u[i]-u[im])/dx, (u[ip]-u[i])/dx)
        uL[i] = u[i] + 0.5*σ*dx
        uR[i] = u[i] - 0.5*σ*dx
    end
    return uL, uR
end

# F[i] is flux at face i+1/2, for a=+1 use upwind = left state
function flux_faces_minmod!(F::Vector{Float64}, u::Vector{Float64}, dx::Float64)
    N = length(u)
    uL, _ = reconstruct_minmod(u, dx)
    @inbounds for i in 1:N
        F[i] = uL[i]              # f(u)=u
    end
    return F
end

"""
rhs_minmod!(du, u, dx, t)
Computes RHS L(u,t) = (F_{i-1/2} - F_{i+1/2})/dx + S(u,t),
with S(u,t) = -λ u and minmod-limited upwind flux for f.
"""
function rhs_minmod!(du::Vector{Float64}, u::Vector{Float64}, dx::Float64, t::Float64)
    N = length(u)
    F = similar(u)
    flux_faces_minmod!(F, u, dx)
    @inbounds for i in 1:N
        im = prev(i,N)
        divF  = (F[im] - F[i]) / dx
        du[i] = divF - λ*u[i]
    end
    return du
end

# ====== Analytical Solution ======
function analytical_solution(x::Vector{Float64}, t::Float64, λ::Float64)
    # For ∂u/∂t + ∂u/∂x = -λu with u(x,0) = sin(2πx)
    # Solution: u(x,t) = e^{-λt} * sin(2π(x - t))
    return exp(-λ * t) .* sin.(2π .* (x .- t))
end

# ====== SSPRK3 Method ======
"""
step_ssprk3_k123!(u, dt, dx, t)
Stage 1:
  k1 = f(y^n,t^n) + S(y^n,t^n)
  y^(1) = y^n + dt*k1
Stage 2:
  k2 = f(y^(1), t^n+dt) + S(y^(1), t^n+dt)
  y^(2) = 3/4 y^n + 1/4 y^(1) + 1/4 dt*k2
Stage 3:
  k3 = f(y^(2), t^n+dt/2) + S(y^(2), t^n+dt/2)
  y^{n+1} = 1/3 y^n + 2/3 y^(2) + 2/3 dt*k3
"""
function step_ssprk3_k123!(u::Vector{Float64}, dt::Float64, dx::Float64, t::Float64)
    k1 = similar(u); k2 = similar(u); k3 = similar(u)
    y1 = similar(u); y2 = similar(u)

    # Stage 1
    rhs_minmod!(k1, u, dx, t)              # k1 = f(y^n,t^n) + S(y^n,t^n)
    @. y1 = u + dt*k1                      # y^(1)

    # Stage 2
    rhs_minmod!(k2, y1, dx, t + dt)        # k2 at (y^(1), t^n+dt)
    @. y2 = (3/4)*u + (1/4)*y1 + (1/4)*dt*k2

    # Stage 3
    rhs_minmod!(k3, y2, dx, t + dt/2)      # k3 at (y^(2), t^n+dt/2)
    @. u  = (1/3)*u + (2/3)*y2 + (2/3)*dt*k3  # y^{n+1}

    return u
end

function evolve_ssprk3(u0::Vector{Float64}; T::Float64=0.125, CFL::Float64=0.8)
    N  = length(u0)
    dx = Lx/N
    dt = CFL*dx                  # a=+1 ⇒ dt = CFL*dx
    nsteps = Int(ceil(T/dt)); dt = T/nsteps
    t = 0.0
    u = copy(u0)
    for _ in 1:nsteps
        step_ssprk3_k123!(u, dt, dx, t)
        t += dt
    end
    return u, t, dt
end

# ====== Predictor-Corrector Method ======
function step_predictor_corrector_minmod!(u::Vector{Float64}, dt::Float64, dx::Float64)
    N  = length(u)
    F  = similar(u)
    Fh = similar(u)
    uh = similar(u)

    # predictor (n → n+1/2)
    flux_faces_minmod!(F, u, dx)
    @inbounds for i in 1:N
        im = prev(i,N)
        divF = (F[im] - F[i]) / dx
        uh[i] = u[i] + 0.5*dt*divF + 0.5*dt*(-λ*u[i])  # S^n
    end

    # corrector - trapezoidal rule for source
    flux_faces_minmod!(Fh, uh, dx)
    @inbounds for i in 1:N
        im = prev(i,N)
        divFh = (Fh[im] - Fh[i]) / dx
        u[i] = u[i] + dt*divFh + 0.5*dt*(-λ*u[i] - λ*uh[i])  # (S^n + S^{n+1/2})/2
    end
    return u
end

function evolve_predictor_corrector(u0::Vector{Float64}; T::Float64=0.125, CFL::Float64=0.8)
    N  = length(u0)
    dx = Lx/N
    dt = CFL*dx
    nsteps = Int(ceil(T/dt)); dt = T/nsteps
    t = 0.0
    u = copy(u0)
    for _ in 1:nsteps
        step_predictor_corrector_minmod!(u, dt, dx)
        t += dt
    end
    return u, t, dt
end

# ====== Main Comparison ======
function main()
    # Parameters
    N = 512
    x = collect(range(0, Lx, length=N+1)[1:end-1])  # Cell centers
    u0 = sin.(2π .* x)
    T_final = 0.1
    CFL = 0.4

    # Compute solutions
    println("Computing SSPRK3 solution...")
    u_ssprk3, t_final, dt_used = evolve_ssprk3(u0; T=T_final, CFL=CFL)
    
    println("Computing Predictor-Corrector solution...")
    u_pc, t_final_pc, dt_pc = evolve_predictor_corrector(u0; T=T_final, CFL=CFL)
    
    # Analytical solution
    u_analytical = analytical_solution(x, t_final, λ)

    # Compute errors
    error_ssprk3 = u_ssprk3 .- u_analytical
    error_pc = u_pc .- u_analytical

    # Error statistics
    println("\n" * "="^50)
    println("ERROR STATISTICS")
    println("="^50)
    println("SSPRK3 Method:")
    println("  Max error: $(maximum(abs.(error_ssprk3)))")
    println("  L2 error:  $(norm(error_ssprk3))")
    println("  RMS error: $(sqrt(sum(error_ssprk3.^2)/length(error_ssprk3)))")
    println("\nPredictor-Corrector Method:")
    println("  Max error: $(maximum(abs.(error_pc)))")
    println("  L2 error:  $(norm(error_pc))")
    println("  RMS error: $(sqrt(sum(error_pc.^2)/length(error_pc)))")
    println("\nSimulation Parameters:")
    println("  N: $N, CFL: $CFL, λ: $λ")
    println("  Time step: $(round(dt_used, sigdigits=6))")
    println("  Final time: $t_final")
    println("="^50)

    # Create plots
    fig = Figure(resolution=(1400, 1000), fontsize=16)

    # Plot 1: Solutions comparison
    ax1 = Axis(fig[1, 1], 
        title="Numerical Solutions Comparison (t = $(round(t_final, digits=3)), λ = $λ)",
        xlabel="x", 
        ylabel="u(x,t)",
        titlesize=20)

    lines!(ax1, x, u_analytical, linewidth = 2.0, label="Analytical", color=:cyan)
    scatter!(ax1, x[1:6:end], u_ssprk3[1:6:end], markersize = 10.0, marker = 'o', label="SSPRK3", color=:black)
    scatter!(ax1, x[1:8:end], u_pc[1:8:end], markersize = 10.0,marker = 'x' , label="Predictor-Corrector", color=:red)
    axislegend(ax1, position=:lb, framevisible=true)

    # Plot 2: Errors comparison
    ax2 = Axis(fig[1, 2], 
        title="Absolute Errors",
        xlabel="x", 
        ylabel="|Error|",
        titlesize=20)

    lines!(ax2, x, abs.(error_ssprk3), linewidth=2.5, label="SSPRK3 Error", color=:blue)
    lines!(ax2, x, abs.(error_pc), linewidth=2.5, label="PC Error", color=:red)
    axislegend(ax2, position=:rt, framevisible=true)

    # Plot 3: Pointwise difference between methods
    ax3 = Axis(fig[2, 1], 
        title="Difference: SSPRK3 - Predictor-Corrector",
        xlabel="x", 
        ylabel="SSPRK3 - PC",
        titlesize=20)

    diff_methods = u_ssprk3 .- u_pc
    lines!(ax3, x, diff_methods, linewidth=2.5, color=:purple)
    hlines!(ax3, [0.0], color=:gray, linestyle=:dash, linewidth=1)

    # Plot 4: Log-scale error magnitude
    ax4 = Axis(fig[2, 2], 
        title="Logarithmic Absolute Errors",
        xlabel="x", 
        ylabel="log10(|Error|)",
        titlesize=20)

    # Avoid log(0) by adding small constant
    epsilon = 1e-15
    lines!(ax4, x, log10.(abs.(error_ssprk3) .+ epsilon), linewidth=2.5, label="SSPRK3", color=:blue)
    lines!(ax4, x, log10.(abs.(error_pc) .+ epsilon), linewidth=2.5, label="PC", color=:red)
    axislegend(ax4, position=:rt, framevisible=true)

    # Add overall title
    Label(fig[0, :], "Advection-Reaction Equation: ∂u/∂t + ∂u/∂x = -λu with Minmod Slope Limiter", 
          fontsize =24, font=:bold)

    # Add statistics text box in a dedicated axis
    ax5 = Axis(fig[3, 1:2], 
        xlabel="", ylabel="",
        xgridvisible=false, ygridvisible=false,
        leftspinevisible=false, rightspinevisible=false,
        bottomspinevisible=false, topspinevisible=false,
        xticklabelsvisible=false, yticklabelsvisible=false,
        xticksvisible=false, yticksvisible=false)
    
    stats_text = "Error Statistics:\n" *
                 "SSPRK3 - L²: $(round(norm(error_ssprk3), sigdigits=5))\n" *
                 "PC     - L²: $(round(norm(error_pc), sigdigits=5))\n" *
                 "CFL: $CFL, N: $N\n" *
                 "λ: $λ, Δt: $(round(dt_used, sigdigits=4))"

    #text!(ax5, 0.5, 0.5, text=stats_text, align=(:center, :center),
    #      fontsize=14, font=:monospace, color=:darkgreen)
    
    text!(ax5, 0.5, 0.5, text=stats_text, align=(:center, :center),
      fontsize=14, font=:regular, color=:darkgreen)
      
      
    xlims!(ax5, 0, 1)
    ylims!(ax5, 0, 1)

    # Adjust layout
    rowsize!(fig.layout, 1, Relative(0.3))
    rowsize!(fig.layout, 2, Relative(0.3))
    rowsize!(fig.layout, 3, Relative(0.1))

    display(fig)
    
    # Save figure
    save("ssprk3_vs_predictor_corrector_comparison.png", fig, resolution=(1400, 1000))
    println("\nFigure saved as 'ssprk3_vs_predictor_corrector_comparison.png'")
    
    return u_ssprk3, u_pc, u_analytical, error_ssprk3, error_pc, x, dt_used
end

# ====== Convergence Test Function ======
function convergence_test()
    println("\n" * "="^60)
    println("CONVERGENCE TEST")
    println("="^60)
    
    T_final = 0.1
    CFL = 0.4
    N_values = [64, 128, 256, 512]
    
    errors_ssprk3 = Float64[]
    errors_pc = Float64[]
    dx_values = Float64[]
    
    for N in N_values
        println("Testing N = $N...")
        x = collect(range(0, Lx, length=N+1)[1:end-1])
        u0 = sin.(2π .* x)
        
        u_ssprk3, t_final, dt = evolve_ssprk3(u0; T=T_final, CFL=CFL)
        u_pc, _, _ = evolve_predictor_corrector(u0; T=T_final, CFL=CFL)
        
        u_analytical = analytical_solution(x, t_final, λ)
        
        error_ssprk3 = norm(u_ssprk3 .- u_analytical) / sqrt(N)
        error_pc = norm(u_pc .- u_analytical) / sqrt(N)
        
        push!(errors_ssprk3, error_ssprk3)
        push!(errors_pc, error_pc)
        push!(dx_values, Lx/N)
    end
    
    # Compute convergence rates
    println("\nConvergence Rates:")
    for i in 2:length(N_values)
        rate_ssprk3 = log(errors_ssprk3[i-1] / errors_ssprk3[i]) / log(2)
        rate_pc = log(errors_pc[i-1] / errors_pc[i]) / log(2)
        println("N $(N_values[i-1]) → $(N_values[i]): SSPRK3 = $(round(rate_ssprk3, digits=3)), PC = $(round(rate_pc, digits=3))")
    end
    
    return dx_values, errors_ssprk3, errors_pc
end

# ====== Run the complete analysis ======
println("Starting SSPRK3 vs Predictor-Corrector Comparison")
println("PDE: ∂u/∂t + ∂u/∂x = -λu with λ = $λ")
println("Initial condition: u(x,0) = sin(2πx)")
println("Domain: [0, $Lx] with periodic boundary conditions")

# Run main comparison
u_ssprk3, u_pc, u_analytical, error_ssprk3, error_pc, x, dt_used = main()

# Run convergence test
dx_vals, err_ssprk3, err_pc = convergence_test()

println("\n" * "="^60)
println("ANALYSIS COMPLETE")
println("="^60)
println("Summary:")
println("- Both methods use Minmod slope limiter for stability")
println("- SSPRK3: 3rd-order time integration, 3 stages per step") 
println("- Predictor-Corrector: 2nd-order time integration, 2 stages per step")
