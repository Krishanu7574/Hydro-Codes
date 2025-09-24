using CairoMakie 
using Roots 
using LinearAlgebra

#### Constants required for simulations #####

const N_L = 800
const Γ = 1.666
const xL = -0.5
const xR = 0.5
const small = 1e-12
const big = 1e12
const CFL = 0.35
const tf = 0.35

x = LinRange(xL,xR,N_L)
dx = (xR - xL) / N_L


function Lorentz(vx)
    v2 = vx*vx
    if v2 > 1 - small
        @warn "v^2 is close to 1.0; clamping"
        v2 = clamp(v2, small, 1 - small)
    end
    return clamp(sqrt(1/(1 - v2)), 1.0, big)
end


function prim_to_cons(u)
    n,vx,pg = u 
    γ = Lorentz(vx)
    e = pg/(Γ - 1) + n
    N_cons = n*γ
    Mx = (γ^2)*(e + pg) * vx
    E = (γ^2)*(e + pg) - pg
    if E < N_cons || E < Mx
        println("Physical constrains violated")
    end
    return [N_cons, Mx, E]
end

function prim_to_flux(u)
    n,vx,pg = u
    γ = Lorentz(vx)
    e = pg/(Γ - 1) + n
    N_cons = n*γ
    Mx = (γ^2)*(e + pg) * vx
    E = (γ^2)*(e + pg) - pg
    return [vx*N_cons, Mx*vx + pg , E*vx + pg*vx]
end

function cons_to_prim(U)
    N_cons , Mx , E = U
    function g(vx) 
        e_val = E - vx * Mx 
        n_val = N_cons * sqrt(1 - vx*vx) 
        P = (Γ - 1)*(e_val - n_val)
        return vx * ( E + P) - Mx
    end

    v_root = find_zero(g,(-1 + small, 1 - small), atol = small)
    if g(v_root) > small 
        println("roots are not in within tolerance ", [v_root, g(v_root)])
    end
    nr = N_cons/(Lorentz(v_root))
    er = E - v_root*Mx
    pg = (Γ - 1)*(er - nr)
    return [nr,v_root,pg]
end

function sound_speed(n,pg)
    e = pg/(Γ - 1) + n
    cs2 = (Γ*(Γ - 1)*(e - n))/ (n + Γ*(e - n))
    if cs2 > 1.0 - small || cs2 < small
        clamp(cs2,small, 1 - small)
    end
    return [cs2, sqrt(cs2)]
end 

function HLLE_bound_speeds(uL,uR)
    nL,vxL,pgL = uL 
    nR,vxR,pgR = uR 

    cs2L , csL = sound_speed(nL,pgL)
    cs2R , csR = sound_speed(nR,pgR)

    v_bar = 0.5*(vxL + vxR)
    cs_bar = 0.5*(csL + csR)

    λL = min(0.0, (v_bar - cs_bar) / (1 - v_bar*cs_bar), (vxL - csL) / ( 1 - vxL*csL) )
    λR = max(0.0, (v_bar + cs_bar) / (1 + v_bar*cs_bar), (vxR + csR) / ( 1 + vxR*csR) )

    return [λL , λR] 
end

function HLLE_flux(uL,uR)
    SL , SR = HLLE_bound_speeds(uL,uR)
    FL = prim_to_flux(uL) 
    FR = prim_to_flux(uR)
    UR = prim_to_cons(uR)
    UL = prim_to_cons(uL)

    if SL >= 0 
        return FL
    elseif SR <= 0
        return FR
    else 
        return (SR .* FL - SL .* FR + SL .* SR .* (UR .- UL)) ./ (SR .- SL)
    end
end

function maximum_signal(u)
    S_max = 0.0
    for i in 1:length(u)
        ni , vxi , pgi = u[i]
        cs2i , csi = sound_speed(ni,pgi)
        SL  = abs((vxi - csi) / (1 - vxi*csi))
        SR  = abs((vxi + csi) / (1 + vxi*csi))
        S_max = max(S_max, SL, SR)
    end
    return S_max
end


function minmod_slope_limiter_c(u)
    NC = length(u)
    uL_prim = deepcopy(u)
    uR_prim = deepcopy(u)
    #u_mid = deepcopy(u)

    for i in 1 : 3 
        for j in 2:NC - 1
            sp = u[j+1][i] - u[j][i]
            sm = u[j][i] - u[j-1][i]

            asp = abs(sp)
            asm = abs(sm)

            if sp*sm > 0.0
                if asp > asm
                    dU = 0.5*sm 
                else
                    dU = 0.5*sp
                end
            else
                dU = 0.0
            end

            uL_prim[j][i] = u[j][i] - dU
            uR_prim[j][i] = u[j][i] + dU
        end
    end
    return uL_prim , uR_prim
end


# Initial condition (Sod shock tube) - FIXED
function initial_condition(x)
    u_prim = Vector{Vector{Float64}}(undef, length(x))
    for i in 1:length(x)
        if x[i] < 0
            u_prim[i] = [10.0,0.0,13*(1/3)] #[1.0, 0.0, 1.0]   #       # # Left state
        else
            u_prim[i] = [1.0,0.0,(2/3)*(10^-6)] #[0.125, 0.0, 0.1]  #   # Right state
        end
    end
    return u_prim
end

#[10.0,0.0,13*(1/3)] : (1.0,0.0,(2/3)*(10^-6))

# Initialize simulation
u_prim = initial_condition(x)

#println([u_prim[300][1],u_prim[300][2],u_prim[300][3] ])

U = [prim_to_cons(u) for u in u_prim] 
U_half_L = deepcopy(U)
U_half_R = deepcopy(U)
U_new = deepcopy(U)



############## SSPRK3 time stepping (3rd-order Strong Stability Preserving Runge-Kutta) ###################

function SSPRK3_step(U, dx, dt, u_prim)
    N = length(U)
    
    # Stage 1
    U1 = deepcopy(U)
    u_prim1 = deepcopy(u_prim)
    
    # Apply slope limiter and compute fluxes
    uL_prim, uR_prim = minmod_slope_limiter_c(u_prim1)
    
    # Compute fluxes at cell interfaces
    F = Vector{Vector{Float64}}(undef, N-1)
    for i in 1:N-1
        F[i] = HLLE_flux(uR_prim[i], uL_prim[i+1])
    end
    
    # Update conservative variables (Stage 1)
    for i in 2:N-1
        U1[i] = U[i] - (dt/dx) .* (F[i] - F[i-1])
    end
    
    # Convert to primitive variables for Stage 2
    u_prim1 = [cons_to_prim(U1[i]) for i in 1:N]
    
    # Stage 2
    U2 = deepcopy(U1)
    uL_prim, uR_prim = minmod_slope_limiter_c(u_prim1)
    
    # Compute fluxes for Stage 2
    F2 = Vector{Vector{Float64}}(undef, N-1)
    for i in 1:N-1
        F2[i] = HLLE_flux(uR_prim[i], uL_prim[i+1])
    end
    
    # Update conservative variables (Stage 2)
    for i in 2:N-1
        U2[i] = 0.75 .* U[i] + 0.25 .* (U1[i] - (dt/dx) .* (F2[i] - F2[i-1]))
    end
    
    # Convert to primitive variables for Stage 3
    u_prim2 = [cons_to_prim(U2[i]) for i in 1:N]
    
    # Stage 3
    U3 = deepcopy(U2)
    uL_prim, uR_prim = minmod_slope_limiter_c(u_prim2)
    
    # Compute fluxes for Stage 3
    F3 = Vector{Vector{Float64}}(undef, N-1)
    for i in 1:N-1
        F3[i] = HLLE_flux(uR_prim[i], uL_prim[i+1])
    end
    
    # Final update (Stage 3)
    for i in 2:N-1
        U3[i] = (1/3) .* U[i] + (2/3) .* (U2[i] - (dt/dx) .* (F3[i] - F3[i-1]))
    end
    
    return U3
end

# Boundary conditions (transmissive/outflow)
function apply_boundary_conditions!(U)
    # Left boundary
    U[1] = U[2]
    # Right boundary  
    U[end] = U[end-1]
end

# Main simulation loop
function run_simulation()
    # Initialize
    u_prim = initial_condition(x)
    U = [prim_to_cons(u) for u in u_prim]
    
    t = 0.0
    iteration = 0
    
    # Create figure for real-time plotting
    fig = Figure(resolution = (1200, 800))
    
    # Create subplots for each primitive variable
    ax1 = Axis(fig[1, 1], title = "Density (n)", xlabel = "x", ylabel = "n")
    ax2 = Axis(fig[1, 2], title = "Velocity (vx)", xlabel = "x", ylabel = "vx")
    ax3 = Axis(fig[2, 1], title = "Pressure (pg)", xlabel = "x", ylabel = "pg")
    ax4 = Axis(fig[2, 2], title = "Lorentz Factor (γ)", xlabel = "x", ylabel = "γ")
    
    # Initial plots
    n_vals = [u[1] for u in u_prim]
    vx_vals = [u[2] for u in u_prim]
    pg_vals = [u[3] for u in u_prim]
    gamma_vals = [Lorentz(u[2]) for u in u_prim]
    
    lines!(ax1, x, n_vals, color = :blue, linewidth = 2, label = "t = 0")
    lines!(ax2, x, vx_vals, color = :red, linewidth = 2, label = "t = 0")
    lines!(ax3, x, pg_vals, color = :green, linewidth = 2, label = "t = 0")
    lines!(ax4, x, gamma_vals, color = :purple, linewidth = 2, label = "t = 0")
    
    axislegend(ax1)
    axislegend(ax2) 
    axislegend(ax3)
    axislegend(ax4)
    
    display(fig)
    
    # Simulation loop
    while t < tf
        # Calculate maximum signal speed for CFL condition
        S_max = maximum_signal(u_prim)
        dt = CFL * dx / S_max
        
        # Ensure we don't overshoot final time
        if t + dt > tf
            dt = tf - t
        end
        
        # Time step
        U_new = SSPRK3_step(U, dx, dt, u_prim)
        
        # Apply boundary conditions
        apply_boundary_conditions!(U_new)
        
        # Update variables
        U = deepcopy(U_new)
        u_prim = [cons_to_prim(U[i]) for i in 1:length(U)]
        
        # Update time and iteration counter
        t += dt
        iteration += 1
        
        # Print progress
        if iteration % 10 == 0
            println("Iteration: $iteration, Time: $t, dt: $dt")
        end
        
        # Update plot every 50 iterations
        if iteration % 50 == 0
            n_vals = [u[1] for u in u_prim]
            vx_vals = [u[2] for u in u_prim]
            pg_vals = [u[3] for u in u_prim]
            gamma_vals = [Lorentz(u[2]) for u in u_prim]
            
            # Clear previous plots and add new ones
            empty!(ax1)
            empty!(ax2)
            empty!(ax3)
            empty!(ax4)
            
            lines!(ax1, x, n_vals, color = :blue, linewidth = 2, label = "t = $(round(t, digits=3))")
            lines!(ax2, x, vx_vals, color = :red, linewidth = 2, label = "t = $(round(t, digits=3))")
            lines!(ax3, x, pg_vals, color = :green, linewidth = 2, label = "t = $(round(t, digits=3))")
            lines!(ax4, x, gamma_vals, color = :purple, linewidth = 2, label = "t = $(round(t, digits=3))")
            
            axislegend(ax1)
            axislegend(ax2)
            axislegend(ax3)
            axislegend(ax4)
            
            # Force display update
            display(fig)
        end
    end
    
    println("Simulation completed! Final time: $t, Iterations: $iteration")
    
    # Final plot
    fig_final = Figure(resolution = (1000, 600))
    
    axf1 = Axis(fig_final[1, 1], title = "Final Density Profile", xlabel = "x", ylabel = "n")
    axf2 = Axis(fig_final[1, 2], title = "Final Velocity Profile", xlabel = "x", ylabel = "vx")
    axf3 = Axis(fig_final[2, 1], title = "Final Pressure Profile", xlabel = "x", ylabel = "pg")
    axf4 = Axis(fig_final[2, 2], title = "Final Lorentz Factor", xlabel = "x", ylabel = "γ")
    
    n_vals = [u[1] for u in u_prim]
    vx_vals = [u[2] for u in u_prim]
    pg_vals = [u[3] for u in u_prim]
    gamma_vals = [Lorentz(u[2]) for u in u_prim]
    
    lines!(axf1, x, n_vals, color = :blue, linewidth = 2)
    lines!(axf2, x, vx_vals, color = :red, linewidth = 2)
    lines!(axf3, x, pg_vals, color = :green, linewidth = 2)
    lines!(axf4, x, gamma_vals, color = :purple, linewidth = 2)
    
    display(fig_final)
    save("Hydro_RHLLE_ONED_prim_SSPRK3.png", fig_final)
    return u_prim, U
end

# Run the simulation
final_prim, final_cons = run_simulation()

# Additional analysis functions
function calculate_conservation_errors(U_initial, U_final, dx)
    # Calculate total mass, momentum, energy conservation errors
    total_mass_initial = sum([U[1] for U in U_initial]) * dx
    total_momentum_initial = sum([U[2] for U in U_initial]) * dx
    total_energy_initial = sum([U[3] for U in U_initial]) * dx
    
    total_mass_final = sum([U[1] for U in U_final]) * dx
    total_momentum_final = sum([U[2] for U in U_final]) * dx
    total_energy_final = sum([U[3] for U in U_final]) * dx
    
    mass_error = abs(total_mass_final - total_mass_initial) / abs(total_mass_initial)
    momentum_error = abs(total_momentum_final - total_momentum_initial) / abs(total_momentum_initial)
    energy_error = abs(total_energy_final - total_energy_initial) / abs(total_energy_initial)
    
    println("Conservation Errors:")
    println("Mass: $(mass_error * 100)%")
    println("Momentum: $(momentum_error * 100)%")
    println("Energy: $(energy_error * 100)%")
    
    return mass_error, momentum_error, energy_error
end

# Calculate conservation errors
U_initial = [prim_to_cons(u) for u in initial_condition(x)]
mass_err, mom_err, energy_err = calculate_conservation_errors(U_initial, final_cons, dx)
