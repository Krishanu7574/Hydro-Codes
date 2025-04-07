using Plots

# Parameters

#γ = 4/3  # Adiabatic index for relativistic gas 5/3 , 4/3 , 1.999


println("Enter the value of γ : ")
γ = parse(Float64, readline())

if γ > 2 || γ <= 1
    println("Unphysical Adiabatic Constant")
    exit()  
end


c = 1.0  # Speed of light 
N = 500  # Number of grid points
L = 1.0  # Domain length
L_mid = L / 2

N_mid = N / 2
dx = L / N  # Grid spacing
dt = 0.0001  # Time step
T = 0.375  # Total simulation time
CFL = 0.001;  # CFL condition
println(L_mid)
println(N_mid)



# Initial condition 
function initial_condition(x)
    if x < L_mid
        ρ =  1.0  # Rest-mass density on the left  1.0, 10.0,2.0
        v =  0.0  # Velocity on the left
        p =  1.0  #13*(1/3)  # Pressure on the left   1.0, 40.0/3 ,1990
    else
        ρ =  0.125  # Rest-mass density on the right 1.0,0.125,1.0
        v =  0.0   # Velocity on the right
        p =  0.1  #(2/3)*(10^-6)   # Pressure on the right (1e-6)*(2/3) , 0.1 , 995
    end
    W = 1 / sqrt(1 - v^2)  # Lorentz factor
    #h = 1 + (γ * p) / ((γ - 1) * ρ)  # Specific enthalpy with equation of state  p = (γ - 1)(e - ρ)
    #h = (γ * p) / ((γ - 1) * ρ)  # Specific Enthalpy with equation of state p = (γ - 1)*p 
    
    
    e = (p / (γ - 1) ) + ρ
    
    D = ρ * W  # Relativistic mass density
    S = (e + p) * W^2 * v  # Momentum density
    E = (e + p) * W^2 - p   # Energy density
    
    if W > 1e8
        print("Warning the speed is very close to speed of light")
    end
    
    if E < D && E < S
        return println("Warning: Physical constraints violated in initial_condition at  x = $x")
    else
        return [D, S, E]
    end
end
x1 = 0.2
x2 = 0.8
println(initial_condition(x1))
println(initial_condition(x2))


function fixed_iter(g,x0,tol)
   i = 0
   x = g(x0)
    while true
        t = abs(x - x0)
        #println(t)
        i = i + 1
        if t < tol
            break 
        end
        x0 = x
        x = g(x0)
    end
    return([x, i, x - g(x)])
end



using Roots

δ = 1e-6




function pressure_nu(U)
    D, S, E = U
    
  δ = 1e-6
    ϵ = 1e-6

    # Define the pressure function as a function of v
    function p(v)
        e = E - v * S
        n = D * sqrt(1 - v^2)
        return (γ - 1) * (e - n)
    end

    # Define the function g(v) = S / (E + p(v))
    g(v) = S / (E + p(v))

    # Define the upper and lower bounds for v
    vu = min(1, S / E + δ)

    dsq = (γ * E)^2 - 4 * (γ - 1) * S^2

    if dsq < 0
        println("Discriminant is negative: unphysical solution")
    end

    d = sqrt(dsq)

    # Avoid division by zero
    if S == 0.0
        S += ϵ
    end

    vl = (1 / (2 * S * (γ - 1))) * (γ * E - d)

    # Initial guess
    v0 = 0.5 * (vl + vu)

    # Solve for v using fixed point iteration
    v_result = fixed_iter(g, v0, ϵ)
    vroot = v_result[1]  # Extract the value of v from the result

    # Now compute final quantities using the value of v
    e_val = E - vroot * S
    n_val = D * sqrt(1 - vroot^2)
    p_val = (γ - 1) * (e_val - n_val)

    return ([p_val, vroot , n_val , e_val])
end

println(" The left quantities ", pressure_nu(initial_condition(x1)))

println(" The right quantities ", pressure_nu(initial_condition(x2)))



function flux(U)
    D, S, E = U
    p = pressure_nu(U)[1]  # Pressure
    v = pressure_nu(U)[2]  # speed
    return [v * D, v * S + p, (p + E) * v]
end

flux(initial_condition(x1))
flux(initial_condition(x2))



δ1 = 1e-12

# RHLLE FLUX

function relativistic_hlle_flux(U_L,U_R)
    p_L = pressure_nu(U_L)[1]
    p_R = pressure_nu(U_R)[1]
    
    
    e_L = pressure_nu(U_L)[4]
    e_R = pressure_nu(U_R)[4]
    
    n_L = pressure_nu(U_L)[3]
    n_R = pressure_nu(U_R)[3]
    
    e_n_L = (e_L - n_L)
    e_n_R = (e_R - n_R)
    			
 
    c_s_L_2F = (γ*(γ - 1)*(e_n_L))/(n_L + γ*(e_n_L))  # Given in paper (also derived)
    c_s_R_2F = (γ*(γ - 1)*(e_n_R))/(n_R + γ*(e_n_R))  # Given in paper (also derived)	

 
 
    c_s_L_2S = (γ * p_L)/(n_L + (p_L + e_L))
    c_s_R_2S = (γ * p_R)/(n_R + (p_R + e_R))
    
    c_s_L_2 = max(c_s_L_2F, c_s_L_2S)
    c_s_R_2 = max(c_s_R_2F, c_s_R_2S)
    
    
    if c_s_L_2 < 0  ||  c_s_R_2 < 0
    	error("Unphysical Sound speeds detected", [c_s_L_2 , c_s_R_2 , p_L , p_R , e_n_L, e_n_R ])
    end
    
    
    c_s_L = sqrt(c_s_L_2)
    c_s_R = sqrt(c_s_R_2)
    
    c_s_M = 0.5*(c_s_L + c_s_R)
    
    v_L = pressure_nu(U_L)[2]
    v_R = pressure_nu(U_R)[2]
    
    v_M = 0.5*(v_L + v_R)
    
    λ_L = min(0, (v_M - c_s_M) / (1 - v_M * c_s_M), (v_L - c_s_L) / (1 - v_L * c_s_L))
    
    λ_R = max(0, (v_M + c_s_M) / (1 + v_M * c_s_M), (v_R + c_s_R) / (1 + v_R * c_s_R))
    
    F_L = flux(U_L)
    F_R = flux(U_R)
    
    if λ_L >= 0
        return F_L
    elseif λ_R <= 0
        return F_R
    else
        return (λ_R * F_L - λ_L * F_R + λ_L * λ_R * (U_R - U_L)) / (λ_R - λ_L)
    end
end







#MINMOD SLOPE LIMITER SCHEME

function minmod_slope_limiter(U)
    NC = length(U) - 2  # Number of cells (excluding ghost cells)
    UL = deepcopy(U)    # Initialize UL with the same structure as U
    UR = deepcopy(U)    # Initialize UR with the same structure as U

    for j in 2:NC+1  # Loop over interior cells (adjusting for 1-based indexing)
        Ujm1 = U[j-1]
        Uj   = U[j]
        Ujp1 = U[j+1]

        # Compute differences for each component (D, S, E)
        sp = Ujp1 .- Uj  # Element-wise difference
        sm = Uj .- Ujm1  # Element-wise difference

        # Apply minmod to each component
        dU = 0.25 .* (sign.(sp) .+ sign.(sm)) .* min.(abs.(sp), abs.(sm))

        # Compute limited states
        Ujp = Uj .+ dU
        Ujm = Uj .- dU

        UL[j] = Ujp
        UR[j-1] = Ujm
    end

    return UL, UR
end

using Base.Threads

# Create grid
x = LinRange(0, L, N)
U = [initial_condition(xi) for xi in x]
U_half_L = deepcopy(U)  # Left states at half time step
U_half_R = deepcopy(U)  # Right states at half time step
U_new = deepcopy(U)     # Final updated states

println(" The value of a left cell after slopelim", minmod_slope_limiter(U[50]))
println("The value of a right cell after slopelim",minmod_slope_limiter( U[200]))


t = 0.0

while t <= T
    
    global t
    # Step 1: Apply minmod slope limiter
    
    UL, UR = minmod_slope_limiter(U)
    

    constraint_violated = false  # Flag for constraint violation

    # Step 2: Update primitive variables at half time step
    @threads for i in 2:N-1
        U_L = UL[i]
        U_R = UR[i]
        
        FL = flux(U_L)
        FR = flux(U_R)
        
        # Update half-step states
        U_half_L[i] = U_L - 0.5 * (dt / dx) * (FR - FL)
        U_half_R[i] = U_R - 0.5 * (dt / dx) * (FR - FL)
        
  

        # Ensure physical validity
        if U_half_L[i][3] < U_half_L[i][1] || U_half_L[i][3] < U_half_L[i][2] ||
           U_half_R[i][3] < U_half_R[i][1] || U_half_R[i][3] < U_half_R[i][2]
            println("Physical constraint violated at half time step of $t and index $i")
            constraint_violated = true
        end
    end

    if constraint_violated
        break
    end

    # Step 3: Apply RHLLE flux for final update
    constraint_violated = false
    @threads for j in 2:N-1
        U_L = U_half_L[j]
        U_R = U_half_R[j]
        
        F_L = relativistic_hlle_flux(U_half_L[j-1], U_half_R[j-1])
        F_R = relativistic_hlle_flux(U_L, U_R)
        
        U_new[j] = U[j] - (dt / dx) * (F_R - F_L)

        # Ensure physical validity
        if U_new[j][3] < U_new[j][1] || U_new[j][3] < U_new[j][2]
            println("Physical constraint violated at time $t and index $j")
            constraint_violated = true
        end
    end

    if constraint_violated
        break
    end

    # Apply transmissive boundary conditions
    U_new[1] = U_new[2]
    U_new[end] = U_new[end-1]

    # Update solution for next time step
    U .= U_new
    
    v_final = [pressure_nu(U[i])[2] for i in 1:N]
    
#    global dt = CFL * dx / maximum(v_final)

    global t = t + dt 
    println("Time step :", t , " dt for next time step :", dt)
end


D_final = [U[i][1] for i in 1:N]

p_final = [pressure_nu(U[i])[1] for i in 1:N]

v_final = [pressure_nu(U[i])[2] for i in 1:N]

#e_final = [pressure_nu(U[i])[3] for i in 1:N];




println("THE CFL LIMIT: ", maximum(v_final) * (dt / dx))


plot1 = scatter(x[1:2:end], D_final[1:2:end], label="Relativistic Density (D)", xlabel="x", ylabel="Value", title="R R P Solution at t = $T, γ = $γ", markersize=2)

plot2 = scatter(x[1:2:end], v_final[1:2:end], label="Velocity (v)", markersize=2)

plot3 = scatter(x[1:2:end], p_final[1:2:end], label="Pressure (p)", markersize=2)

#plot4 = scatter(v_final[1:2:end], p_final[1:2:end], label="Pressure (p)", markersize=2)

plot(plot1, plot2, plot3, layout=(1,3), size=(1200,400))

savefig("RHLLE_1D_SCHEME.svg") 




