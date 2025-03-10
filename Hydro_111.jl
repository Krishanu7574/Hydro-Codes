using Plots

# Parameters

#γ = 4/3  # Adiabatic index for relativistic gas 5/3 , 4/3 , 1.999


println("Enter the value of γ : ")
γ = parse(Float64, readline())

c = 1.0  # Speed of light (set to 1 for simplicity)
N = 250  # Number of grid points
L = 1.0  # Domain length
L_mid = L / 2

N_mid = N / 2
dx = L / N  # Grid spacing
dt = 0.0001  # Time step
T = 0.25  # Total simulation time
CFL = 0.01;  # CFL condition
println(L_mid)
println(N_mid)



# Initial condition 
function initial_condition(x)
    if x < L_mid
        ρ =  1.0  # Rest-mass density on the left  1.0, 10.0,2.0
        v =  0.0  # Velocity on the left
        p =  1.0  # Pressure on the left   1.0, 40.0/3 ,1990
    else
        ρ =  0.125  # Rest-mass density on the right 1.0,0.125,1.0
        v =  0.0    # Velocity on the right
        p =  0.1   # Pressure on the right 0.00001*(2/3) , 0.1 , 995
    end
    W = 1 / sqrt(1 - v^2)  # Lorentz factor
    h = 1 + (γ * p) / ((γ - 1) * ρ)  # Specific enthalpy
    D = ρ * W  # Relativistic mass density
    S = ρ * h * W^2 * v  # Momentum density
    E = ρ * h * W^2 - p   # Energy density
    
    if W > 1e8
        print("Warning the speed is very close to speed of light")
    end
    
    if E < D && E < S
        return println("Warning: Physical constraints violated in initial_condition at  x = $x")
    else
        return [D, S, E]
    end
end
x1 = 0.8
x2 = 0.2
println(initial_condition(x1))
println(initial_condition(x2))




#ROOT FINDING

using Roots

δ = 1e-6

function pressure_nu(U)
    D,S,E = U 
    function g(v)
        return ((v * E - S) + (γ - 1) * S * (1 - v^2)*(1 / γ))^2 - v^2 * (1 - v^2) * ((γ - 1)*D - (1/γ)*D*(γ - 1)^2)^2
        #return  ((v * E - S) + (γ - 1) * S * (1 - v^2)*(1 / γ))^2
    end
    
    vu = min(1,S/E + δ)
    
#    vl = 0.0
    
    dsq = (γ*E)^2 - 4*(γ - 1)*S^2
    
    if dsq < 0
        print("discreminant is negative unphysical solution")
    end
    d = sqrt(dsq)
    
    ϵ = 1e-6
    
    if S == 0.0
        S = (S + ϵ)
    else
        S = S
    end
    
    vl = (1/(2*S*(γ - 1))) * (γ*E - d)

    
    #v0 = 0.5*(1 - D/E)*(vl - vu) + 0.5*(vl + vu) 
    
    v0 =  0.5*(vl + vu) 
    
    
    
#    vroot = find_zero(g, v0, maxiters = 20)

    vroot = find_zero(g, v0, atol=1e-6)
    e = E - vroot*S
    p = (γ - 1)*(E - vroot*S) 
    return [p,vroot,e,v0,vl,d]
end

println(pressure_nu(initial_condition(x1)))

println(pressure_nu(initial_condition(x2)))




# FLUX 
function flux(U)
    D, S, E = U
    p = pressure_nu(U)[1]  # Pressure
    v = pressure_nu(U)[2]
    return [v * D, v * S + p, (p + E) * v]
end

flux(initial_condition(x1))


# RHLLE FLUX

function relativistic_hlle_flux(U_L,U_R)
    p_L = pressure_nu(U_L)[1]
    p_R = pressure_nu(U_R)[1]
    
    c_s_L_2 = γ * p_L / (U_L[1] + (γ/(γ-1))*p_L)
    c_s_R_2 = γ * p_R / (U_R[1] + (γ/(γ-1))*p_R)
    
    if c_s_L_2 < 0  ||  c_s_R_2 < 0
        print("Unphysical Sound speeds")
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




# TIME EVOLUTION

using Base.Threads

# Create grid
x = LinRange(0, L, N)
U = [initial_condition(xi) for xi in x]
U_half_L = deepcopy(U)  # Left states at half time step
U_half_R = deepcopy(U)  # Right states at half time step
U_new = deepcopy(U)     # Final updated states

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
            println("Physical constraint violated at time $t and index $i")
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
    
    global dt = CFL * dx / maximum(v_final)

    global t = t + dt 
    println(t,dt)
end


D_final = [U[i][1] for i in 1:N]

p_final = [pressure_nu(U[i])[1] for i in 1:N]

v_final = [pressure_nu(U[i])[2] for i in 1:N]

e_final = [pressure_nu(U[i])[3] for i in 1:N];




println("THE CFL LIMIT: ", maximum(v_final) * (dt / dx))


scatter(x[1:2:end], D_final[1:2:end], label="(1/10)Relativistic Density (D)", xlabel="x", ylabel="Value", title="Relativistic Riemann Problem Solution at t = $T, γ = $γ", markersize=2)

scatter!(x[1:2:end], v_final[1:2:end], label="Velocity (v)", markersize=2)

scatter!(x[1:2:end], p_final[1:2:end], label="(1/20)Pressure (p)", markersize=2)

savefig("RHLLE_1D_SCHEME.svg") 

