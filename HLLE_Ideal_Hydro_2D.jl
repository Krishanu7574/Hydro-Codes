using CairoMakie
using Roots
using LinearAlgebra

# ============================================================
#  (2+1)D Ideal Relativistic Hydrodynamics — HLLE + SSPRK3
#
#  State vectors (conservative):  U  = [N, Mx, My, E]
#  Primitive variables:           u  = [n, vx, vy, pg]
#
#  N  = n γ                         (baryon number density)
#  Mx = γ²(e+p) vx                  (x-momentum)
#  My = γ²(e+p) vy                  (y-momentum)
#  E  = γ²(e+p) − p                 (energy density)
#  e  = p/(Γ−1) + n                 (total energy density)
#  γ  = 1/√(1 − v²),  v² = vx²+vy²
#
#  Dimensional splitting (Strang):
#      L(dt) = Lx(dt/2) ∘ Ly(dt) ∘ Lx(dt/2)
#  Each 1-D sweep reuses the 1D HLLE flux with the
#  appropriate velocity component treated as the "normal".
# ============================================================

# ── Grid parameters ─────────────────────────────────────────
const Nx   = 200          # cells in x
const Ny   = 200          # cells in y
const Γ    = 1.666          # adiabatic index  
const xL   = -0.5;  const xR = 0.5
const yL   = -0.5;  const yR = 0.5
const CFL  = 0.35
const tf   = 0.35         # final time

const dx = (xR - xL) / Nx
const dy = (yR - yL) / Ny

const small = 1e-12
const big   = 1e12

x = LinRange(xL + dx/2, xR - dx/2, Nx)   # cell centres
y = LinRange(yL + dy/2, yR - dy/2, Ny)

# ── Lorentz factor  (4-velocity safe) ──────────────
function Lorentz(vx, vy)
    v2 = vx*vx + vy*vy
    if v2 >= 1 - small
        @warn "v² ≥ 1; clamping"
        v2 = clamp(v2, small, 1 - small)
    end
    return clamp(sqrt(1 / (1 - v2)), 1.0, big)
end




# ── Primitive → Conservative ─────────────────────────────────
# Returns [N, Mx, My, E]
function prim_to_cons(u)
    n, vx, vy, pg = u
    γ  = Lorentz(vx, vy)
    e  = pg / (Γ - 1) + n
    h  = e + pg            # enthalpy density  (= e + p)
    N  = n * γ
    Mx = γ^2 * h * vx
    My = γ^2 * h * vy
    E  = γ^2 * h - pg
    return [N, Mx, My, E]
end

# ── Primitive → x-direction flux ─────────────────────────────
# Fx = [N vx,  Mx vx + p,  My vx,  (E+p) vx]
function prim_to_flux_x(u)
    n, vx, vy, pg = u
    γ  = Lorentz(vx, vy)
    e  = pg / (Γ - 1) + n
    h  = e + pg
    N  = n * γ
    Mx = γ^2 * h * vx
    My = γ^2 * h * vy
    E  = γ^2 * h - pg
    return [N*vx,  Mx*vx + pg,  My*vx,  (E + pg)*vx]
end

# ── Primitive → y-direction flux ─────────────────────────────
# Fy = [N vy,  Mx vy,  My vy + p,  (E+p) vy]
function prim_to_flux_y(u)
    n, vx, vy, pg = u
    γ  = Lorentz(vx, vy)
    e  = pg / (Γ - 1) + n
    h  = e + pg
    N  = n * γ
    Mx = γ^2 * h * vx
    My = γ^2 * h * vy
    E  = γ^2 * h - pg
    return [N*vy,  Mx*vy,  My*vy + pg,  (E + pg)*vy]
end


# ── Conservative → Primitive (2D root-find) ──────────────────
# Unknowns: |v| (scalar speed).  Direction recovered from Mx,My.

function cons_to_prim(U)
    N, Mx, My, E = U
    M = sqrt(Mx^2 + My^2)

    # Near-rest state
    if M < small
        n = N
        p = (Γ - 1) * (E - n)
        return [n, 0.0, 0.0, p]
    end

    function f(v)
        # 0 < v < 1
        n = N * sqrt(1 - v^2)
        e = E - v * M
        p = (Γ - 1) * (e - n)
        return v * (E + p) - M
    end

    #v = find_zero(f, (small, 1 - small), atol=small)
    v = find_zero(f, 0.5 ; atol=1e-8, rtol=1e-8)

    if abs(f(v)) > 1e-8 # small
        println("root not within tolerance: ", [v, f(v)])
    end
    #println(f(v))
    n = N * sqrt(1 - v^2)
    e = E - v * M
    p = (Γ - 1) * (e - n)

    # Recover direction from momentum vector
    vx = v * Mx / M
    vy = v * My / M

    return [n, vx, vy, p]
end


u1 = [0.125, 0.2 , 0.3, 0.1]
u2 = [1.0, 0.0, 0.0, 30.0]
cons_to_prim(prim_to_cons(u2))



# ── Sound speed ──────────────────────────────────────────────
function sound_speed(n, pg)
    e   = pg / (Γ - 1) + n
    cs2 = Γ * (Γ - 1) * (e - n) / (n + Γ * (e - n))
    cs2 = clamp(cs2, small, 1 - small)
    return cs2, sqrt(cs2)
end

# ── HLLE wave-speed bounds (direction d ∈ {:x,:y}) ──────────
# vn = normal velocity component
function HLLE_speeds(uL, uR, dir::Symbol)
    n_L, vx_L, vy_L, pg_L = uL
    n_R, vx_R, vy_R, pg_R = uR

    vn_L = (dir == :x) ? vx_L : vy_L
    vn_R = (dir == :x) ? vx_R : vy_R

    _, cs_L = sound_speed(n_L, pg_L)
    _, cs_R = sound_speed(n_R, pg_R)

    vn_bar = 0.5 * (vn_L + vn_R)
    cs_bar = 0.5 * (cs_L + cs_R)

    λL = min(0.0,
             (vn_bar - cs_bar) / (1 - vn_bar * cs_bar),
             (vn_L   - cs_L)   / (1 - vn_L   * cs_L))
    λR = max(0.0,
             (vn_bar + cs_bar) / (1 + vn_bar * cs_bar),
             (vn_R   + cs_R)   / (1 + vn_R   * cs_R))
    return λL, λR
end

# ── HLLE numerical flux (direction d) ───────────────────────
function HLLE_flux(uL, uR, dir::Symbol)
    SL, SR = HLLE_speeds(uL, uR, dir)
    FL = (dir == :x) ? prim_to_flux_x(uL) : prim_to_flux_y(uL)
    FR = (dir == :x) ? prim_to_flux_x(uR) : prim_to_flux_y(uR)
    UL = prim_to_cons(uL)
    UR = prim_to_cons(uR)

    if SL >= 0
        return FL
    elseif SR <= 0
        return FR
    else
        return (SR .* FL .- SL .* FR .+ SL .* SR .* (UR .- UL)) ./ (SR - SL)
    end
end

# ── Global maximum signal speed ──────────────────────────────
function maximum_signal(u_prim)
    S_max = 0.0
    for j in 1:Ny, i in 1:Nx
        n, vx, vy, pg = u_prim[i, j]
        _, cs = sound_speed(n, pg)
        v    = sqrt(vx^2 + vy^2)
        # Relativistic bound: max wave speed in any direction
        SR = (v + cs) / (1 + v * cs)
        S_max = max(S_max, SR)
    end
    return S_max
end
# ── Minmod slope limiter (2D, per axis) ─────────────────────
# Returns left/right reconstructed states on each cell edge.
# Operates on the 2D array of primitives.
function minmod_limiter_x(u_prim)
    uL = deepcopy(u_prim)   # right edge of cell (i)
    uR = deepcopy(u_prim)   # left  edge of cell (i+1)

    NVARS = 4
    for j in 1:Ny, i in 2:Nx-1
        for k in 1:NVARS
            sp = u_prim[i+1, j][k] - u_prim[i, j][k]
            sm = u_prim[i, j][k]   - u_prim[i-1, j][k]
            dU = (sp * sm > 0.0) ? (abs(sp) < abs(sm) ? 0.5*sp : 0.5*sm) : 0.0
            uL[i, j][k] = u_prim[i, j][k] + dU   # right interface of cell i
            uR[i, j][k] = u_prim[i, j][k] - dU   # left  interface of cell i (seen from i+1)
        end
    end
    # uL[i] is reconstructed state at x_{i+1/2}⁻
    # uR[i] is reconstructed state at x_{i-1/2}⁺  (but we only need the half-step below)
    # We return: for interface (i+½): left state = uL[i], right state = uR[i+1]
    return uL, uR
end

function minmod_limiter_y(u_prim)
    uL = deepcopy(u_prim)
    uR = deepcopy(u_prim)

    NVARS = 4
    for j in 2:Ny-1, i in 1:Nx
        for k in 1:NVARS
            sp = u_prim[i, j+1][k] - u_prim[i, j][k]
            sm = u_prim[i, j][k]   - u_prim[i, j-1][k]
            dU = (sp * sm > 0.0) ? (abs(sp) < abs(sm) ? 0.5*sp : 0.5*sm) : 0.0
            uL[i, j][k] = u_prim[i, j][k] + dU
            uR[i, j][k] = u_prim[i, j][k] - dU
        end
    end
    return uL, uR
end





# ── Outflow boundary conditions ──────────────────────────────
function apply_bc!(U)
    for j in 1:Ny
        U[1,   j] = U[2,   j]
        U[Nx,  j] = U[Nx-1,j]
    end
    for i in 1:Nx
        U[i,  1 ] = U[i,  2 ]
        U[i,  Ny] = U[i, Ny-1]
    end
end


# ── SSPRK3 step (one full dt, both directions) ───────────────
# Uses Strang splitting:  Lx(dt/2) → Ly(dt) → Lx(dt/2)
# Each Lx/Ly sub-step is a full SSPRK3 in that direction.
# ============================================================
#  Unsplit SSPRK3 for (2+1)D relativistic hydrodynamics
#
#  Replaces: ssprk3_sweep_x, ssprk3_sweep_y, full_step
#
#  The combined RHS  L(U) = −∂Fx/∂x − ∂Fy/∂y
#  is evaluated at each stage with the SAME intermediate state,
#  so x and y fluxes are always consistent with each other.
#  No splitting error.  No half-step bookkeeping.
#
#  SSPRK3 Butcher tableau (Shu-Osher form):
#    U¹       = Uⁿ + dt · L(Uⁿ)
#    U²       = (3/4) Uⁿ + (1/4)[ U¹ + dt · L(U¹) ]
#    Uⁿ⁺¹    = (1/3) Uⁿ + (2/3)[ U² + dt · L(U²) ]
# ============================================================

# ── Combined RHS: both flux directions together ──────────────
# Returns dU[i,j] = −(Fx_{i+½} − Fx_{i−½})/dx
#                   −(Fy_{j+½} − Fy_{j−½})/dy
function rhs_2d(u_prim)
    dU = [zeros(4) for _ in 1:Nx, _ in 1:Ny]

    # ── x-direction contribution ──────────────────────────────
    uL_x, uR_x = minmod_limiter_x(u_prim)

    for j in 1:Ny
        # interface fluxes  F[i] lives at  x_{i+½}
        F = Vector{Vector{Float64}}(undef, Nx - 1)
        for i in 1:Nx-1
            F[i] = HLLE_flux(uL_x[i, j], uR_x[i+1, j], :x)
        end
        for i in 2:Nx-1
            dU[i, j] .+= -(F[i] .- F[i-1]) ./ dx
        end
    end

    # ── y-direction contribution ──────────────────────────────
    uL_y, uR_y = minmod_limiter_y(u_prim)

    for i in 1:Nx
        # interface fluxes  G[j] lives at  y_{j+½}
        G = Vector{Vector{Float64}}(undef, Ny - 1)
        for j in 1:Ny-1
            G[j] = HLLE_flux(uL_y[i, j], uR_y[i, j+1], :y)
        end
        for j in 2:Ny-1
            dU[i, j] .+= -(G[j] .- G[j-1]) ./ dy
        end
    end

    return dU
end





# ── Unsplit SSPRK3 full timestep ─────────────────────────────
function full_step(U, u_prim, dt)

    # ── Stage 1 ───────────────────────────────────────────────
    L0  = rhs_2d(u_prim)
    U1  = [U[i,j] .+ dt .* L0[i,j]  for i in 1:Nx, j in 1:Ny]

    apply_bc!(U1)
    up1 = [cons_to_prim(U1[i,j])     for i in 1:Nx, j in 1:Ny]

    # ── Stage 2 ───────────────────────────────────────────────
    L1  = rhs_2d(up1)
    U2  = [0.75 .* U[i,j] .+
           0.25 .* (U1[i,j] .+ dt .* L1[i,j])
           for i in 1:Nx, j in 1:Ny]

    apply_bc!(U2)
    up2 = [cons_to_prim(U2[i,j])     for i in 1:Nx, j in 1:Ny]

    # ── Stage 3 ───────────────────────────────────────────────
    L2  = rhs_2d(up2)
    U3  = [(1/3) .* U[i,j] .+
           (2/3) .* (U2[i,j] .+ dt .* L2[i,j])
           for i in 1:Nx, j in 1:Ny]

    apply_bc!(U3)
    return U3
end



# NOTE: ssprk3_sweep_x and ssprk3_sweep_y are no longer needed
# and can be removed from your file.


# ── Initial condition: 2D cylindrical blast wave ─────────────
# (generalisation of the Sod tube to 2D)
# Inner high-pressure region: r < 0.25
# Outer low-pressure region:  r ≥ 0.25
function initial_condition()
    u_prim = Matrix{Vector{Float64}}(undef, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        r = sqrt(x[i]^2 + y[j]^2)
        if r < 0.25
            # high density/pressure inner state
            u_prim[i, j] = [1.0, 0.0, 0.0, 30.0]
        else
            # low density/pressure outer state
            u_prim[i, j] = [1.0,  0.0, 0.0, 1.0]
        end
    end
    return u_prim
end

# ── Simulation ───────────────────────────────────────────────
function run_simulation()
    u_prim = initial_condition()
    U      = [prim_to_cons(u_prim[i,j]) for i in 1:Nx, j in 1:Ny]

    t         = 0.0
    iteration = 0

    println("Starting (2+1)D RHLLE simulation  Nx=$(Nx)  Ny=$(Ny)  tf=$(tf)")

    while t < tf
        S_max = maximum_signal(u_prim)
        dt    = CFL * min(dx, dy) / S_max
        dt    = min(dt, tf - t)

        U = full_step(U, u_prim, dt)
        apply_bc!(U)
        u_prim = [cons_to_prim(U[i,j]) for i in 1:Nx, j in 1:Ny]

        t         += dt
        iteration += 1

        if iteration % 20 == 0
            println("  iter=$(iteration)  t=$(round(t,digits=5))  dt=$(round(dt,digits=7))")
        end
    end

    println("Done!  iterations=$(iteration)  t=$(t)")
    return u_prim, U
end

# ── Run ──────────────────────────────────────────────────────
u_prim_final, U_final = run_simulation()

# ── Plot final state ─────────────────────────────────────────
n_field   = [u_prim_final[i,j][1] for i in 1:Nx, j in 1:Ny]
vx_field  = [u_prim_final[i,j][2] for i in 1:Nx, j in 1:Ny]
vy_field  = [u_prim_final[i,j][3] for i in 1:Nx, j in 1:Ny]
pg_field  = [u_prim_final[i,j][4] for i in 1:Nx, j in 1:Ny]
γ_field   = [Lorentz(u_prim_final[i,j][2], u_prim_final[i,j][3]) for i in 1:Nx, j in 1:Ny]
v_field   = sqrt.(vx_field.^2 .+ vy_field.^2)

fig = Figure(resolution = (1200, 900))

ax1 = Axis(fig[1,1], title="Density n",         xlabel="x", ylabel="y", aspect=1)
ax2 = Axis(fig[1,2], title="Pressure p",         xlabel="x", ylabel="y", aspect=1)
ax3 = Axis(fig[2,1], title="Speed |v|",          xlabel="x", ylabel="y", aspect=1)
ax4 = Axis(fig[2,2], title="Lorentz factor γ",   xlabel="x", ylabel="y", aspect=1)

xv = collect(x)
yv = collect(y)

heatmap!(ax1, xv, yv, n_field,  colormap=:inferno)
heatmap!(ax2, xv, yv, pg_field, colormap=:plasma)
heatmap!(ax3, xv, yv, v_field,  colormap=:viridis)
heatmap!(ax4, xv, yv, γ_field,  colormap=:magma)

Colorbar(fig[1,1][1,2], ax1.scene.plots[1])
Colorbar(fig[1,2][1,2], ax2.scene.plots[1])
Colorbar(fig[2,1][1,2], ax3.scene.plots[1])
Colorbar(fig[2,2][1,2], ax4.scene.plots[1])

display(fig)
#save("Hydro_RHLLE_2D_SSPRK3.png", fig)
println("Saved: Hydro_RHLLE_2D_SSPRK3_1.png")

# ── Conservation diagnostics ─────────────────────────────────
function conservation_errors(u_prim_0, U_final)
    U0 = [prim_to_cons(u_prim_0[i,j]) for i in 1:Nx, j in 1:Ny]
    dA = dx * dy

    mass0  = sum(U[1] for U in U0)  * dA
    mom0x  = sum(U[2] for U in U0)  * dA
    mom0y  = sum(U[3] for U in U0)  * dA
    ener0  = sum(U[4] for U in U0)  * dA

    massF  = sum(U[1] for U in U_final) * dA
    momFx  = sum(U[2] for U in U_final) * dA
    momFy  = sum(U[3] for U in U_final) * dA
    enerF  = sum(U[4] for U in U_final) * dA

    println("\nConservation errors:")
    println("  Mass:       $(abs(massF-mass0)/abs(mass0)*100) %")
    println("  Momentum x: $(abs(momFx-mom0x)/(abs(mom0x)+small)*100) %")
    println("  Momentum y: $(abs(momFy-mom0y)/(abs(mom0y)+small)*100) %")
    println("  Energy:     $(abs(enerF-ener0)/abs(ener0)*100) %")
end

conservation_errors(initial_condition(), U_final)
