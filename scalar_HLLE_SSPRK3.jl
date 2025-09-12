using CairoMakie

# ---------------- grid & arrays ----------------
const ng = 2
N   = 400
Nx  = N + 2ng
xL, xR = -1.0, 1.0
dx = (xR - xL)/N
x  = range(xL + dx/2, xR - dx/2, length=N)

u_in = zeros(Float64,Nx)       # initial condition (with ghosts)
u    = zeros(Float64,Nx)       # evolved solution (with ghosts)
Ff   = zeros(Float64,Nx-1)     # face-centered fluxes



# ---------------- slope limiter ----------------
minmod(a,b) = (sign(a)==sign(b)) ? sign(a)*min(abs(a),abs(b)) : 0.0
function slopes_minmod(u)
    s = similar(u)
    @inbounds for j in 2:length(u)-1
        s[j] = minmod(u[j]-u[j-1], u[j+1]-u[j])
    end
    s[1]=0.0; s[end]=0.0
    return s
end

# ---------------- superbee slope limiter ----------------
function superbee(a, b)
    if sign(a) == sign(b)
        sa = sign(a)
        return sa * max(0.0,
                        min(2abs(a), abs(b)),
                        min(abs(a), 2abs(b)))
    else
        return 0.0
    end
end


# ---------------- MC slope limiter ----------------
function mc(a, b)
    if sign(a) == sign(b)
        sa = sign(a)
        return sa * min(2abs(a), 2abs(b), 0.5*(abs(a) + abs(b)))
    else
        return 0.0
    end
end



function slopes_mc(u)
    s = similar(u)
    @inbounds for j in 2:length(u)-1
        s[j] = mc(u[j]-u[j-1], u[j+1]-u[j])
    end
    s[1]=0.0; s[end]=0.0
    return s
end



function slopes_superbee(u)
    s = similar(u)
    @inbounds for j in 2:length(u)-1
        s[j] = superbee(u[j]-u[j-1], u[j+1]-u[j])
    end
    s[1]=0.0; s[end]=0.0
    return s
end


# ---------------- physics: u_t + a u_x = 0 ----------------
const a = 1.0
flux(u) = a*u

# HLLE flux solver
hlle(UL, UR, FL, FR, λL, λR) = (λL ≥ 0)  ? FL :
                               (λR ≤ 0)  ? FR :
                               (λR*FL - λL*FR + λL*λR*(UR - UL)) / (λR - λL)

# Characteristic speeds for advection
wavespeeds_advection(uL,uR) = (min(0.0,a), max(0.0,a))

# ---------------- boundary fill (outflow) ----------------
function fill_bc!(u)
    u[1:ng]        .= u[ng+1]
    u[end-ng+1:end].= u[end-ng]
end

# ---------------- RHS operator (semi-discrete) ----------------
function rhs_hlle(u, dx)
    fill_bc!(u)
    #s = slopes_minmod(u)
    #s = slopes_superbee(u)
    s = slopes_mc(u)
    dudt = zeros(size(u))

    @inbounds for i in 1:length(u)-1
        uL = u[i]   + 0.5*s[i]
        uR = u[i+1] - 0.5*s[i+1]
        λL, λR = wavespeeds_advection(uL,uR)
        FL, FR = flux(uL), flux(uR)
        Ff[i] = hlle(uL, uR, FL, FR, λL, λR)
    end

    @inbounds for j in ng+1:length(u)-ng
        dudt[j] = -(Ff[j] - Ff[j-1]) / dx
    end

    return dudt
end

# ---------------- SSPRK3 integrator ----------------
function ssprk3_step(u, dt, dx)
    k1 = rhs_hlle(u, dx)
    u1 = u .+ dt .* k1

    k2 = rhs_hlle(u1, dx)
    u2 = (3/4) .* u .+ (1/4) .* (u1 .+ dt .* k2)

    k3 = rhs_hlle(u2, dx)
    return (1/3) .* u .+ (2/3) .* (u2 .+ dt .* k3)
end

# ---------------- initial condition ----------------
u[ng+1:ng+N] .= (x .≥ -0.4) .& (x .≤ -0.1)  # block pulse
u_in .= u
fill_bc!(u)

# ---------------- time loop ----------------
t, tf = 0.0, 0.6
CFL = 0.5
cmax = abs(a)

while t < tf
    dt = CFL * dx / cmax
    if t + dt > tf
        dt = tf - t
    end
    u .= ssprk3_step(u, dt, dx)
    global t += dt
end

# ---------------- plot ----------------
x_list = range(xL, xR, length(u))
# plot(x_list, u, label="final t=$(round(tf,digits=2))")
# plot!(x_list, u_in, label="initial t=0", linestyle=:dash)
# xlabel!("x"); ylabel!("u")
# title!("1D Advection: HLLE + slope limiter + SSPRK3")
# plot!(x_list, slopes(u),label = "Slope_limiter_ufinal")

fig = Figure(resolution = (800, 400))
ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "u", title = "1D Advection: HLLE + slope limiter + SSPRK3")
lines!(ax1, x_list, u, label="final t=$(round(tf,digits=2))", color = :blue)
lines!(ax1, x_list, u_in, label="initial t=0", linestyle=:dash, color = :red)
axislegend(ax1; position = :lt)
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel = "Slope", title = "Slope of u_final")
lines!(ax2, x_list, slopes(u), label = "Slope of u_final", color = :green)
axislegend(ax2; position = :lt)
fig