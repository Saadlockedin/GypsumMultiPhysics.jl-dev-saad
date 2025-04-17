using GypsumMultiPhysics
using Printf
using LaTeXStrings
using GLMakie, GridVisualize
using Plots 
using XLSX
using ExtendableGrids

FE=create_FE()

## generate the grid
nx      = 60;
Lx      = 0.012; # [m] thickness of gympson board
y       = 0.6;  # [m]
z       = 0.6;  # [m]
Volume  = Lx*y*z; # [m3]
ρ = 1;
T = 2; 
SpecInd = [ρ,T]
m = simplexgrid(0:Lx/nx:Lx); # create grid
gridplot(m;Plotter=GLMakie,resolution = (600, 250), legend =:rt)

function storage!(f, u, node, data)
    f[ρ] =  Φ* u[ρ]                     # Storage for mass concentration of kth gas phase component
    f[T] = ρs* cs +Φ *(ρ*cp) * u[T]   # Storage for temperature

    return nothing
end

function flux!(f, u, edge, data)
    # Convective and diffusive flux for ρ


    # Convective and diffusive flux for T


    return nothing
end

function reaction!(f, u, node, data)
    # Source term for mass concentration
    f[ρ] = Qₘ

    # Source term for temperature
    f[T] = Φ*u[Pₜ]+Qₜ

    return nothing
end


# time parameters
incident_HFlux = 40000; # kw
fname2         = "TypeSB_40K.xlsx";
time_max       = 1400;    # min
dt             = 0.05; 
savetime       = 10;       # save data at each 30 seconds
urelax         = 1;
## generate the grid
nx      = 60;
Lx      = 0.012; # [m] thickness of gympson board
y       = 0.6;  # [m]
z       = 0.6;  # [m]
Volume  = Lx*y*z; # [m3]    

m = simplexgrid(0:Lx/nx:Lx); # create grid
gridplot(m;Plotter=GLMakie,resolution = (600, 250), legend =:rt)

function storage!(f, u, node, data)
    f[c₁] =  ε * ρₐ * u[c₁]                     # Storage for oxygen concentration
    f[T] = (1 -  ε) * ρ_c *  C_pc * u[T]   # Storage for temperature
    f[c₂] = u[c₂]                                   # Storage for adsorbed oxygen
    return nothing
end

function flux!(f, u, edge, data)
    # Convective and diffusive flux for c1
    vel_c1 = project(edge, ρₐ *  V)
    Bplus =  ε *  Dₐ * ρₐ * fbernoulli(vel_c1 / ( ε *  Dₐ * ρₐ))
    Bminus =  ε *  Dₐ * ρₐ * fbernoulli(-vel_c1 / ( ε *  Dₐ * ρₐ))
    f[c₁] = Bminus * u[c₁, 1] - Bplus * u[c₁, 2]

    # Convective and diffusive flux for T
    vel_T = project(edge, ρₐ *  C_pa *  V)
    Bplus =  λₛ * fbernoulli(vel_T /  λₛ)
    Bminus =  λₛ * fbernoulli(-vel_T /  λₛ)
    f[T] = Bminus * u[T, 1] - Bplus * u[T, 2]

    return nothing
end

function reaction!(f, u, node, data)
    # Reaction rate for adsorbed oxygen
    r_ads =  -Bₐ*u[c₁]^(n) * f_c2 *abs(exp(- Eₐ /  R * (1 / u[T] - 1 /  T₀)))

    # Source term for oxygen concentration
    f[c₁] = - (1 -  ε) * ρ_c * r_ads

    # Source term for temperature
    f[T] = (1 -  ε) * ρ_c *  ΔH /  MW_O2 * r_ads

    # Source term for adsorbed oxygen
    f[c₂] = r_ads
    return nothing
end

function bcondition!(f, u, bnode, data)
    boundary_dirichlet!(f, u, bnode, species = c₁, region = 1, value =  c1₀)
    boundary_dirichlet!(f, u, bnode, species = c₁, region = 2, value = 0.0)
    boundary_dirichlet!(f, u, bnode, species = T, region = 1, value =  T₀)
    boundary_dirichlet!(f, u, bnode, species = T, region = 2, value =  T₀)
    return nothing
end

data = (
    Bₐ = Bₐ,
    C_pc = C_pc,
    C_pa = C_pa,
    Dₐ = Dₐ,
    ε = ε,
    Eₐ = Eₐ,
    R = R,
    V = V,
    λₛ = λₛ,
    ΔH = ΔH,
    MW_O2 = MW_O2,
    T₀ = T₀,
    c1₀ = c1₀,
    n = n,
    f_c2 = f_c2,
)

sys = VoronoiFVM.System(grid;data = data,
flux       = flux!,
storage    = storage!,
reaction   = reaction!,
bcondition = bcondition!,
unknown_storage = :sparse,
assembly = :edgewise, species = SpecInd)


inival = unknowns(sys)
inival[c₁, :] .=  c1₀  # Initial oxygen concentration
inival[T, :] .=  T₀    # Initial temperature
inival[c₂, :] .= 0.0       # Initial adsorbed oxygen


# Simulation Parameters
tend = Days * 24 * 60 * 60  # Total simulation time (s)
problem = ODEProblem(sys, inival, (0.0, tend))

# Solver
sol = solve(problem, Rosenbrock23(); dt = Δt, reltol = 1e-5, abstol = 1e-5)
#sol=solve(problem,ImplicitEuler(),  

                                   #force_dtmin=true,
                                   #adaptive=true,
                                   #reltol=1.0e-3,
                                   #abstol=1.0e-3,
                                   #initializealg=NoInit(),
                                  #dtmin=1e-6,

                                   #)

sol=reshape(sol,sys)

# Visualization setup
vis = GridVisualizer(Plotter = GLMakie, resolution = (1200, 400), layout = (1, 3), legend = :rb)

# Final time visualization
scalarplot!(vis[1, 1], grid, sol[end][c₁, :], label = "Oxygen Concentration (c₁)")
scalarplot!(vis[1, 2], grid, sol[end][T, :], label = "Temperature (T)")
scalarplot!(vis[1, 3], grid, sol[end][c₂, :], label = "Adsorbed Oxygen (c₂)")
reveal(vis)

# Temporal Evolution
time_points = [139, 278, 471] .* 24 .* 60 .* 60  # Days in seconds

for t in time_points
    # Plot oxygen concentration
    scalarplot!(vis[1, 1], grid, sol(t)[c₁, :], 
                label = "c₁ at day $(t / (24 * 60 * 60))", clear = false)

    # Plot temperature
    scalarplot!(vis[1, 2], grid, sol(t)[T, :], 
                label = "T at day $(t / (24 * 60 * 60))", clear = false)

    # Plot adsorbed oxygen
    scalarplot!(vis[1, 3], grid, sol(t)[c₂, :], 
                label = "c₂ at day $(t / (24 * 60 * 60))", clear = false)
end

# Load MATLAB data
mat_data = matread("Coal_Stockpile_Case_A_Final.mat")

# Extract MATLAB experimental data
A_temp_139 = mat_data["A_temp_139"]
A_temp_278 = mat_data["A_temp_278"]
A_temp_471 = mat_data["A_temp_471"]

A_concent_139 = mat_data["A_concent_139"]
A_concent_278 = mat_data["A_concent_278"]
A_concent_471 = mat_data["A_concent_471"]

A_adsorb_139 = mat_data["A_adsorb_139"]
A_adsorb_278 = mat_data["A_adsorb_278"]
A_adsorb_471 = mat_data["A_adsorb_471"]

# Define distance range (Ensure it matches the solution length)
distance = range(0, stop=10, length=60)

# Extract Julia simulation data at time points
sol_139 = sol(139 * 24 * 60 * 60)
sol_278 = sol(278 * 24 * 60 * 60)
sol_471 = sol(471 * 24 * 60 * 60)

# Plot Temperature
p1 = Plots.plot(distance, sol_139[T, 1:end-1], label="Current Model, t= 139th day", linewidth=2, color=:red)
 Plots.plot!(p1, distance, sol_278[T, 1:end-1], label="Current Model, t= 278th day", linewidth=2, color=:blue)
 Plots.plot!(p1, distance, sol_471[T, 1:end-1], label="Current Model, t= 471th day", linewidth=2, color=:black)

 Plots.plot!(p1, A_temp_139[:,1], A_temp_139[:,2], linestyle=:dash, linewidth=2, color=:magenta, label="Kumaran et. al (2019), t=139th day")
 Plots.plot!(p1, A_temp_278[:,1], A_temp_278[:,2], linestyle=:dash, linewidth=2, color=:blue, label="Kumaran et. al (2019), t=278th day")
 Plots.plot!(p1, A_temp_471[:,1], A_temp_471[:,2], linestyle=:dash, linewidth=2, color=:black, label="Kumaran et. al (2019), t=471th day")

Plots.xlabel!(p1, "Distance (m)")
Plots.ylabel!(p1, "Temperature (K)",guidefont=(6, "times"))
Plots.title!(p1, "Temperature Distribution in Coal Pile")

# Plot Oxygen Concentration
p2 = Plots.plot(distance, sol_139[c₁, 1:end-1], label="Current Model, t= 139th day", linewidth=2, color=:red)
 Plots.plot!(p2, distance, sol_278[c₁, 1:end-1], label="Current Model, t= 278th day", linewidth=2, color=:blue)
 Plots.plot!(p2, distance, sol_471[c₁, 1:end-1], label="Current Model, t= 471th day", linewidth=2, color=:black)

 Plots.plot!(p2, A_concent_139[:,1], A_concent_139[:,2], linestyle=:dash, linewidth=2, color=:magenta, label="Kumaran et. al (2019), t=139th day")
 Plots.plot!(p2, A_concent_278[:,1], A_concent_278[:,2], linestyle=:dash, linewidth=2, color=:blue, label="Kumaran et. al (2019), t=278th day")
 Plots.plot!(p2, A_concent_471[:,1], A_concent_471[:,2], linestyle=:dash, linewidth=2, color=:black, label="Kumaran et. al (2019), t=471th day")

Plots.xlabel!(p2, "Distance (m)")
Plots.ylabel!(p2, "Oxygen Concentration (#v/v)",guidefont=(6, "times"))
Plots.title!(p2, "Oxygen Concentration in Coal Pile")

# Plot Oxygen Adsorbed
p3 = Plots.plot(distance, sol_139[c₂, 1:end-1], label="Current Model, t= 139th day", linewidth=2, color=:red)
 Plots.plot!(p3, distance, sol_278[c₂, 1:end-1], label="Current Model, t= 278th day", linewidth=2, color=:blue)
 Plots.plot!(p3, distance, sol_471[c₂, 1:end-1], label="Current Model, t= 471th day", linewidth=2, color=:black)

 Plots.plot!(p3, A_adsorb_139[:,1], A_adsorb_139[:,2], linestyle=:dash, linewidth=2, color=:magenta, label="Kumaran et. al (2019), t=139th day")
 Plots.plot!(p3, A_adsorb_278[:,1], A_adsorb_278[:,2], linestyle=:dash, linewidth=2, color=:blue, label="Kumaran et. al (2019), t=278th day")
Plots.plot!(p3, A_adsorb_471[:,1], A_adsorb_471[:,2], linestyle=:dash, linewidth=2, color=:black, label="Kumaran et. al (2019), t=471th day")

Plots.xlabel!(p3, "Distance (m)")
Plots.ylabel!(p3, "Oxygen Adsorbed in Coal (kg O₂/kg coal)", guidefont=(6, "times"))
Plots.title!(p3, "Oxygen Adsorbed in Coal Pile")

# Combine all plots in a single window
Plots.plot(p1, p2, p3, layout=(3,1), size=(600, 600))

