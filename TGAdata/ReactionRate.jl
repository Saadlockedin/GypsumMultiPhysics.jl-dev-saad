using CSV
using DataFrames 
using MAT
using LinearAlgebra
using Plots

mat_file_path = "TGAdata/TGA_GW.mat"
mat_data = matread(mat_file_path)
println("Variables in the .mat file: ", keys(mat_data))

tga_data= mat_data["gyp_90_GWT2"]


t = tga_data[:,1]
T = tga_data[:,2].+273
m = tga_data[:,3]
m0=m[1]
minf=m[end]
α=abs.(m0.-m)/(m0-minf)
da_dT = zeros(length(t))

# Forward difference for the first point
da_dT[1] = (α[2] - α[1]) / (T[2] - T[1])

# Central difference for the interior points
for i in 2:length(t)-1
    da_dT[i] = (α[i+1] - α[i-1]) / (T[i+1] - T[i-1])
end

# Backward difference for the last point
da_dT[end] = (α[end] - α[end-1]) / (T[end] - T[end-1])
dal_dT=abs.(da_dT)

#=da_dT = zeros(size(t))
for i in 2:length(t) 
    da_dT[i] = (α[i] - α[i-1]) / (T[i] - T[i-1]) 
end
=# 
dT_dt = zeros(length(t))
# Forward difference for the first point
dT_dt[1] = (T[2] - T[1]) / (t[2] - t[1])
# Central difference for the interior points
for i in 2:length(t)-1
    dT_dt[i] = (T[i+1] - T[i-1]) / (t[i+1] - t[i-1])
end
# Backward difference for the last point
dT_dt[end] = (T[end] - T[end-1]) / (t[end] - t[end-1])
dT_dt
plot(t,T)
β=90/60 # Heating rate °C/s 
x=log.(β.*dal_dT)          
                      
inv_T=(1 ./T)
#plot((1 ./T),x)     
            
alpha_levels =[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#alpha_levels = range(0.1, stop=0.9, step=0.01)

# Arrays for storing the fitted parameters
E_values = Float64[]
A_values = Float64[]
n_values = Float64[]

# Perform linear regression for each conversion level α
for α_level in alpha_levels
#α_level = alpha_levels[2]
    # Find the indices where α is close to the desired level
    idx = findall(x -> isapprox(x, α_level, atol=0.05), α)

    if !isempty(idx)
        # Get corresponding data points for the fitting
        ln_dalpha_dt_fit = dal_dT[idx]
        inv_T_fit = inv_T[idx]
        alpha_fit = α[idx]
        
        # Prepare the design matrix X for linear regression: [1, 1/T, ln(α)]
        X = hcat(ones(length(ln_dalpha_dt_fit)), inv_T_fit, log.(1 .-alpha_fit))
        
        # Perform linear regression to solve X * [ln(A); -E/R; n] = ln(dα/dt)
        Y = alpha_fit
        
        # Solve for the parameters [ln(A), -E/R, n]
        params = X \ Y
        
        # Extract the parameters
        ln_A = params[1]
        E_R = -params[2]  # -E/R (because of the negative sign in the equation)
        n = params[3]
        
        # Calculate A from ln(A)
        A = exp(ln_A)
        E = E_R * 8.314  # Activation energy in J/mol, assuming R = 8.314 J/mol·K
        
        # Store the results
        push!(E_values, E)
        push!(A_values, A)
        push!(n_values, n)
    end
end

Activation_energy= E_values
Pre_exponential= A_values
Reaction_order=n_values

plot(alpha_levels, E_values, label="Activation Energy (E)", xlabel="α (Conversion)", ylabel="E (J/mol)", title="Activation Energy vs Conversion")
plot()
plot(A_values,E_values)


