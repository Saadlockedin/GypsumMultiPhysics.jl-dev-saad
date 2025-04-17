"""
Materials Properties Dictionary (`mp`)

This dictionary (`mp`) stores material properties for various materials. Each material is a key,
and its properties are stored in nested dictionaries.

"σₑ"
"σᵢ"

## Abbreviations
- `"σₑ"`: Electrical conductivity
- `"k"`: Thermal conductivity
- `"Cp"`: Specific heat capacity
- `"α"`: Thermal expansion coefficient
- `"ρ"`: Density
- `"σᵢ"`: Ionic conductivity

## Example Usage
To access and compute a property for a material at a specific temperature `T`, use:
```julia
map(mp["Ni"]["σₑ"], 1000)  # Compute electrical conductivity for Ni at T = 1000 K
"""
mp = Dict()

""" Properties of Nickel (Ni):

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, depends on T
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 
"sϵ": solid phase volume fraction 

"""
varschem = matread("D:\\GSU\\M ENG\\Research\\Gypsum\\GypsumMultiPhysics.jl\\MATLABcodes\\ChemistryDatabase_new.mat");
function extract(d)
    v = collect(values(d));
    k = collect(keys(d));
    expr = quote end
    for (k, v) in d
        push!(expr.args, :($(Symbol(k)) = $v))
    end
    eval(expr)
    return
end
extract(varschem);

interpk_gyp      = LinearInterpolation(k_gyp[:,1] .+ 273.0,k_gyp[:,2]);

LatentHeat = sort!(LatentHeat,dims = 1); 
interpLatentHeat = LinearInterpolation(LatentHeat[:,1],LatentHeat[:,2]);

Heat_Flow = sort!(Heat_Flow, dims = 1);
interpHeat_Flow  = LinearInterpolation(Heat_Flow[:,1] .+ 273.0,Heat_Flow[:,2]);

CpEff_gyp = sort!(CpEff_gyp,dims = 1); 
interpCpEff_gyp = LinearInterpolation(CpEff_gyp[:,1] .+ 273.0,CpEff_gyp[:,2]); # the interpolation function of the effective CpEff

interpMg  = LinearInterpolation(Mg[:,1],Mg[:,2]*1e-5);
interpMgv = LinearInterpolation(Mg[:,1],Mg[:,3]*1e-5);

interpCpVapor = LinearInterpolation(CpVapor[:,1],CpVapor[:,2]); # Specific heat capacity of water steam [J/kg.K].
interpkwater  = LinearInterpolation(kwater[:,1],kwater[:,2]);
interpkair    = LinearInterpolation(kair[:,1],kair[:,2]);

interpWT12 = LinearInterpolation(waterTable[:,1] .+ 273.0,waterTable[:,2]*1000.0); 
interpWT21 = LinearInterpolation(waterTable[:,2]/100.0,waterTable[:,1] .+ 273.0); 
interpWT13 = LinearInterpolation(waterTable[:,1] .+ 273,waterTable[:,3]); #,BP[i,j,k]]

mp["Gyp"] = Dict(
    "LH" => (T) -> interpLatentHeat(T), #Latent Heat
    "k" => (T) -> interpk_gyp(T), #Thermal conductivity
    "HF" => (T) -> interpHeat_Flow(T), #Heat Flow
    "Cp" => (T) -> interpCpEff_gyp(T), #Cp
    "MugA" =>  (T) -> interpMg(T),
    "MugV" =>  (T) -> interpMgv(T),
    "CpV" =>  (T) -> interpCpVapor(T), # Specific heat capacity of water steam [J/kg.K]
    "kw" =>  (T) -> interpkwater(T),
    "ka" =>  (T) -> interpkair(T),
    "MWair" => 28.97, # gram/mole molecular weight of water 
    "MWh2o" => 18.02, # gram/molde molecular weight of air
    "Ru"    => 83.144598, # bar.cm3/(mol.K)  Gas constant
    "Pamb"  => 1.01325,   # bar Pressure
    "Pvapor"    => 0.0093,           # Water initial pressure at ambient x = 0; 0.4*17.5*0.00133322
    "Yv_E"      => 0.1239,                           # Mass fraction of water vapor at fire side
    "Xv"      => 0.0092,                  # ambient water mole fraction at ambient side .0093/1.01325
    "MWmix"     => 28.8691,      # Molecular weight of air-water mixture  Xv*MWh2o + (1-Xv)*MWair
    "Yv_W"      =>  0.0057,                 # Mass fraction of water vapor Xv*MWh2o./MWmix
    "phi_vapor" => 0.0067, #Initial density of water in the computational domain ; [kg/m3] 1000.*Yv_W.*Pamb.*MWmix./(Ru.*Temp_W);
    "phi_air" => 1.1583, # Initial density of air in the computational domain 1000*(1-Yv_W).*Pamb.*MWmix/(Ru.*Temp_W);
    "hm" => 0.00955, # [m/s] convective mass transfer coefficient for species ransport [m/s]
    "hcE" => 25, # [W/m2.K]
    "hcW" => 9, # [W/m2.K]
    "Φ"=> 0.85, # Porosity
    "τ"=> 1.869, #tortuosity
    "D_ab"=> 2.56e-5, # [m2/s] Binary diffusion coefficient 
    "D_eff"=> 0.1164, # Φ*D_ab/τ 
    "Uᴰ"=> 0.6e-5, #darcy velocity 
    "μ"=> 1.862e-05, # dynamic viscosity of air
    "K"=> 1.00e-15, #Permeability tensor"
    "ϵ" => 0.85, # emissivity
    "σ" => 5.67e-8, # stefan boltzmann constant
    "Tgas"=> 303, # Tgas is the mean temperature at the gas phase on the fire side
    "rho_in"=> 810, # initial density, kg/m3
    "H20_unbound_frac" => 0.02, # unbound water present
    "H20_bound_frac" => 0.20927, # bound water present
    "Mug" => 0.1854e-4,
    "kgass" => 0.0265, # ka'.*phi_air./(phi_vapor+phi_air) + kw'.*phi_vapor./(phi_vapor+phi_air)
    "alphaT" => 0.7973, # (interpk_gyp.*(1-Φ)+kgass*Φ)./RhoCpeff;
    "RhoCpeff" => 8.11e5, # rho_in'.*1000 + Φ.*(phi_vapor.*CpV+phi_air.*CpA) ; CpA=1000
    "H2o_initial" => 185.7087, #(H20_unbound_frac + H20_bound_frac)*rho_in
    "H2o_st_1" => 143.3315, #(H20_unbound_frac + 0.75*H20_bound_frac)*rho_in
    
    )

mp["Ni"] = Dict(
    "σₑ" => (T) -> 9.5e7 / T * exp(-1150.0 / T), 
    "k" => (T) -> 1.380217e4 * T^(-0.927164) * exp(1.145850e-3 * T - 33.13548 / T),
    "Cp" => (T) -> 6.59e-5 * T^2 - 2.16e-2 * T + 506.0,
    "α" => (T) -> 4.4e-9 * T + 11.01e-6,
    "ρ" => 8900.0,
    "sϵ" => 0.5
)

""" Properties of AISI 441:

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, depends on T
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 

"""
mp["AISI 441"] = Dict(
    "σₑ" => (T) -> 0.4862 * T^2 - 1356.70 * T + 1.77e6,
    "k" => (T) -> 7.0142e-6 * T^2 - 7.3196e-3 * T + 25.477,
    "Cp" => (T) -> -2.82e6 / (T^2) + 0.136 * T + 472.0,
    "α" => (T) -> 3.22e-9 * T + 9.94e-6,
    "ρ" => 5960.0
)

""" Properties of Alloy-718:

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, depends on T
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 

"""
mp["Alloy-718"] = Dict(
    "σₑ" => (T) -> -1.13e-3 * T^3 + 4.2061 * T^2 - 5.257e3 * T + 2.914e6,
    "k" => (T) -> 7.72e-6 * T^2 + 7.35e-3 * T + 8.87,
    "Cp" => (T) -> 0.0937 * T + 515.4,
    "α" => (T) -> 3.73e-9 * T + 10.34e-6,
    "ρ" => 8190.0
)

""" Properties of Pt:

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, depends on T
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 

"""
mp["Pt"] = Dict(
    "σₑ" => (T) -> 2.5772 * T^2 - 7.662e3 * T + 7.875e6,
    "k" => (T) -> 19.24495 * T^0.188281 * exp(8.41185e-5 * T + 52.08891 / T),
    "Cp" => (T) -> 0.0134 * T + 131.41,
    "α" => (T) -> 1.35e-15 * T^3 - 3.21e-12 * T^2 + 4.94e-9 * T + 7.68e-6,
    "ρ" => 21452.0
)

""" Properties of Au:

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, constant
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 

"""
mp["Au"] = Dict(
    "σₑ" => (T) -> 1.41e10 / T - 2.956e6,
    "k" => (T) -> 49.21633 * T^0.328633 * exp(-5.871e-4 * T + 47.35758 / T),
    "Cp" => 125.6,
    "α" => (T) -> 4.25e-9 * T + 11.45e-6,
    "ρ" => 19300.0
)

# Lanthanum Strontium Manganite (LSM)
""" Properties of Lanthanum Strontium Manganite (LSM):

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, constant
"Cp": Specific heat capacity, constant
"α": Thermal expansion coefficient, constant
"ρ": Density (constant) 
"σᵢ": ionic conductivity, depends on T
"sϵ": solid phase volume fraction

"""

mp["LSM"] = Dict(
    "σₑ" => (T) -> 3.48e7 / T * exp(-0.1 / (k_B * T)),
    "k" => 4.0,
    "Cp" => 573.0,
    "α" => 12.0e-6,
    "ρ" => 5120.0,
    "σᵢ" => (T) -> 1.508e4 / T * exp(-0.5 / (k_B * T)),
    "sϵ" => 0.5
)

# Yttria-Stabilized Zirconia (YSZ)
""" Properties of Yittria Stabilized Zirconia (YSZ):

"σₑ": Electrical conductivity, constant
"k": Thermal conductivity, constant
"Cp": Specific heat capacity, depends on T 
"α": Thermal expansion coefficient, constant
"ρ": Density (constant) 
"σᵢ": ionic conductivity, depends on T
"sϵ": solid phase volume fraction

"""
mp["YSZ"] = Dict(
    "σₑ" => 1e-12,
    "k" => 2.0,
    "Cp" => (T) -> -5.0e-5 * T^2 + 0.185 * T + 476.0,
    "α" => 10.3e-6,
    "ρ" => 5960.0,
    "σᵢ" => (T) -> 1.0033e7 / T * exp(-0.85 / (k_B * T)),
    "sϵ" => 0.5
)

# Mica
""" Properties of Au:

"σₑ": Electrical conductivity, constant
"k": Thermal conductivity, constant
"Cp": Specific heat capacity, constant
"α": Thermal expansion coefficient, constant

"""
mp["Mica"] = Dict(
    "σₑ" => 3000.0,
    "k" => 0.52,
    "Cp" => 770.0,
    "α" => 10.3e-6
)

# Gadolinium-Doped Ceria (GDC)
""" Properties of YGadolinium-Doped Ceria (GDC):

"σₑ": Electrical conductivity, constant
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, depends on T 
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 
"σᵢ": ionic conductivity, depends on T
"sϵ": solid phase volume fraction

"""
mp["GDC"] = Dict(
    "σₑ" => (T) -> 4.6e-1,
    "k" => (T) -> 2e-5 * T^2 - 0.0259 * T + 11.79,
    "Cp" => (T) -> 2.6957 * T - 108.7,
    "α" => (T) -> -5e-16 * T^4 + 1e-12 * T^3 - 1e-9 * T^2 + 4e-7 * T - 6e-5,
    "ρ" => 7000.0,
    "σᵢ" => (T) -> 0.0186 * T - 8.685,
    "sϵ" => 1.0 
)

# Lanthanum Strontium Cobalt Ferrite (LSCF)
""" Properties of Lanthanum Strontium Cobalt Ferrite (LSCF):

"σₑ": Electrical conductivity, depends on T
"k": Thermal conductivity, depends on T
"Cp": Specific heat capacity, constant 
"α": Thermal expansion coefficient, depends on T
"ρ": Density (constant) 
"σᵢ": ionic conductivity, depends on T
"sϵ": solid phase volume fraction


"""
mp["LSCF"] = Dict(
    "σₑ" => (T) -> 6.3889 * T - 777.78,
    "k" => (T) -> 5e-6 * T^2 - 0.0083 * T + 6.15,
    "Cp" => 120.0,
    "α" => (T) -> 2.2727 * T - 8e-6,
    "ρ" => 2100.0,
    "σᵢ" => (T) -> 1.3983 * T + 0.8602,
    "sϵ" => 1.0 
)


# Define the material assignment function
function assign_materials!(electrode::GypsumMultiPhysics.AutoUpdateDict{String, Any}, materials::Vector{String}, materials_dict::Dict{Any, Any})
    # Assign the material names to "mat"
    electrode["mat"] = materials
    
    # Copy properties for each material
    for material in materials
        if haskey(materials_dict, material)
            for (prop, value) in materials_dict[material]
                # Store each material's property under a unique key
                electrode["$material.$prop"] = value
            end
        else
            error("Material $material not found in the provided materials dictionary.")
        end
    end
end



