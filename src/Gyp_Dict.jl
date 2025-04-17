"""
# Solid Oxide Fuel Cell (SOFC) / Solid Oxide Electrolysis Cell (SOEC) Modeling

This script provides data structures and functions for defining the electrode and electrolyte 
assemblies in SOFC/SOEC systems. It allows the user to specify:
- Initial species concentrations and cell component thicknesses.
- Customizable species in the fuel channel (e.g., H₂, CH₄, H₂O, N₂, Ar, CO, CO₂).
- Species in the air channel (e.g., O₂, N₂).
- Layered structures for anode, cathode, and electrolyte.
- Size of the components.

## Main Features
1. Species indexing and identification for fuel, air, and electrolyte components.
2. Grid creation for discretizing spatial domains.
3. Microstructural property definitions for porosity, tortuosity, pore diameter, and particle diameter.

### Required Packages
- `LessUnitful` for unit handling.
"""

@unitfactors mol bar atm R cm mm m s Pa K


#k_B = k;
k_B = 8.617333262145e-5 #eV/K
tr_file = "lib/transport.dat"
therm_file = "lib/therm.dat"
#=
const F = N_A * e # Faraday constant
const neH2 = 2; # number of electrons for hydrogen oxidation
const neO2 = 4; 
const neCO = 2;
=#
const Ru = R;


#const T = 800+273K
#const P = 1atm


mutable struct AutoUpdateDict{K, V} <: AbstractDict{K, V}
    data::Dict{K, V}
    update_funcs::Dict{K, Function}  # Dictionary of update functions keyed by dependent field

    # Constructor
    function AutoUpdateDict{K, V}() where {K, V}
        new{K, V}(Dict{K, V}(), Dict{K, Function}())
    end
end

# Implement dictionary interface
Base.getindex(d::AutoUpdateDict, key) = get(d.data, key, nothing)

Base.setindex!(d::AutoUpdateDict, value, key) = begin
    d.data[key] = value
    if haskey(d.update_funcs, key)
        d.update_funcs[key]()  # Call the update function if it exists
    end
end

Base.length(d::AutoUpdateDict) = length(d.data)
Base.keys(d::AutoUpdateDict) = keys(d.data)
Base.pairs(d::AutoUpdateDict) = pairs(d.data)
Base.iterate(d::AutoUpdateDict, state...) = iterate(d.data, state...)



    #=
    ## Microstructural data 
    - porosity
    - tortuosity
    - pore diameter
    - particle diameter
    pm = Properties(ϵ, τ, pore_dia, part_dia)

    FE["pm"] is the dictionary that store the microstructural data in fuel electrode
    AE["pm"] is the dictionary that store the microstructural data in air electrode
    pmL is for the lower bound values of microstructural data
    pmH is for the higher bound values of microstructural data
    =#
"""
create_FE(; dist_func)

Create and initialize the fuel electrode (FE) structure for SOFC/SOEC systems.

# Arguments:
- `dist_func`: A distribution function (optional, default: η -> 0.5 * (1 + tanh(10 * (η - 1)))).

# Returns:
- `FE`: A dictionary containing fuel electrode properties, species data, grid, and microstructural properties.
- `SpecIndFE`: A dictionary mapping species names to their indices in the FE region.
"""

    function create_FE(; dist_func = η -> 0.5 * (1 + tanh(10 * (η - 1))))
        FE = AutoUpdateDict{String, Any}()
        
        # Define default values
        FE["Sp"] = ["H2O", "N2", "O2"]
        #FE["tData"] = create_transport_data(FE["Sp"], tr_file)
        #FE["tObj"] = create_thermo(FE["Sp"], therm_file)
        FE["ImF"] = [0.0092, 0.772824, 0.217976] # initial mole fraction put formula instead of values 
        #FE["mF"]  = copy(FE["ImF"])  # mole fractions 
        FE["slt"] = 12mm #surface layer thickness
        #FE["alt"] = 100μm
        #FE["pmL"] = Properties(0.30, 1.5)
        #FE["pmH"] = Properties(0.40, 5.0)
        FE["dxL"] = 1e-4 # lowest grid size
        FE["dxH"] = 1e-4 # highest grid size, if dxL = dxH a uniform grid will be generated. 
        FE["mat"] = ["Gyp"] #Layers
        #FE["matvf"]  = [0.5, 0.5] # the sum of the solid phases volume fraction should be 1.  
        #FE["Niϵ"] = 0.5 # volume fraction of solid phase of Nickel 
        FE["T"]   = 1273K # initial temperature inside the fuel electrode 
        FE["P"]   = 1atm # initial pressure inside the fuel electrode

        assign_materials!(FE, FE["mat"], mp)
        

        #= Grid and microstructure update function
        function update_grid_and_microstructure()
            FE["X"] = geomspace(0, FE["alt"] + FE["slt"], FE["dxH"], FE["dxL"])
            grid = simplexgrid(FE["X"])
            FEGrid = grid[Coordinates]/FE["alt"] 
            FEGrid = -(FEGrid .- FEGrid[end])
            FEfdist = map(dist_func, FEGrid)    
            # Update microstructural properties
            FE["ϵ"] = FE["pmL"].ϵ * (1 .- FEfdist) .+ FE["pmH"].ϵ .* FEfdist
            FE["τ"] = FE["pmL"].τ * (1 .- FEfdist) .+ FE["pmH"].τ .* FEfdist
            FE["pore_dia"] = FE["pmL"].pore_dia * (1 .- FEfdist) .+ FE["pmH"].pore_dia .* FEfdist
            FE["part_dia"] = FE["pmL"].part_dia * (1 .- FEfdist) .+ FE["pmH"].part_dia .* FEfdist

            FE["lTPB"]  = similar(FE["ϵ"]) .+ 5.2147e+12
            FE["Cdl"]   = similar(FE["ϵ"]) .+ 7.919
            # Carman-Kozeny Equation is used for permiability
            FE["K"] = ((FE["ϵ"].^3)./(1 .-FE["ϵ"]).^2).*(FE["part_dia"].^2)./FE["τ"].^2 # permiability

            FE["σₑ"] = (1 .- FE["ϵ"])*FE["Ni.σₑ"](FE["T"])*FE["Ni.sϵ"]
            FE["σᵢ"] = (1 .- FE["ϵ"])*FE["YSZ.σᵢ"](FE["T"])*(1-FE["Ni.sϵ"])

         #   FE["σₑ"] = (1 .- FE["ϵ"])*FE["Ni.σₑ"](T)*FE["Ni.sϵ"]
         #   FE["σᵢ"] = (1 .- FE["ϵ"])*FE["YSZ.σᵢ"](T)*(1-FE["Ni.se"])

            FE["grid"] = grid
            FE["tData"] = create_transport_data(FE["Sp"], tr_file)
            FE["tObj"] = create_thermo(FE["Sp"], therm_file)

            # Create dictionary to map species to index
            spec = vcat(FE["Sp"], ["ϕₑ", "ϕᵢ"])
            FE["SpecInd"] = Dict(zip(spec, 1:length(spec)))
        end
    
        
    
        FE.update_funcs["alt"] = update_grid_and_microstructure
        FE.update_funcs["slt"] = update_grid_and_microstructure
        FE.update_funcs["pmL"] = update_grid_and_microstructure
        FE.update_funcs["pmH"] = update_grid_and_microstructure
        FE.update_funcs["dxL"] = update_grid_and_microstructure
        FE.update_funcs["dxH"] = update_grid_and_microstructure
        FE.update_funcs["Sp"] = update_grid_and_microstructure
    
        # Initial grid and microstructure setup
        update_grid_and_microstructure()

       
    =#
        return FE
    end