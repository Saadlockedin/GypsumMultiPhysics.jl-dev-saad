module GypsumMultiPhysics

using VoronoiFVM
using ExtendableGrids
using AMGCLWrap
using ExtendableSparse
using GridVisualize: gridplot, gridplot!
using DiffusionFlux
using IdealGas
using TransportProperties
using RxnHelperUtils
using LinearAlgebra
using DataFrames
using DataStructures
using LinearSolve
using StaticArrays
using OrdinaryDiffEq
using LessUnitful
using MAT
using Interpolations

# -------------------------------------------------------------------------
# Flux
# -------------------------------------------------------------------------
# Nernst potential calculations
include("FicksFlux.jl")
export Flux_Ficks_Array_Darcy!
export Flux_Ficks_Array!
include("Gyp_Dict.jl")
export create_FE
include("Gyp_mp.jl")
export mp, assign_materials!


end