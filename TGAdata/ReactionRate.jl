using CSV
using DataFrames 
using MAT
using Plots

mat_file_path = "TGAdata/TGA_GW.mat"
mat_data = matread(mat_file_path)
println("Variables in the .mat file: ", keys(mat_data))

tga_data= mat_data["gyp_90_GWT2"]


t = tga_data[:,1]
T = tga_data[:,2]
m = tga_data[:,3]
m0=m[1]
minf=m[end]
α=(m0.-m)/(m0-minf)
da_dT = zeros(size(t))
for i in 2:1138
    da_dT[i] = (α[i] - α[i-1]) / (T[i] - T[i-1])
end
#da_dT
da_dT=da_dT.+ abs(minimum(da_dT)) .+ 1e-10
β=10/60 # Heating rate °C/s
x=log.(β*da_dT)

plot((1/T)',x)