
#--------------------------------------------------------------------------
# Numerical solution of the non-dimensional steady convection-diffusion
# equation for heat & mass transfer
#
#       PHI,t, + (u*PHI),x = sig*PHI,xx + C1*PHI + C2*PHI^2 + Sphi
#
#       Calcination inside the gypson board
#                   
#                               Hayri Sezer
#                  Georgia Southern University 02/22/2022
#--------------------------------------------------------------------------

## add used libraries here 

using Plots; # use plotting in julia library
gr(); # this is a backend used in julia plotting library
using JFVM;
#using JFVMvis; # finite volume model developed by simukalde 
using MAT; # used to import matlab data
using Interpolations; # interpolarion library from Julia 
using XLSX;
using LaTeXStrings;
using Parameters;
using Trapz

function gypsumfvm3d()
## time parameters
dt       = .1
# time_max = 60*24+20; # min()
time_max = 400; 
savetime = 2; # save data at each 30 seconds
urelax   = 0.1;
## molecular weight of water and air; Pressure & Gas constant
MWair = 28.97; # gram/mole
MWh2o = 18.02; # gram/molde
Ru    = 83.144598; # bar.cm3/(mol.K) :
Pamb  = 1.01325;   # bar()
## generate the grid()
nx      = 60;
ny      = 20; 
nz      = 20; 
Lx = 0.012; # [m] thickness of gympson board
Ly      = 0.6;  # [m]
Lz      = 0.6;  # [m]
Volume = Lx*Ly*Lz; # [m3]
m = createMesh3D(nx,ny,nz,Lx,Ly,Lz); # create grid
## Upload the experimental temperature from the fire side 
Tfire_exp = matread("Tfire_exp.mat");
# load Tfire_exp; # experimentally measured temperature at fire side.
Tfire = Tfire_exp["Tfire_exp"];
#load hicks_c12_t2
#hicks = hicks_c12_t2
# load TFire2
# Tfire = TFire2
## use the interpolation to get the temperature corresponding to a specific time 

#Temp_E = interp1q[Tfire[:,1]*60,Tfire[:,2],0]+273;# Fire side    [K]
interpTfire = LinearInterpolation(Tfire[:,1]*60, Tfire[:,2]);
Temp_E = interpTfire(0)+273.0; 


#Temp_E = interp1q[hicks[:,1],hicks[:,2],0];# Fire side    [K]
Temp_W = 29.0+273.0; #interp1q[Tfire[:,1]*60,Tfire[:,2],0]+273; # Ambient side [K]
## initial condition
BCTemp = createBC(m); # all the boundaries are flux zero by default. 
BCP = createBC(m); 
Temp = createCellVariable(m,Temp_W); 
Temp_old = Temp;
## initiliaze the drag velocity
U_D = arithmeticMean(createCellVariable(m,0.0));
# U_Dx = arithmeticMean[U_D]; 
## EKU experimental data for model comparison
fname1  = "70kw_flametrial_9min.xlsx"; 
xf = XLSX.readxlsx("70kw_flametrial_9min.xlsx");
sh = xf["Untitled"];
shh = sh[:];
timeR = shh[2:end,2];
termo_couples = shh[2:end,3:2:27];
#timeR = data[:,1]
#termo_couples = data[:,2:2:26]
## ---- Here we should cut the initial data from the experiments. Because in the experiments the heater was not on for initial time of 2 minutes approximately. 
cutinT = 145; 
timeRR = timeR[findall(timeR.>cutinT)]; 
termo_couples = termo_couples[findall(timeR.>cutinT),:]; 
timeRR = timeRR.-timeRR[1]; 

termo_couples = convert(Array{Float64}, termo_couples);
timeR = timeRR; 
#
#
# Create figure()
plot(timeR,termo_couples)

plot(timeR,termo_couples, xaxis=("time [s]"), yaxis=(L"Temperature [^oC]"))

## initial condition for species transport
P = createCellVariable(m,Pamb); 
P_old     = P;
Pvapor    = 0.4*17.5*0.00133322;               # Water initial pressure at ambient x = 0
Yv_E      = 0.1239;                            # Mass fraction of water vapor at fire side
Xv_W      = Pvapor./P.value;                         # ambient water mole fraction at ambient side
Xv        = Xv_W;
MWmix     = Xv_W.*MWh2o + (1 .-Xv_W).*MWair;       # Molecular weight of air-water mixture
Yv_W      = Xv_W.*MWh2o./MWmix;                  # Mass fraction of water vapor
## Here initialize the cell variables for the vapor
BCphiV  = createBC(m); # all the boundaries are flux zero by default. 
BCphiA  = createBC(m);
phi_vapor = createCellVariable(m,0);     # Initial density of water in the computational domain ; [kg/m3]
phi_vapor.value .= 1000.0 .*Yv_W.*Pamb.*MWmix./(Ru.*Temp_W);
phi_air   = createCellVariable(m,0); # Initial density of air in the computational domain
phi_air.value .= 1000.0 .*(1 .-Yv_W).*Pamb.*MWmix/(Ru.*Temp_W);
phi_vapor_old = phi_vapor;
phi_air_old   = phi_air; 

## End of cell variables
hm        = 9.55e-3; # [m/s] convective mass transfer coefficient for species ransport [m/s]
hcE        = 25;    # [W/m2.K]
hcW        = 9;    # [W/m2.K]
Pore      = 0.85;  # porosity
K         = 1e-15; # Permiability [m2]
DAB       = 2.56e-5; # [m2/s] Binary diffusion coefficient
Tau       = 1.869; # tortuosity
Deff     = createCellVariable(m,Pore.*DAB./Tau); 
DeffFace = arithmeticMean(Deff); # diffusion should be a face variable. 
# Deff[1:nx]= Pore*DAB/Tau; # effective diffusion of gases.
e         = 0.85; # emissivity
sigma     = 5.67*10^(-8); # stefan boltzmann constant
Tgas      = 303.0; # Tgas is the mean temperature at the gas phase on the fire side.
rho_in = 810.0; #initial density; kg/m3
H20_unbound_frac = 0.02; #unbound water present
H20_bound_frac = 0.20927;
## Here upload the data for the reaction rates & thermo-physical properties
# data is obtained from
#
#  K. Ghazi Wakili; E. Hugi; L. Wullschleger; T. Frank; Gypsum board in fire e
#   modeling & experimental validation. Journal of Fire Sciences 25 [2007]
#
# Data is constructed by linear interpolation
#
#--------------------------------------------------------------------------

## load data from data base
# load database_gyp  # --- old data base 

#load ChemistryDatabase.mat    # --- New data base 
# This code is used to extract the dictionaries as global variables. 
# Any .mat file can be extracted in Julia. 

varschem  = matread("ChemistryDatabase_new.mat");
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


# Here extract data from matlab .mat file as a dictionary in Julia. Julia dictionary can be converted to variables 

#= 
--------------------------------------------------------------------------------------------------
rho_gyp   = varschem["rho_gyp"];
CpEff_gyp = varschem["CpEff_gyp"];
Heat_Flow = varschem["Heat_Flow"];
k_gyp     = varschem["k_gyp"];
W_loss    = varschem["W_loss"];
-----------------------------------------------------------------------------------------------------
=#

# load LatentHeat
# load rho_gyp; # load the density data, rho_gyp [kg/m3], T [C]
# load CpEff_gyp; # load the specific heat data; CpEff_gyp [J/Kg.K], T [C]
# load Heat_Flow; # load the heat flow endothermic reaction of gypson board; HF [mW], T [C]
# load k_gyp; # load the thermal conductivity data, k_gyp [W/m.K], T [C]
# load W_loss; # load weight loss: W_loss = [#loss]
# load W_Loss_rate; # weight loss rate percentage.
# load HF_Lid05mm; # load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry [DSC] of pure CaSO4.2H2O
# load HF_Lid1mm; # load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry [DSC] of pure CaSO4.2H2O
# load HF_NoLid; # load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry [DSC] of pure CaSO4.2H2O
# load Mg; # load data for the dynamic viscosity of air
# load MLR
# load kgas
# # load MLRnew
# load CpVapor
## ------------------------------------------------------------------------
HeatingRate = 50; # 20 degree/min is converted to degree/second 
## ------------------------------------------------------------------------
str1 = string(HeatingRate);
fname = "HeatingRate"*str1*".xlsx";

str2 = "W_loss = W_loss_White_";
str3 = str2*str1*";";
eval(str3); 


# W_loss = W_loss_White_20[:,[1,2]] 
# W_Loss_rate = W_loss_White_20[:,[1,3]]; 
W_loss[:,2] = 100*W_loss[:,2]/W_loss[1,2];
rho_gyp = W_loss; 
rho_gyp[:,2] = rho_in*W_loss[:,2]/100;
Tempdata = rho_gyp[:,1].+273;
# 
# p = polyfit(Tempdata,rho_gyp[:,2],10); 
# 
# rho_gypFit[:,2] = smooth[Tempdata,rho_gyp[:,2],0.1,"rloess"]
# 
# # rho_gypFit[:,2] = polyval(p,Tempdata); 
rhodiff    = -diff(rho_gyp[:,2])./diff(Tempdata)*HeatingRate/60;
# rhodiffFit = -diff(rho_gypFit[:,2])./diff(Tempdata)*HeatingRate/60
# 
# figure()
# plot(Tempdata[1:end-1], rhodiff,'o')
# hold on 
# plot(Tempdata[1:end-1], rhodiffFit,'-')
# figure() 
# plot(Tempdata, rho_gyp[:,2],'o')
# hold on 
# plot(Tempdata, rho_gypFit[:,2],'-')

# rho_gyp[:,2] = rho_gypFit[:,2]; # uncomment this line when smoothing is needed
## calculate the effective latent heat of gypson board dehydration.

L1ind = findall(CpEff_gyp[:,1] .>100);
L1ind1 = findall(CpEff_gyp[:,1] .< 250);
L1ind = intersect(L1ind,L1ind1);
L1 = 5.6e+5; #

#--- Numerical integration 
LatentHeat = zeros(Float64, size(CpEff_gyp,1),2); 

L1 = (trapz(CpEff_gyp[L1ind,1],CpEff_gyp[L1ind,2])-CpEff_gyp[1,2]*150);
LatentHeat[1:size(CpEff_gyp,1),1] = CpEff_gyp[:,1] .+ 273.0;
LatentHeat[1:L1ind[1]-1,2] .= 1e+25;
LatentHeat[L1ind,2] .= L1;

#L2ind = find(CpEff_gyp(:,1) >600 & CpEff_gyp(:,1) <780);
L2ind = findall(CpEff_gyp[:,1] .>600);
L2ind2 = findall(CpEff_gyp[:,1] .< 780);
L2ind = intersect(L2ind,L2ind2);
L2 = (trapz(CpEff_gyp[L2ind,1],CpEff_gyp[L2ind,2])-CpEff_gyp[1,2]*180);

LatentHeat[L1ind[end]+1:L2ind[1]-1,2] .= 1e+25;
LatentHeat[L2ind,2] .= L2;
LatentHeat[L2ind[end]+1:end,2] .= 1e+25;





# inds = find(CpEff_gyp[:,1] <151)
# LatentHeat[inds,2] = 1e+25

# plot(LatentHeat[:,1],LatentHeat[:,2])
## here calculate the MLR
# for i = 2:length(rho_gyp)-1
#     
#     MLRnew[i,2] = -(rho_gyp[i+1,2]-rho_gyp[i-1,2])/(rho_gyp[i+1,1]-rho_gyp[i-1,1])/3
#     
# end
# MLRnew[1,2] = 0
# MLRnew[2] = 0
# MLRnew[length(rho_gyp),2] = MLRnew[length(rho_gyp)-1,2]
# 
# MLRnew[:,1] = rho_gyp[:,1]+273


# plot(MLRnew[:,1],MLRnew[:,2])


##

plot(LatentHeat[:,1],LatentHeat[:,2])

# plot(W_Loss_rate[:,1],W_Loss_rate[:,2],"-o")

interprho_gyp = LinearInterpolation(rho_gyp[:,1] .+ 273,rho_gyp[:,2]);
weight = Volume*interprho_gyp(Temp_W);
volume = 10.7e-6/interprho_gyp(Temp_W);
# W_loss_rate[:,2] = 0.01*(W_Loss_rate[:,2]).*interp1q[rho_gyp[:,1]+273,rho_gyp[:,2],Temp_W]; # this is the source term for vapor


# plot(W_Loss_rate[:,1], W_loss_rate[:,2])
# W_loss_rate[:,2] = 0.01*(W_Loss_rate[:,2]).*63.6173; # this is the source term for vapor
# W_loss_rate[:,1] = W_Loss_rate[:,1]+273

# plot(CpEff_gyp[:,1],CpEff_gyp[:,2])
# plot(k_gyp[:,1],k_gyp[:,2])
# plot(rho_gyp[:,1],rho_gyp[:,2])
# plot(Heat_Flow[:,1],Heat_Flow[:,2])
# plot(HF_NoLid[:,1]+273,HF_NoLid[:,2])
## Here create the cell variables 
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





    CpEffgyp  = createCellVariable(m,interpCpEff_gyp(Temp.value[10,10,10])); # Temperature values are converted to Kelvin.
    HF        = createCellVariable(m,interpHeat_Flow(Temp.value[10,10,10]));# Temperature values are converted to Kelvin. & mW is converted to W.
    kgyp      = createCellVariable(m,interpk_gyp(Temp.value[10,10,10]));
    kgypFace  = arithmeticMean(kgyp); 
    Mug       = createCellVariable(m,interpMg(Temp.value[10,10,10]));
    CpV       = createCellVariable(m,1.996*1000); # Specific heat capacity of water steam [J/kg.K].
    CpA       = createCellVariable(m,1000);     # Specific heat capcity of air [J/kg.K].
    kw        = createCellVariable(m,0.1);
    ka        = createCellVariable(m,0.1);
    MugA      = createCellVariable(m,0.1);
    MugV      = createCellVariable(m,0.1);
    LHeat     = createCellVariable(m,0.1);
    SphiTSource = createCellVariable(m,0.1);
    SphiVSource = createCellVariable(m,0.1);
    SphiASource = createCellVariable(m,0.1);
    MugFaceCell = createCellVariable(m,0.1);
    alphaT      = createCellVariable(m,0.00001);
    RhoCpeff    = createCellVariable(m,100);
    P           = createCellVariable(m,1.0132463948867976);
    LowHF = 10000.0; # kw of incident low heat flux; 
    HighHF = 40000.0; # kw of incident high heat flux
    
    HFl = zeros(1,ny,nz); 
    for i = 1:nz
        HFl[1,1:ny,i] = range(LowHF,HighHF,length=ny)
    end    
    incident_HFlux = HFl; # uniform incident_HFlux
    
## end of the thermo-physical data
## Reaction rate coefficients
## 
# range1 = [310,415]" ## [K]" temperature interval where we perform the fitting
# Tind = range1[1]<Tempdata & Tempdata<range1[2]
# [C_low,TableData1] = Reaction_Rate[rho_gyp,rhodiff,Tind,Tempdata,1]
# C_low = C_low'
# 
# range2 = [700,1000]'
# Tind = range2[1]<Tempdata & Tempdata<range2[2]
# [C_up,TableData2] = Reaction_Rate[rho_gyp,rhodiff,Tind,Tempdata,2]
# C_up = C_up'
# str4 = "coeff_"; 
# str4 = [str4,str1]; 
# table = table[range1,range2,C_low,C_up]; 
# 
# # Table = join[TableData1 , TableData2]
# writetable[TableData1,fname,"Sheet','data','Range','A1"]
# writetable[TableData2,fname,"Sheet','data','Range','D1"]
# writetable[table,fname,"Sheet',1,'Range','G1"]
# save(str4, "C_low', 'C_up")
## 

HeatingRate = HeatingRate/60.0; # --- to convert C/min to C/s; 

dummy  = matread("polyHR.mat");
extract(dummy);

# --- coefficients are calculated by average rates from each heating Rate 
# C_low = ([meanL1,meanL2])
# C_up  = ([meanU1,meanU2])
# --- coefficients are calulated by using polynomial fitting 
C_low = zeros(Float64,2);
C_up  = zeros(Float64,2);
C_low[1] = exp(P1[1]*HeatingRate+P1[2]);
C_low[2] = exp(P2[1]*HeatingRate+P2[2]);
C_up[1]  = exp(P3[1]*HeatingRate+P3[2]);
C_up[2]  = exp(P4[1]*HeatingRate+P4[2]);

# load CoeffHeatingRate.mat
# C_low[1] = exp(p1[1]*HeatingRate+p1[2])
# C_low[2] = exp(p2[1]*HeatingRate+p2[2])
# C_up = [coeffACup,coeffBCup]



#--- original data from Wakili et al. 
#  C_up  = [0.0873877575974930,336.990116244978]
#  C_low = [177262519.522505,10405.0274800779]

# ### --- trial_20 c/min()
# C_up  = [25881.19013,19483.98907]
# C_low = [608206297.2,11568.9044]

#--------------------------------------------------------------------------
#time loop
time  = 0.;
ii    = 1;
kk    = 0;

SumM   = zeros(Float64,10000,1:nx+2,1:ny+2,1:nz+2);
Ssum   = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
SphiA  = zeros(Float64,1:nx+2,1:ny+2,1:nz+2); # source term for air
rhogyp = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);

CpEffgypi = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);

HFi    = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
kgypi  = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
LHeati = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
MugAi  = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
MugVi  = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
CpVi   = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
kwi    = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
kai    = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
SphiV  = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
HeatingRate = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
dummyHR = ones(Float64,1:nx+2,1:ny+2,1:nz+2)*0.1;

Psat    = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
BP      = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
Rho_Sat = zeros(Float64,1:nx+2,1:ny+2,1:nz+2);

hfg          =  zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
SphiV_recond =  zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
SphiVnew     =  zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
Q_recond     =  zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
SphiT        =  zeros(Float64,1:nx+2,1:ny+2,1:nz+2);
timeSave = zeros(Float64, 10000);
phi_vapor_EE = zeros(Float64, 10000);
# load all_data.mat

while (time .< time_max)
ii += 1;
time          = time + dt;
Temp_old      = Temp;
phi_vapor_old = phi_vapor;
phi_air_old   = phi_air;
## --------------------------------------------------------------------
Yv     = phi_vapor.value./(phi_vapor.value+phi_air.value);
## --- Boundary condition for the temperature field -------------------

TempFaceEnd = (Temp.value[end,2:end-1,2:end-1]+Temp.value[end-1,2:end-1,2:end-1])/2.0;
heat_flux   = reshape(incident_HFlux,20,20) .- e*sigma*(TempFaceEnd.^4 .- Tgas.^4); #+hcE*(Tgas-Temp[end])
Temp_E      =  interpTfire(time)+273.0;
#--- Here take the temperature at the face value ----------------------
TempFace1 = (Temp.value[1,2:end-1,2:end-1]+Temp.value[2,2:end-1,2:end-1])/2;
dTempdx_W = (hcW*(Temp_W .- TempFace1) .+ e*sigma*(Temp_W.^4 .- TempFace1.^4))./kgypFace.xvalue[1,:,:];
dTempdx_E = heat_flux./kgypFace.xvalue[end,:,:];
#--- Boundary condition for Temperature -------------------------------
BCTemp.left.a  .= 1;
BCTemp.left.b  .= 0; 
BCTemp.left.c  = reshape(dTempdx_W,ny,nz);
BCTemp.right.a .= 1;
BCTemp.right.b .= 0;
BCTemp.right.c = reshape(dTempdx_E,ny,nz);
#----------------------------------------------------------------------
## ----------- Boundary conditions for the water vapor ----------------
phi_vapor_W   = 1000*Yv_W[1,2:end-1,2:end-1].*Pamb.*MWmix[1,2:end-1,2:end-1]/(Ru*Temp_W); # Ambient side [K]
phi_vapor_E   = 1000*Yv_E.*Pamb.*28.869/(Ru*Temp_E); # Fire side   [K]

PhiFaceW = (phi_vapor.value[1,2:end-1,2:end-1]+phi_vapor.value[2,2:end-1,2:end-1])/2;
PhiFaceE = (phi_vapor.value[end,2:end-1,2:end-1]+phi_vapor.value[end-1,2:end-1,2:end-1])/2;

#     dphi_vaporx_W = (hm*(phi_vapor_W-phi_vapor.value[1])- phi_vapor.value[1].*U_D.xvalue[1])./DeffFace.xvalue[1]
#     dphi_vaporx_E = (hm*(phi_vapor_E-phi_vapor.value[end]) - phi_vapor.value[end].*U_D.xvalue[end])./DeffFace.xvalue[end]

dphi_vaporx_W = (hm*(phi_vapor_W-PhiFaceW)- PhiFaceW.*U_D.xvalue[1,:,:])./DeffFace.xvalue[1,:,:];
dphi_vaporx_E = (hm*(phi_vapor_E .- PhiFaceE) - PhiFaceE.*U_D.xvalue[end,:,:])./DeffFace.xvalue[end-1,:,:];
#----------------------------------------------------------------------
#------- Boundary condition for water vapor ---------------------------
#     BCphiV.left.a  = 1; BCphiV.left.b  = 0; BCphiV.left.c  = dphi_vaporx_W
#     BCphiV.right.a = 1; BCphiV.right.b = 0; BCphiV.right.c = dphi_vaporx_E

BCphiV.left.a  .= 0;
BCphiV.left.b  .= 1;
BCphiV.left.c  = reshape(phi_vapor_W,ny,nz);
BCphiV.right.a .= 0;
BCphiV.right.b .= 1;
BCphiV.right.c .= phi_vapor_E;

#----------------------------------------------------------------------
## ---------- Boundary conditions for the air -------------------------
# MmixB = phi_vapor[1]*Ru*Temp[1]/(1000*Pamb) + (1- phi_vapor[1]*Ru*Temp[1]/(1000*Pamb*MWh2o))*MWair
YvW   = phi_vapor.value[1,2:end-1,2:end-1]*Ru.*Temp.value[1,2:end-1,2:end-1]./(1000*Pamb*MWmix[1,2:end-1,2:end-1]);
YvE   = phi_vapor.value[end,2:end-1,2:end-1]*Ru.*Temp.value[end,2:end-1,2:end-1]./(1000*Pamb*MWmix[end,2:end-1,2:end-1]);

phi_air_W = 1000.0*(1.0 .- YvW).*Pamb.*MWmix[1,2:end-1,2:end-1]./(Ru*TempFace1); # Ambient side [K]
phi_air_E = 1000.0*(1.0 .- YvE).*P.value[end,2:end-1,2:end-1].*MWmix[end,2:end-1,2:end-1]./(Ru*TempFaceEnd); # Fire side   [K]

PhiAFaceE = (phi_air.value[end,2:end-1,2:end-1] .+ phi_air.value[end-1,2:end-1,2:end-1])/2.0;
PhiAFaceW = (phi_air.value[1,2:end-1,2:end-1]+phi_air.value[2,2:end-1,2:end-1])/2.0;

dphi_airx_W = hm*(phi_air_W .- PhiAFaceW)./DeffFace.xvalue[1,:,:];
dphi_airx_E = 0.0;

#----------------------------------------------------------------------
#------- Boundary condition for Air -----------------------------------
BCphiA.left.a  .= 0;
BCphiA.left.b  .= 1.0;
BCphiA.left.c  = reshape(phi_air_W,ny,nz);
BCphiA.right.a .= 0.0;
BCphiA.right.b .= 1.0;
BCphiA.right.c = reshape(phi_air_E,ny,nz);
#----------------------------------------------------------------------
# MM[ii] = MmixB
## -------------- end of boundary conditions ------------------------------

#----- update thermo-physical properties & reaction rates ---------------
for i = 2:nx+1
    for j = 2:ny+1
        for k = 2:nz+1
            CpEffgypi[i,j,k]  = interpCpEff_gyp(Temp.value[i,j,k]); # Temperature values are converted to Kelvin.
            HFi[i,j,k]        = interpHeat_Flow(Temp.value[i,j,k]);# Temperature values are converted to Kelvin. & mW is converted to W.
            kgypi[i,j,k]      = interpk_gyp(Temp.value[i,j,k]);
            LHeati[i,j,k]     = interpLatentHeat(Temp.value[i,j,k]); 

            MugAi[i,j,k]      = interpMg(Temp.value[i,j,k]);
            MugVi[i,j,k]      = interpMgv(Temp.value[i,j,k]);
            CpVi[i,j,k]       = 1000.0*interpCpVapor(Temp.value[i,j,k]); # Specific heat capacity of water steam [J/kg.K].
            kwi[i,j,k]        = interpkwater(Temp.value[i,j,k]);
            kai[i,j,k]        = interpkair(Temp.value[i,j,k]);
        end
    end
end

CpEffgyp.value  .= CpEffgypi; # Temperature values are converted to Kelvin.
HF.value        .= HFi;# Temperature values are converted to Kelvin. & mW is converted to W.
kgyp.value      .= kgypi;
LHeat.value     .= LHeati;

MugA.value      .= MugAi;
MugV.value      .= MugVi;
CpV.value       .= CpVi; # Specific heat capacity of water steam [J/kg.K].
kw.value        .= kwi;
ka.value        .= kai;

Mug       = Xv.*MugVi+(1.0 .- Xv).*MugAi;
kgypFace  = arithmeticMean(kgyp);

## --------- effective rho*cp for temperature field -----------------------
rhogyp     .= reshape(rho_in .- sum(SumM,dims=1),nx+2,ny+2,nz+2);
RhoCpeff.value   .= rhogyp.*1000.0 .+ Pore.* (phi_vapor.value .* CpV.value .+ phi_air.value .* CpA.value);
kgass      = kai.*phi_air.value./(phi_vapor.value+phi_air.value) + kwi.*phi_vapor.value./(phi_vapor.value+phi_air.value);


alphaT.value .= (kgyp.value .* (1.0 .- Pore) .+ kgass .* Pore)./RhoCpeff.value; 
#alphaT = (kgyp * (1.0 - Pore) + kgass * Pore)/RhoCpeff;

alphaTFace    = arithmeticMean(alphaT);
  
## ---------- effective velocity field for temperature --------------------
Gradphi_vapor = gradientTerm(phi_vapor);   
Gradphi_air   = gradientTerm(phi_air);

CpVFace      = arithmeticMean(CpV); 
CpAFace      = arithmeticMean(CpA); 
RhoCpeffFace = arithmeticMean(RhoCpeff); 

JDiffV = arithmeticMean(RhoCpeff);
JDiffA = arithmeticMean(RhoCpeff);
JPresV = arithmeticMean(RhoCpeff);
JPresA = arithmeticMean(RhoCpeff);
#UvelTemp = arithmeticMean(RhoCpeff);
JDiffV = -DeffFace * Gradphi_vapor;
JPresV = arithmeticMean(phi_vapor)*U_D; 
JDiffA = DeffFace * Gradphi_air; 
JPresA = arithmeticMean(phi_air) * U_D;

#= 

JDiffV.xvalue .= -DeffFace.xvalue .* Gradphi_vapor.xvalue;
JDiffV.yvalue .= -DeffFace.yvalue .* Gradphi_vapor.yvalue;
JDiffV.zvalue .= -DeffFace.zvalue .* Gradphi_vapor.zvalue;

JPresV.xvalue .= arithmeticMean(phi_vapor).xvalue .* U_D.xvalue; 
JPresV.yvalue .= arithmeticMean(phi_vapor).yvalue .* U_D.yvalue; 
JPresV.zvalue .= arithmeticMean(phi_vapor).zvalue .* U_D.zvalue; 

JDiffA.xvalue .= -DeffFace.xvalue .* Gradphi_air.xvalue; 
JDiffA.yvalue .= -DeffFace.yvalue .* Gradphi_air.yvalue; 
JDiffA.zvalue .= -DeffFace.zvalue .* Gradphi_air.zvalue; 

JPresA.xvalue .= arithmeticMean(phi_air).xvalue .* U_D.xvalue;
JPresA.yvalue .= arithmeticMean(phi_air).yvalue .* U_D.yvalue;
JPresA.xvalue .= arithmeticMean(phi_air).xvalue .* U_D.xvalue;

=#

UvelTemp = (CpVFace * (JDiffV + JPresV) + CpAFace * (JDiffA + JPresA))/RhoCpeffFace; 
      
## Source term for the vapor [kg/m3.s]; Mass production consumption of the kth gas phase component per unit volume.
H2o_initial = (H20_unbound_frac + H20_bound_frac)*rho_in;
H2o_st_1 = (H20_unbound_frac + 0.75*H20_bound_frac)*rho_in;
## --- Here update heating rate to calculate Arhenious coefficients 

HeatingRate .= reshape(max(vec(dummyHR),vec((Temp.value .- Temp_old.value)/dt)),nx+2,ny+2,nz+2);      
C_low1 = exp.(P1[1]*HeatingRate .+ P1[2]);
C_low2 = exp.(P2[1]*HeatingRate .+ P2[2]);
C_up1  = exp.(P3[1]*HeatingRate .+ P3[2]);
C_up2  = exp.(P4[1]*HeatingRate .+ P4[2]);
    
##    
TTemp = Temp.value;

for i = 1:nx+2
    for j = 1:ny+2
        for k = 1:nz+2
            mask1 = TTemp[i,j,k]<=873;
            SphiV[i,j,k] = mask1*(C_low1[i,j,k]*
            exp(-C_low2[i,j,k]./TTemp[i,j,k])
            *(H2o_st_1-Ssum[i,j,k]))+(1-mask1)*
            (C_up1[1]*exp(-C_up2[2]./(TTemp[i,j,k]))
            *(H2o_initial-Ssum[i,j,k]));      
        end
    end
end


#push!(SumM,SphiV*dt);
SumM[ii,:,:,:] = SphiV*dt;
Ssum           = reshape(sum(SumM,dims=1),nx+2,ny+2,nz+2);
## Source term for the temperature [W/m3]
STP   =(Pore*(P.value-P_old.value)*1e+5/dt)./RhoCpeff.value; # source term from pressure for temperature
QT    = SphiV./RhoCpeff.value.*-1.0*450.0*1000.0;
## Water re-condensation
# determine the boiling point
      
   
for i = 2:nx+1
    for j = 2:ny+1
        for k = 2:nz+1
            if Temp.value[i,j,k]<150.0+273.0;
                Psat[i,j,k]= interpWT12(Temp.value[i,j,k]); #Pa
                BP[i,j,k]      = interpWT21(P.value[i,j,k]);
                Rho_Sat[i,j,k] = Psat[i,j,k].*MWmix[i,j,k]./(Ru.*100*Temp.value[i,j,k]); # kg/m3
            else()
                Rho_Sat[i,j,k] = phi_vapor.value[i,j,k];
            end
        end
    end
end



for i = 1:nx+2
    for j = 1:ny+2
        for k = 1:nz+2
            if Temp.value[i,j,k] .< BP[i,j,k] && phi_vapor.value[i,j,k] .> Rho_Sat[i,j,k];
                hfg[i,j,k] = 1000.0*interpWT13(BP[i,j,k]);
                #SphiV_recond[i] = abs(U_D[i])*(phi_vapor[i] - Rho_Sat[i])/dxu[i]
                SphiV_recond[i,j,k] = (phi_vapor.value[i,j,k] - Rho_Sat[i,j,k])/dt;
                SphiVnew[i,j,k]     = SphiV[i,j,k] - SphiV_recond[i,j,k];
                Q_recond[i,j,k]     = SphiV_recond[i,j,k] * hfg[i,j,k];
            else()
                Q_recond[i,j,k]     = 0;
                SphiV_recond[i,j,k] = 0;
                SphiVnew[i,j,k]     = SphiV[i,j,k];
            end
            SphiT[i,j,k] = STP[i,j,k]+QT[i,j,k]+Q_recond[i,j,k]/RhoCpeff.value[i,j,k];
       end
    end
end

##
SphiTSource.value .= SphiT;
SphiVSource.value .= SphiVnew;
SphiASource.value .= SphiA;
# SphiT[nx/2,ny/2,nz/2]
# SphiVnew[nx/2,ny/2,nz/2]
# SphiA[nx/2,ny/2,nz/2]
# ---- Here reiterate the inner values ------------------------------------
#--- Temperature 
MbcTemp, RHSbcTemp = boundaryConditionTerm(BCTemp);
MtTemp, RHStTemp   = transientTerm(Temp_old, dt,1); # create the sparse matrix for the transient term 
MdiffTemp          = diffusionTerm(alphaTFace); 
RHS_source         = constantSourceTerm(SphiTSource);

MbcphiV, RHSbcphiV = boundaryConditionTerm(BCphiV);
MtphiV, RHStphiV   = transientTerm(phi_vapor_old, dt,1); # create the sparse matrix for the transient term 
MdiffphiV          = diffusionTerm(DeffFace/Pore); 
RHS_s_phiV         = constantSourceTerm(SphiVSource);

MbcphiA, RHSbcphiA = boundaryConditionTerm(BCphiA);
MtphiA, RHStphiA   = transientTerm(phi_air_old, dt,1); # create the sparse matrix for the transient term 
MdiffphiA          = diffusionTerm(DeffFace/Pore); 
RHS_s_phiA         = constantSourceTerm(SphiASource);


## ------------- solve for temperature --------------------------------
#--- create matrixes---------------------------------------------------

MconvTemp = convectionUpwindTerm(UvelTemp);   # create the sparse matrix for convective term    
MTemp     = MtTemp-MdiffTemp+MbcTemp+MconvTemp; 
RHSTemp   = RHSbcTemp+RHStTemp+RHS_source; 
Temp      = solvePDE(m,MTemp,RHSTemp); #--- solve the PDE 
#     Temp    = solvePDE[m,full(MTemp),RHSTemp]
#     Temp    = Temp*urelax+ Temp_old*(1-urelax)
## ----------- solve for species transport -----------------------------
MconvphiV = convectionUpwindTerm(U_D/Pore);   # create the sparse matrix for convective term      
MphiV     = MtphiV-MdiffphiV+MbcphiV+MconvphiV; 
RHSphiV   = RHSbcphiV+RHStphiV+RHS_s_phiV; 
phi_vapor = solvePDE(m,MphiV,RHSphiV); #--- solve the PDE    
#     phi_vapor    = solvePDE[m,full(MphiV),RHSphiV]; #--- solve the PDE  
 phi_vapor = phi_vapor*urelax+phi_vapor_old*(1.0-urelax);
#--- solve for phi_air here   
MconvphiA  = convectionUpwindTerm(U_D/Pore);   # create the sparse matrix for convective term     
MphiA      = MtphiA-MdiffphiA+MbcphiA+MconvphiA; 
RHSphiA    = RHSbcphiA+RHStphiA+RHS_s_phiA; 
phi_air    = solvePDE(m,MphiA,RHSphiA); #--- solve the PDE 

TempFace = arithmeticMean(Temp);
phVFace  = arithmeticMean(phi_vapor);
phAFace  = arithmeticMean(phi_air);
T2Dtop = Temp.value[end,:,:]; 
## ---------- calculate pressure ------------------------------------------
Yv           = phi_vapor.value./(phi_vapor.value+phi_air.value);
Xv           = Yv.*MWair./(Yv*MWair .+ MWh2o-Yv*MWh2o);
MWmix        = Xv*MWh2o + (1.0 .- Xv)*MWair; # Molecular weight of air-water mixture [gram]
P.value      .= phi_vapor.value .* (Ru*Temp.value)./(Yv.*MWmix)/1000.0; # [bar]
P.value[1,:,:]   = 2*Pamb .- P.value[2,:,:];
P.value[end,:,:] = 2*Pamb .- P.value[end-1,:,:];
Pface        = arithmeticMean(P);
delP         = gradientTerm(P);

MugFaceCell.value .= Mug;
MugFace           = arithmeticMean(MugFaceCell);

delP.xvalue[end,:,:] = delP.xvalue[end-1,:,:];
delP.xvalue[1,:,:]   = delP.xvalue[2,:,:];
U_D                  = -1e+5*delP*(K/MugFace);  

P_old = P;

# phimon[ii] = phi_vapor[end];
    ## save data for post processing    
    if (ii-1)%(savetime/dt) .== 0
        kk = kk+1
#         TT[kk,:] = TempFace.xvalue
#         PP[kk,:] = P.value
#         phi_vapor_save[kk,:] = phVFace.xvalue
#         phi_air_save[kk,:]   = phAFace.xvalue
        timeSave[kk]         = time;
        println("time =  ",timeSave[kk], "\n")
        println("temperature at x = 8 mm:  ",TempFace.xvalue[41,Int(ny/2),Int(nz/2)],"\n") 
        phi_vapor_EE[kk] = phi_vapor_E;
    end
    
#  figure(1);visualizeCells[Temp]; 

end  ###### end of while loop 

figure;visualizeCells[Temp]
figure()
x = m.cellcenters.x


end
gypsumfvm3d()

#=
## Plot the profiles from 3D results 

T1 = Temp.value[2:end-1,ny/2+1,nz/2+1]
T2 = Temp.value[2:end-1,2,nz/2+1]
T3 = Temp.value[2:end-1,ny,nz/2+1]
plot(x,T1,"-o")
hold on
plot(x,T2,"-^")
plot(x,T3,"-.")

# plot the temperature at x = 0.008 m

##  Temperature
load Tempexp
# figure("Name','Temperature profiles','Color",[1 1 1])
# 
# # Create axes()
# axes1 = axes("Position",...
#     [0.13 0.137543032775649 0.819464012251149 0.799194056426229])
# hold(axes1,"on")
# 
# # Create multiple lines using matrix input to plot()
# plot1 = plot(timeSave,TT[:,[60,40,20,1]]-273,"LineWidth",3)
# set(plot1[1],"DisplayName','x = 0 mm, Numerical','Color",[0 0 1])
# set(plot1[2],"DisplayName','x = 4 mm, Numerical','LineStyle','--",...
#     "Color",[1 0 0])
# set(plot1[3],"DisplayName','x = 8 mm, Numerical','LineStyle','-.",...
#     "Color",[0 1 0])
# set(plot1[4],"DisplayName','x = 12 mm, Numerical','LineWidth",2,...
#     "LineStyle',':")
# 
# # timeR = timeR - 163
# plot2 = plot(timeR[timeR>0 & timeR<timeSave[end]],termo_couples[timeR>0 & timeR<timeSave[end],[12,4,8,11]])
# set(plot2[1],"DisplayName','x = 0 mm, experimental','Color",[0 0 1])
# set(plot2[2],"DisplayName','x = 4 mm, experimental','LineStyle','--",...
#     "Color",[1 0 0])
# set(plot2[3],"DisplayName','x = 8 mm, experimental','LineStyle','-.",...
#     "Color",[0 1 0])
# set(plot2[4],"DisplayName','x = 12 mm, experimental','LineWidth",2,...
#     "LineStyle',':")


## For Report 07/28/2021
# Create figure()
figure1 = figure("Name','Temperature profiles','Color",[1 1 1])

# Create axes()
axes1 = axes("Parent",figure1,...
    "Position",[0.13 0.137543032775649 0.819464012251149 0.799194056426229])
hold(axes1,"on")

# Create multiple lines using matrix input to plot()
plot1 = plot(timeSave,TT[:,[61,41,21,1]]-273,"LineWidth',3,'Color',[0 0 0],'Parent",axes1)
set(plot1[1],"DisplayName','x = 0 mm, Numerical")
set(plot1[2],"DisplayName','x = 4 mm, Numerical','LineStyle','--")
set(plot1[3],"DisplayName','x = 8 mm, Numerical','LineStyle','-.")
set(plot1[4],"DisplayName','x = 12 mm, Numerical','LineStyle',':")

# Create multiple lines using matrix input to plot()
plot2 = plot(timeR[timeR>0 & timeR<timeSave[end]],termo_couples[timeR>0 & timeR<timeSave[end],[13,4,8,11]],"LineWidth',3,'Color',[0 0 1],'Parent",axes1)
set(plot2[1],"DisplayName','x = 0 mm, experimental")
set(plot2[2],"DisplayName','x = 4 mm, experimental','LineStyle','--")
set(plot2[3],"DisplayName','x = 8 mm, experimental','LineStyle','-.")
set(plot2[4],"DisplayName','x = 12 mm, experimental','LineStyle',':")

# Create ylabel()
ylabel("Temperature [^oC]','FontSize",14)

# Create xlabel()
xlabel("Time [S]','FontSize",14)

box(axes1,"on")
# Set the remaining axes properties
set(axes1,"XGrid','on','XMinorTick','on','XTick",...
    [0:100:2000],...
    "YGrid';'on';'YMinorTick';'on';'YTick";...
    [0:50:1000])
# Create legend()
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.157065919642549 0.675890323353678 0.210441585955756 0.237552600020612],...
    "EdgeColor",[1 1 1])

# ----- Save Data to excel for Shijin -------------------------------------

fname2 = [fname1[1:4], "_num_exp_data.xlsx"]
Time_num = timeSave; 
num_data = TT[:,[60,40,20,1]]-273; 
Time_exp = timeR[timeR>0 & timeR<timeSave[end]]
exp_data = termo_couples[timeR>0 & timeR<timeSave[end],[12,4,8,11]]
TableDatanum  = table[Time_num',num_data]; 
TableDataexp  = table[Time_exp, exp_data]; 
writetable[TableDatanum,fname2,"Sheet','numdata','Range','A1"]
writetable[TableDataexp,fname2,"Sheet','expdata','Range','A1"]






# # Create plot()
# plot(T8mmexp[:,1],T8mmexp[:,2],"DisplayName','x = 8 mm, Experimental [8]",...
#     "MarkerEdgeColor",[0 1 0],...
#     "Marker';'o";...
#     "LineWidth";2;...
#     "LineStyle';'none";...
#     "Color",[0 1 0])
# 
# # Create plot()
# plot(T4mmexp[:,1],T4mmexp[:,2],"DisplayName','x = 4 mm, Experimental [8]','Marker','^",...
#     "LineWidth";2;...
#     "LineStyle';'none";...
#     "Color",[1 0 0])
# 
# # Create plot()
# plot(T0mmexp[:,1],T0mmexp[:,2],"DisplayName','x =0 mm, Experimental [8]','Marker','x",...
#     "LineWidth";2;...
#     "LineStyle';'none";...
#     "Color",[0 0 1])
# 
# # Create ylabel()
# ylabel("Temperature [^oC]','FontWeight','bold','FontName','Times New Roman")
# 
# # Create xlabel()
# xlabel("Time [min]','FontWeight','bold','FontName','Times New Roman")
# 
# # Uncomment the following line to preserve the Y-limits of the axes()
# ylim(axes1,[0 1000])
# box(axes1,"on")
# # Set the remaining axes properties
# set(axes1,"FontName','Times New Roman','FontSize',14,'FontWeight','bold")
# # Create legend()
# legend1 = legend(axes1,"show")
# set(legend1,...
#     "Position",[0.133795525966869 0.57979351931209 0.450244687266684 0.355932193287349],...
#     "FontSize";13;...
#     "FontWeight','normal")


# validation_ppt
figure1 = figure("Name','Temperature profiles','Color",[1 1 1])

# Create axes()
axes1 = axes("Parent",figure1,...
    "Position",[0.13 0.137543032775649 0.819464012251149 0.799194056426229])
hold(axes1,"on")

# Create multiple lines using matrix input to plot()
plot1 = plot(timeSave/60,TT[:,1:20:end]-273,"LineWidth',3,'Parent",axes1)
set(plot1[1],"DisplayName','x = 0 mm, Numerical [Present]','Color",[0 0 1])
set(plot1[2],"DisplayName','x = 4 mm, Numerical [Present]','LineStyle','--",...
    "Color",[1 0 0])
set(plot1[3],"DisplayName','x = 8 mm, Numerical [Present]','LineStyle','-.",...
    "Color",[0 1 0])
set(plot1[4],"DisplayName','x = 12 mm, Numerical [Present]','LineWidth",2,...
    "LineStyle',':")

# Create plot()
plot(T8mmexp[:,1],T8mmexp[:,2],"DisplayName','x = 8 mm, Experimental*",...
    "MarkerEdgeColor",[0 1 0],...
    "Marker';'o";...
    "LineWidth";2;...
    "LineStyle';'none";...
    "Color",[0 1 0])

# Create plot()
plot(T4mmexp[:,1],T4mmexp[:,2],"DisplayName','x = 4 mm, Experimental*','Marker','^",...
    "LineWidth";2;...
    "LineStyle';'none";...
    "Color",[1 0 0])

# Create plot()
plot(T0mmexp[:,1],T0mmexp[:,2],"DisplayName','x = 0 mm, Experimental*','Marker','x",...
    "LineWidth";2;...
    "LineStyle';'none";...
    "Color",[0 0 1])

# Create ylabel()
ylabel("Temperature [^oC]','FontWeight','bold','FontName','Times New Roman")

# Create xlabel()
xlabel("Time [min]','FontWeight','bold','FontName','Times New Roman")

# Uncomment the following line to preserve the Y-limits of the axes()
ylim(axes1,[0 1000])
box(axes1,"on")
# Set the remaining axes properties
set(axes1,"FontName','Times New Roman','FontSize',14,'FontWeight','bold")
# Create legend()
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.142195177893964 0.719485315552211 0.803418783677949 0.204831927012996],...
    "NumColumns";2;...
    "FontSize";15;...
    "EdgeColor",[1 1 1],...
    "FontWeight','normal")




# phi vapor


figure4 = figure("Name','Water vapor profiles','Color",[1 1 1])


axes1 = axes("Parent",figure4,...
    "Position",[0.157679626287364 0.160934883673525 0.77341578007306 0.781473240117578])
hold(axes1,"on")

plot1 = plot(timeSave/60,phi_vapor_save[:,1:20:end],"LineWidth',2,'Parent",axes1)
set(plot1[1],"DisplayName','x = 0 mm")
set(plot1[2],"DisplayName','x = 4 mm','LineStyle','--")
set(plot1[3],"DisplayName','x = 8 mm','LineStyle','-.")
set(plot1[4],"DisplayName','x = 12 mm','LineStyle',':")

ylabel("Vapor Desnity [Kg/m^3]','FontWeight','bold",...
    "FontName','Times New Roman")

xlabel("Time [min]','FontWeight','bold','FontName','Times New Roman")

box(axes1,"on")
set(axes1,"FontName','Times New Roman','FontSize',16,'FontWeight','bold")
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.567667057133842 0.664373373014608 0.324646648164352 0.253143139180639],...
    "FontSize";16;...
    "EdgeColor",[1 1 1])


# air_pressure
figure3 = figure("Name','Air profiles','Color",[1 1 1])

axes1 = axes("Parent",figure3,...
    "Position",[0.157679626287364 0.151263701362736 0.747320373712636 0.781473242544421])
hold(axes1,"on")

plot1 = plot(timeSave/60,phi_air_save[:,1:20:end],"LineWidth',2,'Parent",axes1)
set(plot1[1],"DisplayName','x = 0 mm")
set(plot1[2],"DisplayName','x = 4 mm','LineStyle','--")
set(plot1[3],"DisplayName','x = 8 mm','LineStyle','-.")
set(plot1[4],"DisplayName','x = 12 mm','LineStyle',':")

ylabel("Air Desnity [Kg/m^3]','FontWeight','bold",...
    "FontName','Times New Roman")

xlabel("Time [min]','FontWeight','bold','FontName','Times New Roman")

box(axes1,"on")
set(axes1,"FontName','Times New Roman','FontSize',16,'FontWeight','bold")

legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.564209207562917 0.652414268289552 0.324646648164352 0.253143139180639],...
    "FontSize";16;...
    "EdgeColor",[1 1 1])



#

figure2 = figure("Name','Pressure profiles','Color",[1 1 1])
axes1 = axes("Parent",figure2)
hold(axes1,"on")
plot1 = plot(timeSave/60,PP[:,1:20:end],"LineWidth',2,'Parent",axes1)
set(plot1[1],"DisplayName','x = 0 mm")
set(plot1[2],"DisplayName','x = 4 mm','LineStyle','--")
set(plot1[3],"DisplayName','x = 8 mm','LineStyle','-.")
set(plot1[4],"DisplayName','x = 12 mm','LineStyle',':")
ylabel("Pressure [Bar]")
xlabel("Time [min]")
box(axes1,"on")
set(axes1,"FontSize',12,'FontWeight','bold")
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.133890919319763 0.722842266561374 0.206484638078221 0.181919637801392])


##phi_vapor_propagation

# plot the vapor density at different time
figure()
time_m = [2.5*60/savetime 5*60/savetime 7.5*60/savetime 10*60/savetime 12.5*60/savetime 15*60/savetime]
figure()
axes1 = axes("Position",...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182])
hold(axes1,"on")

plot1 = plot(x*1000,phi_vapor_save[time_m,nx:-1:1],"LineWidth",2)
set(plot1[1],"DisplayName','t = 2.5 min")
set(plot1[2],"DisplayName','t = 5 min','LineStyle','--")
set(plot1[3],"DisplayName','t = 7.5 min','LineStyle',':")
set(plot1[4],"DisplayName','t = 10 min','LineStyle','-.")
set(plot1[5],"DisplayName','t = 12.5 min",...
    "MarkerFaceColor",[0.466666668653488 0.674509823322296 0.18823529779911],...
    "MarkerEdgeColor';'none";...
    "MarkerSize";4;...
    "Marker';'diamond";...
    "LineStyle','--")
set(plot1[6],"DisplayName','t = 15 min",...
    "MarkerFaceColor",[0.301960796117783 0.745098054409027 0.933333337306976],...
    "MarkerSize";4;...
    "Marker';'hexagram";...
    "LineStyle';'--";...
    "LineWidth",0.5)

ylabel("Vapor Desnity [Kg/m^3]")

xlabel("Distance [mm]")

# Uncomment the following line to preserve the X-limits of the axes()
xlim(axes1,[0 12])
# Uncomment the following line to preserve the Y-limits of the axes()
ylim(axes1,[0 0.8])
box(axes1,"on")
# Set the remaining axes properties
set(axes1,"FontName','Times New Roman','FontSize',14,'FontWeight','bold")
# Create legend()
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.289598115563532 0.718704518549598 0.485815593589705 0.20170453986661],...
    "NumColumns";2;...
    "FontSize";14;...
    "FontWeight','normal")

############
###figure()

time_m = [3*60/savetime 6*60/savetime 9*60/savetime 12*60/savetime 15*60/savetime 18*60/savetime]
xmm = x*1000
phivapor1 = phi_vapor_save[time_m[1],:]
phivapor2 = phi_vapor_save[time_m[2],:]
phivapor3 = phi_vapor_save[time_m[3],:]
phivapor4 = phi_vapor_save[time_m[4],:]
phivapor5 = phi_vapor_save[time_m[5],:]
phivapor6 = phi_vapor_save[time_m[6],:]
#save("phivapor5kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6")
save("phivapor10kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6")

figure()
plot(x*1000,phi_vapor_save[time_m[1],:])
hold on
plot(x*1000,phi_vapor_save[time_m[2],:])
plot(x*1000,phi_vapor_save[time_m[3],:])
plot(x*1000,phi_vapor_save[time_m[4],:])
plot(x*1000,phi_vapor_save[time_m[5],:])
plot(x*1000,phi_vapor_save[time_m[6],:])

figure()
axes1 = axes("Position",...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182])
hold(axes1,"on")

# Create multiple lines using matrix input to plot()
plot1 = plot(x*1000,phi_vapor_save[time_m,nx:-1:1],"LineWidth",2)
set(plot1[1],"DisplayName','t = 3 min")
set(plot1[2],"DisplayName','t = 6 min','LineStyle','--")
set(plot1[3],"DisplayName','t = 9 min','LineStyle',':")
set(plot1[4],"DisplayName','t = 12 min','LineStyle','-.")
set(plot1[5],"DisplayName','t = 15 min",...
    "MarkerFaceColor",[0.466666668653488 0.674509823322296 0.18823529779911],...
    "MarkerEdgeColor';'none";...
    "MarkerSize";4;...
    "Marker';'diamond";...
    "LineStyle','--")
set(plot1[6],"DisplayName','t = 18 min",...
    "MarkerFaceColor",[0.301960796117783 0.745098054409027 0.933333337306976],...
    "MarkerSize";4;...
    "Marker';'hexagram";...
    "LineStyle';'--";...
    "LineWidth",0.5)

# Create ylabel()
ylabel("Vapor Desnity [Kg/m^3]")

# Create xlabel()
xlabel("Distance [mm]")

# Uncomment the following line to preserve the X-limits of the axes()
xlim(axes1,[0 12])
# Uncomment the following line to preserve the Y-limits of the axes()
ylim(axes1,[0 0.8])
box(axes1,"on")
# Set the remaining axes properties
set(axes1,"FontName','Times New Roman','FontSize',14,'FontWeight','bold")
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.289598115563532 0.718704518549598 0.485815593589705 0.20170453986661],...
    "NumColumns";2;...
    "FontSize";14;...
    "FontWeight','normal")

#  Temperature_no_exp

figure12 = figure("Name','Temperature profiles_no_exp','Color",[1 1 1])

# Create axes()
axes1 = axes("Parent",figure12,...
    "Position",[0.169484471131572 0.151078834588102 0.779548920784069 0.786369298191981])
hold(axes1,"on")

# Create multiple lines using matrix input to plot()
plot1 = plot(timeSave/60,TT[:,1:20:end]-273,"LineWidth",3)
set(plot1[1],"DisplayName','x = 0 mm','Color",[0 0 1])
set(plot1[2],"DisplayName','x = 4 mm','LineStyle','--','Color",[1 0 0])
set(plot1[3],"DisplayName','x = 8 mm','LineStyle','-.','Color",[0 1 0])
set(plot1[4],"DisplayName','x = 12 mm','LineStyle',':")

# Create ylabel()
ylabel("Temperature [^oC]','FontWeight','bold','FontName','Times New Roman")

# Create xlabel()
xlabel("Time [min]','FontWeight','bold','FontName','Times New Roman")

# Uncomment the following line to preserve the Y-limits of the axes()
 ylim(axes1,[0 1000])
box(axes1,"on")
# Set the remaining axes properties
set(axes1,"FontName','Times New Roman','FontSize",18)
# Create legend()
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.250001893604248 0.777662520627923 0.521968355183023 0.12136929126696],...
    "NumColumns";2;...
    "FontSize";18;...
    "EdgeColor",[0.972549021244049 0.972549021244049 0.972549021244049])




#####dehydration
# plot the source term with respect to distance at different time


figure()
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[1],:]))
hold on
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[2],:]))
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[3],:]))
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[4],:]))
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[5],:]))
plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[6],:]))

set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[1],:])),"DisplayName','t = 2.5 min','LineWidth",dims=2)
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[2],:])),"DisplayName','t = 5 min','LineWidth',3,'LineStyle','--")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[3],:])),"DisplayName','t = 7.5 min','LineWidth',3,'LineStyle',':")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[4],:])),"DisplayName','t = 10 min','LineWidth',3,'LineStyle','-.")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[5],:])),"DisplayName','t = 12.5 min','MarkerSize',4,'Marker','diamond",...
    "LineWidth";2;...
    "LineStyle','--")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[6],:])),"DisplayName','t = 15 min','MarkerSize',4,'Marker','+",...
    "LineWidth";2;...
    "LineStyle','--")



##########

#time_m = [2.5*60/dt 5*60/dt 7.5*60/dt 10*60/dt 12.5*60/dt 15*60/dt]
time_m = [3*60/dt 6*60/dt 9*60/dt 12*60/dt 15*60/dt 18*60/dt]
# xmm = x*1000
# calc_perc1 = (100/H2o_initial)*sum(Sum[1:time_m[1],:])
# calc_perc2 = (100/H2o_initial)*sum(Sum[1:time_m[2],:])
# calc_perc3 = (100/H2o_initial)*sum(Sum[1:time_m[3],:])
# calc_perc4 = (100/H2o_initial)*sum(Sum[1:time_m[4],:])
# calc_perc5 = (100/H2o_initial)*sum(Sum[1:time_m[5],:])
# calc_perc6 = (100/H2o_initial)*sum(Sum[1:time_m[6],:])
# #save("calc_percent5kw.mat','xmm','calc_perc1','calc_perc2','calc_perc3','calc_perc4','calc_perc5','calc_perc6")
# save("calc_percent10kw.mat','xmm','calc_perc1','calc_perc2','calc_perc3','calc_perc4','calc_perc5','calc_perc6")

# Create figure()
figure6 = figure("Color",[1 1 1])

# Create axes()
axes1 = axes("Parent",figure6,...
    "Position",[0.114825493171472 0.165389714221258 0.62113808801214 0.77453565891307])
hold(axes1,"on")

# Create multiple lines using matrix input to plot()
#plot1 = plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m,:]),"Parent",axes1)
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[1],nx:-1:1])),"DisplayName','t = 3 min','LineWidth",dims=2)
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[2],nx:-1:1])),"DisplayName','t = 6 min','LineWidth',3,'LineStyle','--")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[3],nx:-1:1])),"DisplayName','t = 9 min','LineWidth',3,'LineStyle',':")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[4],nx:-1:1])),"DisplayName','t = 12 min','LineWidth',3,'LineStyle','-.")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[5],nx:-1:1])),"DisplayName','t = 15 min','MarkerSize',4,'Marker','diamond",...
    "LineWidth";2;...
    "LineStyle','--")
set(plot(x*1000,(100/H2o_initial)*sum(Sum[1:time_m[6],nx:-1:1])),"DisplayName','t = 18 min','MarkerSize',4,'Marker','+",...
    "LineWidth";2;...
    "LineStyle','--")

# Create ylabel()
ylabel("Percentage of dehydration [#]','FontWeight','bold",...
    "FontName','Times New Roman")

# Create xlabel()
xlabel("Distance [mm]','FontWeight','bold','FontName','Times New Roman")

# Uncomment the following line to preserve the X-limits of the axes()
xlim(axes1,[0 12])
# Uncomment the following line to preserve the Y-limits of the axes()
ylim(axes1,[0 100])
box(axes1,"on")
# Set the remaining axes properties
set(axes1,"FontName','Times New Roman','FontSize',14,'FontWeight','bold")
# Create legend()
legend1 = legend(axes1,"show")
set(legend1,...
    "Position",[0.77529302456479 0.402985074626866 0.207890739661609 0.490878938640132],...
    "EdgeColor",[1 1 1])


###############combined
time_m = [3*60/savetime 6*60/savetime 9*60/savetime 12*60/savetime 15*60/savetime 18*60/savetime]
time_m2 = [3*60/dt 6*60/dt 9*60/dt 12*60/dt 15*60/dt 18*60/dt]


figure()
plot(Tfire[:,1],Tfire[:,2])

=#


export mp, assign_materials!
T_test = 450
LH_value = GDict["LH"](T_test)
FE = create_FE()
T = FE["T"]
Sp = FE["Sp"]
LH=FE["ImF"]
#Parameters
Φ= 0.85; # Porosity
t=1400; # Time
τ=1.869; #tortuosity
D_ab=2.56e-5; # [m2/s] Binary diffusion coefficient
D_eff= Φ*D_ab/τ
Uᴰ=0.6e-5; #darcy velocity 
μ=1.862e-05 # dynamic viscosity of air
K=1.00e-15 #Permeability tensor
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
#end
