%--------------------------------------------------------------------------
% Numerical solution of the non-dimensional steady convection-diffusion
% equation for heat and mass transfer
%
%       PHI,t, + (u*PHI),x = sig*PHI,xx + C1*PHI + C2*PHI^2 + Sphi
%
%       Calcination inside the gypson board 
%
%                    Hayri Sezer 
%                  WCU 11/30/2018
%--------------------------------------------------------------------------
clc
clear all 
close all 


%--------------------------------------------------------------------------
const_1   = 0.0;
const_2   = 0;

id_scheme = 1;  % 1 use upwind, 2 use centeral differencing scheme
              % 3 exponential schem, 4 expo with scource term 
              %---- option 4 is not yet added.
%implicitness parameter
alpha    = 1; % alpha = 1; fully implicit -- alpha = 0; fully explicit
% time parameters 
dt       = 0.05;
% time_max = 60*24+20; % min
time_max = 300;
savetime = 2; % save data at each 30 seconds 
urelax   = 0.95;
% time_max = 1*great;

%Choose solver
% If ide_solve = 0 Gauss-Seidel, = 1 TDMA
id_solve   = 1;
%urelax_fac = 0.0;

maxiter = 50;
nx      = 61;

% molecular weight of water and air 
MWair = 28.97; % gram/mole
MWh2o = 18.02; % gram/molde
Ru    = 83.144598; % bar.cm3/(mol.K) : 
Pamb  = 1.01325;   % bar 

% generate the grid 
xfirst = 0.;
xtot   = 0.012; % [m] thickness of ympson board
y      = 1.25;  % [m]
z      = 1.05;  % [m] 

Volume = xtot*y*z; % [m3]

[x,dxu,dxp] = Mesh_1D(xfirst,xtot,nx);
    

load Tfire_exp; % experimentally measured temperature at fire side.
Tfire = Tfire_exp;
%load hicks_c12_t2;
%hicks = hicks_c12_t2;
% load TFire2
% Tfire = TFire2;
Temp_E = interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273;% Fire side    [K];

%% Load heat flux data from FDS 
load HTFlux_FDS

heat_flux = interp1q(HTFlux_FDS(:,1),HTFlux_FDS(:,5),0)*1000;% Fire side    [K];

%Temp_E = interp1q(hicks(:,1),hicks(:,2),0);% Fire side    [K];
Temp_W = 20+273; %interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273; % Ambient side [K];


%initial condition
for i=1:nx
    Temp(i) = Temp_W;
end
Temp_old = Temp; 
% initiliaze the drag velocity 
U_D(1:nx) = 0.0;
%% initial condition for species transport 
P(1:nx)   = Pamb; % Initial condition for Pressure; 
P_old     = P;
Pvapor    = 0.4*17.5*0.00133322;               % Water initial pressure at ambient x = 0;
Yv_E      = 0.1239;                            % Mass fraction of water vapor at fire side 
Xv_W      = Pvapor./P;                         % ambient water mole fraction at ambient side 
Xv        = Xv_W;
MWmix     = Xv_W*MWh2o + (1-Xv_W)*MWair;       % Molecular weight of air-water mixture
Yv_W      = Xv_W*MWh2o/MWmix;                  % Mass fraction of water vapor 
phi_vapor = 1000*Yv_W.*Pamb.*MWmix/(Ru*Temp_W);     % Initial density of water in the computational domain ; [kg/m3]
phi_air   = 1000*(1-Yv_W).*Pamb.*MWmix/(Ru*Temp_W); % Initial density of air in the computational domain 
hm        = 9.55e-3; % [m/s] convective mass transfer coefficient for species ransport [m/s] 
hc        = 10;    % [W/m2.K]
Pore      = 0.68;  % porosity
K         = 1e-15; % Permiability [m2] 
DAB       = 2.56e-5; % [m2/s] Binary diffusion coefficient 
Tau       = 1.869; % tortuosity
Deff(1:nx)= Pore*DAB/Tau; % effective diffusion of gases. 
e         = 0.9; % emissivity
sigma     = 5.67*10^(-8); % stefan boltzmann constant 
rho_in = 810; %initial density, kg/m3
H20_unbound_frac = 0.02; %unbound water present
H20_bound_frac = 0.20927;
%% Here upload the data for the reaction rates and thermo-physical properties 
% data is obtained from 
% 
%  K. Ghazi Wakili, E. Hugi, L. Wullschleger, T. Frank, Gypsum board in fire e
%   modeling and experimental validation. Journal of Fire Sciences 25 (2007)
% 
% Data is constructed by linear interpolation 
% 
%--------------------------------------------------------------------------
%% load data from data base
load database_gyp
%% calculate the effective latent heat of gypson board dehydration.

L1ind = find(CpEff_gyp(:,1) >100 & CpEff_gyp(:,1) < 250);
L1 = 5.6e+5; %
L1 = (trapz(CpEff_gyp(L1ind,1),CpEff_gyp(L1ind,2))-CpEff_gyp(1,2)*150);
LatentHeat(1:length(CpEff_gyp),1) = CpEff_gyp(:,1)+273;
LatentHeat(1:L1ind(1)-1,2) = 1e+25;
LatentHeat(L1ind,2) = L1;
L2ind = find(CpEff_gyp(:,1) >600 & CpEff_gyp(:,1) <780);
L2    = (trapz(CpEff_gyp(L2ind,1),CpEff_gyp(L2ind,2))-CpEff_gyp(1,2)*180);
LatentHeat(L1ind(end)+1:L2ind(1)-1,2)=1e+25;
LatentHeat(L2ind,2) = L2;
LatentHeat(L2ind(end)+1:end,2) = 1e+25;


% inds = find(CpEff_gyp(:,1) <151);
% LatentHeat(inds,2) = 1e+25;

plot(LatentHeat(:,1),LatentHeat(:,2))
%% here calculate the MLR 
for i = 2:length(rho_gyp)-1
    
    MLRnew(i,2) = -(rho_gyp(i+1,2)-rho_gyp(i-1,2))/(rho_gyp(i+1,1)-rho_gyp(i-1,1))/3;
    
end
MLRnew(1,2) = 0; 
MLRnew(2) = 0; 
MLRnew(length(rho_gyp),2) = MLRnew(length(rho_gyp)-1,2);

MLRnew(:,1) = rho_gyp(:,1)+273;


% plot(MLRnew(:,1),MLRnew(:,2))


%%

plot(LatentHeat(:,1),LatentHeat(:,2))

plot(W_Loss_rate(:,1),W_Loss_rate(:,2),'-o')
weight = Volume*interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W);
volume = 10.7e-6/interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W);
W_loss_rate(:,2) = 0.01*(W_Loss_rate(:,2)).*interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W); % this is the source term for vapor 
% W_loss_rate(:,2) = 0.01*(W_Loss_rate(:,2)).*63.6173; % this is the source term for vapor 
W_loss_rate(:,1) = W_Loss_rate(:,1)+273;


plot(CpEff_gyp(:,1),CpEff_gyp(:,2))
plot(k_gyp(:,1),k_gyp(:,2))
plot(rho_gyp(:,1),rho_gyp(:,2))
plot(Heat_Flow(:,1),Heat_Flow(:,2))
plot(HF_NoLid(:,1)+273,HF_NoLid(:,2))
for i = 1: nx
    rhogyp(i)    = interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp(i));  % temperature values are converted to Kelvin.  
    CpEffgyp(i)  = interp1q(CpEff_gyp(:,1)+273,CpEff_gyp(:,2),Temp(i)); % Temperature values are converted to Kelvin. 
    HF(i)        = interp1q(Heat_Flow(:,1)+273,Heat_Flow(:,2),Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W. 
    kgyp(i)      = interp1q(k_gyp(:,1)+273,k_gyp(:,2),Temp(i));  % temperature values are converted to Kelvin. 
    Wgyp(i)      = interp1q(W_loss(:,1)+273,W_loss(:,2),Temp(i));  % temperature values are converted to Kelvin. 
    HFLid05mm(i) = interp1q(HF_Lid05mm(:,1)+273,HF_Lid05mm(:,2)*1e-3,Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W. 
    HFLid1mm(i)  = interp1q(HF_Lid1mm(:,1)+273,HF_Lid1mm(:,2)*1e-3,Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W. 
%     HF(i)   = interp1q(HF_NoLid(:,1)+273,HF_NoLid(:,2),Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W. 
    Mug(i)       = interp1q(Mg(:,1),Mg(:,2)*1e-5,Temp(i));
    
end

CpV(1:nx)   = 1.996*1000; % Specific heat capacity of water steam [J/kg.K].
CpA(1:nx)   = 1*1000;     % Specific heat capcity of air [J/kg.K].
%% end of the thermo-physical data
%% Reaction rate coefficients; 
C_up  = [0.0873877575974930,336.990116244978];
C_low = [177262519.522505,10405.0274800779];
%--------------------------------------------------------------------------
%time loop
time  = 0.;
ii    = 1; 
kk    = 0; 
Sum(1,1:nx)  = 0; 
Ssum(1,1:nx) = 0; 
% load all_data.mat

while (time < time_max)
ii = ii+1;
time          = time + dt;
Temp_old      = Temp;
phi_vapor_old = phi_vapor;
phi_air_old   = phi_air; 
%%
Yv     = phi_vapor./(phi_vapor+phi_air);
%% ---------- Boundary condition for the temperature field ----------------
% heat_flux = 10000.0; %W/m^2
heat_flux = interp1q(HTFlux_FDS(:,1),HTFlux_FDS(:,5),time)*1000;
%Temp_Fire = 800+273; %K
Temp_E    = interp1q(Tfire(:,1)*60,Tfire(:,2),time)+273;
%Temp_E = interp1q(hicks(:,1),hicks(:,2),time);% Fire side    [K];
dTempdx_W = (hc*(Temp_W-Temp(1))+e*sigma*(Temp_W^4-Temp(1)^4))./kgyp(1);
% dTempdx_W = hc*(Temp_W-Temp(1))./kgyp(1);

% Temp_W    = 20+273; %interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273; % Ambient side [K];
%dTempdx_E = (hc*(Temp_Fire-Temp(nx))+e*sigma*(Temp_Fire^4-Temp(nx)^4))./kgyp(nx);
%dTempdx_E = (hc*(Temp_E-Temp(nx))+e*sigma*(Temp_E^4-Temp(nx)^4))./kgyp(nx);
dTempdx_E = heat_flux/kgyp(nx);
BcW       = 0; % BcE = 1, drichlet, BcE = 0, flux 
BcE       = 0; % BcE = 1, drichlet, BcE = 0, flux 
BCT       = BC(BcW,BcE,Temp_W,Temp_E,dTempdx_W,dTempdx_E);
%% ----------- Boundary conditions for the water vapor --------------------
% Xv     = Yv.*MWair./(Yv*MWair+MWh2o-Yv*MWh2o);
% MWmix  = Xv*MWh2o + (1-Xv)*MWair; % Molecular weight of air-water mixture [gram]
phi_vapor_W   = 1000*Yv_W(1).*Pamb.*MWmix(1)/(Ru*Temp_W); % Ambient side [K];
% if Temp(end)<100+273
%     phi_vapor_E   = phi_vapor_W;
% else
    phi_vapor_E   = 1000*Yv_E.*Pamb.*28.869/(Ru*Temp_E); % Fire side   [K];
% end

dphi_vaporx_W = (hm*(phi_vapor_W-phi_vapor(1))- phi_vapor(1).*U_D(1))./Deff(1);
dphi_vaporx_E = (hm*(phi_vapor_E-phi_vapor(nx)) - phi_vapor(nx).*U_D(nx))./Deff(nx);
% dphi_vaporx_E = 0; 
BcW       = 0; % BcE = 1, drichlet, BcE = 0, flux 
BcE       = 0; % BcE = 1, drichlet, BcE = 0, flux 
BCV = BC(BcW,BcE,phi_vapor_W,phi_vapor_E,dphi_vaporx_W,dphi_vaporx_E);
%% ---------- Boundary conditions for the air -----------------------------
% MmixB = phi_vapor(1)*Ru*Temp(1)/(1000*Pamb) + (1- phi_vapor(1)*Ru*Temp(1)/(1000*Pamb*MWh2o))*MWair;
YvW   = phi_vapor(1)*Ru*Temp(1)/(1000*Pamb*MWmix(1));
YvE   = phi_vapor(end)*Ru*Temp(end)/(1000*Pamb*MWmix(end));
phi_air_W = 1000*(1-YvW).*Pamb.*MWmix(1)/(Ru*Temp(1)); % Ambient side [K];
% phi_air_WW   = 1000*(1-Yv_W(1)).*Pamb.*MWmix(1)/(Ru*Temp_W); % Ambient side [K];
phi_air_E   = 1000*(1-YvE).*P(end).*MWmix(end)/(Ru*Temp(end)); % Fire side   [K];
dphi_airx_W = hm*(phi_air_W-phi_air(1))./Deff(1);
dphi_airx_E = 0;
BcW         = 1; % BcE = 1, drichlet, BcE = 0, flux 
BcE         = 1; % BcE = 1, drichlet, BcE = 0, flux 
BCA = BC(BcW,BcE,phi_air_W,phi_air_E,dphi_airx_W,dphi_airx_E);
% MM(ii) = MmixB;
%% -------------- end of boundary conditions ------------------------------

%----- update thermo-physical properties and reaction rates ---------------
  for i = 1:nx             
        rhogyp(i)    = interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp(i));  
        CpEffgyp(i)  = interp1q(CpEff_gyp(:,1)+273,CpEff_gyp(:,2),Temp(i)); 
        HF(i)        = max(0,interp1q(Heat_Flow(:,1)+273,Heat_Flow(:,2)-20,Temp(i))); 
%          HF(i)       = interp1q(HF_NoLid(:,1)+273,HF_NoLid(:,2),Temp(i));
        kgyp(i)      = interp1q(k_gyp(:,1)+273,k_gyp(:,2),Temp(i));   
      W_loss_Rate(i) = interp1q(W_loss_rate(:,1),W_loss_rate(:,2),Temp(i)); 
      LHeat(i)       = interp1q(LatentHeat(:,1),LatentHeat(:,2),Temp(i));
%       MLRS(i)        = interp1q(MLR(:,1),MLR(:,2),Temp(i));
       MLRS(i)       = interp1q(MLRnew(:,1),MLRnew(:,2),Temp(i));
      MugA(i)        = interp1q(Mg(:,1),Mg(:,2)*1e-5,Temp(i));
      MugV(i)        = interp1q(Mg(:,1),Mg(:,3)*1e-5,Temp(i));
      Mug(i)         = Xv(i).*MugV(i)+(1-Xv(i))*MugA(i);
      CpV(i)         = 1000*interp1q(CpVapor(:,1),CpVapor(:,2),Temp(i));
      kw(i)          = interp1q(kwater(:,1),kwater(:,2),Temp(i));
      ka(i)          = interp1q(kair(:,1),kair(:,2),Temp(i));
  end
  
%% --------- effective rho*cp for temperature field -----------------------
% RhoCpeff(1:nx) = 1e+6/2;
% RhoCpeff = rhogyp.* CpEffgyp;
% alphaT   = kgyp./RhoCpeff/5;
rhogyp   = rho_in - sum(Sum); 
%RhoCpeff = rhogyp.* CpEffgyp + Pore.*(phi_vapor.*CpV+phi_air.*CpA);
RhoCpeff = rhogyp.*1000 + Pore.*(phi_vapor.*CpV+phi_air.*CpA);
%RhoCpeff = (1-Pore).*rhogyp.*1000 + Pore.*(phi_vapor.*CpV+phi_air.*CpA);
%RhoCpeff = rhogyp.*1000 + (phi_vapor.*CpV+phi_air.*CpA);
% RhoCpeff = 810.*1000;
kgass = ka.*phi_air./(phi_vapor+phi_air) + kw.*phi_vapor./(phi_vapor+phi_air);
%  alphaT = (kgyp*(1-Pore)+kgass*Pore)./(810.*1000);
 alphaT = (kgyp*(1-Pore)+kgass*Pore)./RhoCpeff;
%% ---------- effective velocity field for temperature --------------------
for i = 2:nx-1
   Divphi_vapor(i) =  (phi_vapor(i+1)-phi_vapor(i-1))./(2*dxp(i));
end
% Divphi_vapor(1:nx-1) = (phi_vapor(2:end)-phi_vapor(1:end-1))./dxp(1:nx-1);
Divphi_vapor(1)      = Divphi_vapor(2);
Divphi_vapor(nx)     = Divphi_vapor(nx-1);
% Divphi_air(1:nx-1)   = (phi_air(2:end)-phi_air(1:end-1))./dxp(1:nx-1);
for i = 2:nx-1
   Divphi_air(i) =  (phi_air(i+1)-phi_air(i-1))./(2*dxp(i));
end
Divphi_air(1)      = Divphi_air(2);
Divphi_air(nx)     = Divphi_air(nx-1);

UvelTemp             = (CpV.*(phi_vapor.*U_D-Deff.*Divphi_vapor)+ ...
    CpA.*(phi_air.*U_D-Deff.*Divphi_air))./RhoCpeff;
UvelPhi = (CpV.*(phi_vapor.*U_D+Deff.*Divphi_vapor)+ ...
    CpA.*(phi_air.*U_D+Deff.*Divphi_air));


%Calculate the fluxes here with central differencing for accuracy. 

for i = 2:nx-1
    convFlux(i) = UvelPhi(i).*(Temp(i+1)-Temp(i-1))./(2*dxp(i));
    condflux(i) = kgyp(i).*(Temp(i+1)-Temp(i-1))./(2*dxp(i));
%     condFlux(i) = kgyp(i).*(condflux(i+1)-condflux(i-1))./(2*dxp(i));
end
convFlux(nx)    = convFlux(nx-1);
convFlux(1)     = convFlux(2);
condflux(nx)    = condflux(nx-1)+kgyp(nx).*(Temp(nx)-Temp(nx-1))/dxp(nx-1);
condflux(1)     = condflux(2);

for i = 2:nx-1
%     condFlux(i) = kgyp(i).*(condflux(i+1)-condflux(i-1))./(2*dxp(i));
    condFlux(i) = (condflux(i+1)-condflux(i-1))./(2*dxp(i));
end

condFlux(nx)    = condFlux(nx-1);
condFlux(1)     = condFlux(2);
%totalFlux       = abs(condFlux+convFlux);
totalFlux       = abs(condFlux-convFlux);

% convFlux(1:nx-1) = UvelPhi(1:nx-1).*(Temp(2:end)-Temp(1:end-1))./dxp(1:end-1);
% convFlux(nx)     = convFlux(nx-1);
% condflux(1:nx-1) = kgyp(1:nx-1).*(Temp(2:end)-Temp(1:end-1))./dxp(1:end-1);
% condflux(nx)     = condflux(nx-1);
% condFlux(1:nx-1) = (condflux(2:end)-condflux(1:end-1))/dxp(1:end-1);
% condFlux(nx)     = condFlux(nx-1);
%  totalFlux       = abs(condFlux+convFlux);
 
%% Source term for the vapor [kg/m3.s]; Mass production consumption of the kth gas phase component per unit volume. 
% SphiV =(Temp-Temp_old).*max(0,W_loss_Rate)/Pore/dt;
% SphiV(1:nx)  = max(0,W_loss_Rate)/3/12;
% SphiV(1:nx)  = 1.5*max(0,MLRS)/Pore; % totalFlux./LHeat;
%MW_Gypsum = 
%36.03/(36.03+136.14)= 0.20927 my phone is out of charge
%okay

H2o_initial = (H20_unbound_frac + H20_bound_frac)*rho_in;
H2o_st_1 = (H20_unbound_frac + 0.75*H20_bound_frac)*rho_in;

% Coeff = [100122842.376381,10309.7773942331];


%Coeff = [100051446.185291,10176.6765835859];

% SphiV(1:nx)  = totalFlux./LHeat;
% % SphiVSum(ii) = sum(SphiV)/100; 
 SphiA(1:nx)  = 0; % source term for air 


% for i = 1:nx
%     if Temp(i)<650
% 
%         if sum(Sum(:,i))<0.17*810 
% 
%             SphiV(i) = totalFlux(i)/LHeat(i);
%             %SphiV(i) = max(0,MLRS(i));
%             else
% 
%             SphiV(i) = 0; 
%         end
% 
%     else   
% 
%         if sum(Sum(:,i))<0.23*810 
%              SphiV(i) = totalFlux(i)/LHeat(i);
% %            SphiV(i) = max(0,MLRS(i));
%         else
%                 SphiV(i) = 0; 
% 
%         end
% 
%     end
% end

% % 
% Sum(ii,:) = SphiV*dt;

% for i = 1:nx
%     if sum(Sum(:,i))<0.23*810 
%         SphiV(i) = totalFlux(i)/LHeat(i);
%     else
%         SphiV(i) = 0;
%     end
%     
% end
%       

%SphiV(1:nx)  = max(0,MLRS);

for jj = 1:length(Temp)
    if Temp(jj) <= 873
         SphiV(jj) = C_low(1)*exp(-C_low(2)./Temp(jj)).*(H2o_st_1-Ssum(jj));
    else
        SphiV(jj)  = C_up(1)*exp(-C_up(2)./(Temp(jj)-873)).*(H2o_initial-Ssum(jj));
    end
end
%SphiV(i) = Coeff(1)*exp(-Coeff(2)./Temp)*(H2o_st_1-sum(Sum(:i))*dt);
Sum(ii,:)    = SphiV*dt;
Ssum         = sum(Sum);
%% Source term for the temperature [W/m3]
STP   =(Pore*(P-P_old)*1e+5/dt)./RhoCpeff; % source term from pressure for temperature 
P_old = P; 
%QT = 0; % yeah
QT = SphiV./RhoCpeff.*-1*450*1000;
%QT = SphiV./RhoCpeff.*-1*L1;
%RhoCpeff = (1-Pore).*rhogyp.* CpEffgyp + Pore.*(phi_vapor.*CpV+phi_air.*CpA);
% for jj = 1:length(Temp)
%     if Temp(jj)<600+273
%        QT(jj)    =SphiV(jj)./RhoCpeff(jj).*-1*L1;
%     else
%        QT(jj)    =SphiV(jj)./RhoCpeff(jj).*-1*L2;
%     end
% end
         % -HF/volume/1000./RhoCpeff; % heat production/consumption per unit volume;
%QT    = -1 * SphiV(1:nx).*LHeat;
SphiT = STP+QT;

%% ------------- solve for temperature ------------------------------------
Temp = sclarTrans_eq(nx,dxp,dxu,alpha,dt,UvelTemp,const_1,const_2,alphaT,...
    id_scheme,Temp_old,BCT,SphiT);
Temp = Temp*(1-urelax)+ Temp*urelax;

%% ----------- solve for species transpor ---------------------------------
phi_vapor = sclarTrans_eq(nx,dxp,dxu,alpha,dt,U_D/Pore,const_1,const_2,Deff/Pore,...
    id_scheme,phi_vapor_old,BCV,SphiV);
phi_vapor = phi_vapor*(1-urelax)+phi_vapor_old*urelax; 

phi_air   = sclarTrans_eq(nx,dxp,dxu,alpha,dt,U_D/Pore,const_1,const_2,Deff/Pore,...
    id_scheme,phi_air_old,BCA,SphiA);
phi_air = phi_air*(1-urelax)+phi_air_old*urelax; 

%% ---------- calculate pressure ------------------------------------------
Yv     = phi_vapor./(phi_vapor+phi_air);
Xv     = Yv.*MWair./(Yv*MWair+MWh2o-Yv*MWh2o);
MWmix  = Xv*MWh2o + (1-Xv)*MWair; % Molecular weight of air-water mixture [gram]
P      = phi_vapor.*(Ru.*Temp)./(Yv.*MWmix)/1000; % [bar]
% P      = (phi_vapor+phi_air).*(Ru.*Temp)./(MWmix)/1000; % [bar]
P(1)   = Pamb; 
P(end) = Pamb;

P = P*(1-urelax)+ P_old*urelax;

for i = 2:nx-1
    delP(i) = (P(i+1)-P(i-1))/(2*dxp(i));
end
delP(1)  = delP(2);
delP(nx) = delP(nx-1);

U_D = -1e+5*delP.*(K./Mug);
% delP         = 1e+5*(P(2:end)-P(1:end-1))./dxp(1:nx-1);   % convert bar to pascal;
% U_D  = -(K./Mug(1:nx-1)).*delP;
% U_D(nx)      = U_D(nx-1);
%  U_D (1:nx)= 0; 
% dphi_vaporx_E
% phi_vapor_E
% bond = phi_vapor(nx)
phimon(ii) = phi_vapor(end);
%% save data for post processing 
    
    if mod(ii,savetime/dt) == 0    
        kk = kk+1;
        TT(kk,:) = Temp;
        PP(kk,:) = P;
        phi_vapor_save(kk,:) = phi_vapor; 
        phi_air_save(kk,:)   = phi_air; 
        timeSave(kk)         = time; 
        fprintf('time =  %12.8f\n',timeSave(kk));
        fprintf('temperature at x = 8 mm:  %12.8f\n',Temp(41));
        phi_vapor_EE(kk) = phi_vapor_E;
    end
end

figure
plot(x,Temp,'-o')


% plot the temperature at x = 0.008 m; 

%%  Temperature
load Tempexp
figure('Name','Temperature profiles','Color',[1 1 1]);

% Create axes
axes1 = axes('Position',...
    [0.13 0.137543032775649 0.819464012251149 0.799194056426229]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(timeSave/60,TT(:,1:20:end)-273,'LineWidth',3);
set(plot1(1),'DisplayName','x = 0 mm, Numerical (Present)','Color',[0 0 1]);
set(plot1(2),'DisplayName','x = 4 mm, Numerical (Present)','LineStyle','--',...
    'Color',[1 0 0]);
set(plot1(3),'DisplayName','x = 8 mm, Numerical (Present)','LineStyle','-.',...
    'Color',[0 1 0]);
set(plot1(4),'DisplayName','x = 12 mm, Numerical (Present)','LineWidth',2,...
    'LineStyle',':');


% Create plot
plot(T8mmexp(:,1),T8mmexp(:,2),'DisplayName','x = 8 mm, Experimental [8]',...
    'MarkerEdgeColor',[0 1 0],...
    'Marker','o',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create plot
plot(T4mmexp(:,1),T4mmexp(:,2),'DisplayName','x = 4 mm, Experimental [8]','Marker','^',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);

% Create plot
plot(T0mmexp(:,1),T0mmexp(:,2),'DisplayName','x =0 mm, Experimental [8]','Marker','x',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 0 1]);

% Create ylabel
ylabel('Temperature (^oC)','FontWeight','bold','FontName','Times New Roman');

% Create xlabel
xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1000]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.133795525966869 0.57979351931209 0.450244687266684 0.355932193287349],...
    'FontSize',13,...
    'FontWeight','normal');


%%%%%%%%validation_ppt
figure1 = figure('Name','Temperature profiles','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.137543032775649 0.819464012251149 0.799194056426229]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(timeSave/60,TT(:,1:20:end)-273,'LineWidth',3,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm, Numerical (Present)','Color',[0 0 1]);
set(plot1(2),'DisplayName','x = 4 mm, Numerical (Present)','LineStyle','--',...
    'Color',[1 0 0]);
set(plot1(3),'DisplayName','x = 8 mm, Numerical (Present)','LineStyle','-.',...
    'Color',[0 1 0]);
set(plot1(4),'DisplayName','x = 12 mm, Numerical (Present)','LineWidth',2,...
    'LineStyle',':');

% Create plot
plot(T8mmexp(:,1),T8mmexp(:,2),'DisplayName','x = 8 mm, Experimental*',...
    'MarkerEdgeColor',[0 1 0],...
    'Marker','o',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create plot
plot(T4mmexp(:,1),T4mmexp(:,2),'DisplayName','x = 4 mm, Experimental*','Marker','^',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);

% Create plot
plot(T0mmexp(:,1),T0mmexp(:,2),'DisplayName','x = 0 mm, Experimental*','Marker','x',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 0 1]);

% Create ylabel
ylabel('Temperature (^oC)','FontWeight','bold','FontName','Times New Roman');

% Create xlabel
xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 1000]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.142195177893964 0.719485315552211 0.803418783677949 0.204831927012996],...
    'NumColumns',2,...
    'FontSize',15,...
    'EdgeColor',[1 1 1],...
    'FontWeight','normal');




%% phi vapor
figure4 = figure('Name','Water vapor profiles','Color',[1 1 1]);
axes1 = axes('Parent',figure4);
hold(axes1,'on');
plot1 = plot(timeSave/60,phi_vapor_save(:,1:20:end),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm');
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--');
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.');
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');
ylabel('Vapor Desnity (Kg/m^3)');
xlabel('Time (min)');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',12,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.133890919319763 0.722842266561374 0.206484638078221 0.181919637801392]);

%%
figure3 = figure('Name','Air profiles','Color',[1 1 1]);
axes1 = axes('Parent',figure3);
hold(axes1,'on');
plot1 = plot(timeSave/60,phi_air_save(:,1:20:end),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm');
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--');
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.');
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');
ylabel('Air Desnity (Kg/m^3)');
xlabel('Time (min)');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',12,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.133890919319763 0.722842266561374 0.206484638078221 0.181919637801392]);


%%

figure2 = figure('Name','Pressure profiles','Color',[1 1 1]);
axes1 = axes('Parent',figure2);
hold(axes1,'on');
plot1 = plot(timeSave/60,PP(:,1:20:end),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm');
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--');
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.');
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');
ylabel('Pressure (Bar)');
xlabel('Time (min)');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',12,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.133890919319763 0.722842266561374 0.206484638078221 0.181919637801392]);


%%%phi_vapor_propagation

% plot the vapor density at different time 
figure 
time_m = [1*60/savetime 2*60/savetime 3*60/savetime 4*60/savetime 5*60/savetime];
figure
axes1 = axes('Position',...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(x*1000,phi_vapor_save(time_m,:),'LineWidth',2);
set(plot1(1),'DisplayName','t = 2.5 min');
set(plot1(2),'DisplayName','t = 5 min','LineStyle','--');
set(plot1(3),'DisplayName','t = 7.5 min','LineStyle',':');
set(plot1(4),'DisplayName','t = 10 min','LineStyle','-.');
set(plot1(5),'DisplayName','t = 12.5 min',...
    'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerEdgeColor','none',...
    'MarkerSize',4,...
    'Marker','diamond',...
    'LineStyle','--');
% set(plot1(6),'DisplayName','t = 15 min',...
%     'MarkerFaceColor',[0.301960796117783 0.745098054409027 0.933333337306976],...
%     'MarkerSize',4,...
%     'Marker','hexagram',...
%     'LineStyle','--',...
%     'LineWidth',0.5);

% Create ylabel
ylabel('Vapor Desnity (Kg/m^3)');

% Create xlabel
xlabel('Distance (mm)');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 0.8]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.289598115563532 0.718704518549598 0.485815593589705 0.20170453986661],...
    'NumColumns',2,...
    'FontSize',14,...
    'FontWeight','normal');

%%%%%%%%%%%%%
%%%%figure

time_m = [1*60/savetime 2*60/savetime 3*60/savetime 4*60/savetime 5*60/savetime];
% xmm = x*1000;
% phivapor1 = phi_vapor_save(time_m(1),:);
% phivapor2 = phi_vapor_save(time_m(2),:);
% phivapor3 = phi_vapor_save(time_m(3),:);
% phivapor4 = phi_vapor_save(time_m(4),:);
% phivapor5 = phi_vapor_save(time_m(5),:);
% phivapor6 = phi_vapor_save(time_m(6),:);
% %save('phivapor5kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6');
% save('phivapor10kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6');

% figure
% plot(x*1000,phi_vapor_save(time_m(1),:))
% hold on 
% plot(x*1000,phi_vapor_save(time_m(2),:))
% plot(x*1000,phi_vapor_save(time_m(3),:))
% plot(x*1000,phi_vapor_save(time_m(4),:))
% plot(x*1000,phi_vapor_save(time_m(5),:))
% plot(x*1000,phi_vapor_save(time_m(6),:))
% Create axes
figure
axes1 = axes('Position',...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(x*1000,phi_vapor_save(time_m,:),'LineWidth',2);
set(plot1(1),'DisplayName','t = 3 min');
set(plot1(2),'DisplayName','t = 6 min','LineStyle','--');
set(plot1(3),'DisplayName','t = 9 min','LineStyle',':');
set(plot1(4),'DisplayName','t = 12 min','LineStyle','-.');
set(plot1(5),'DisplayName','t = 15 min',...
    'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerEdgeColor','none',...
    'MarkerSize',4,...
    'Marker','diamond',...
    'LineStyle','--');
% set(plot1(6),'DisplayName','t = 18 min',...
%     'MarkerFaceColor',[0.301960796117783 0.745098054409027 0.933333337306976],...
%     'MarkerSize',4,...
%     'Marker','hexagram',...
%     'LineStyle','--',...
%     'LineWidth',0.5);

% Create ylabel
ylabel('Vapor Desnity (Kg/m^3)');

% Create xlabel
xlabel('Distance (mm)');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 0.8]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.289598115563532 0.718704518549598 0.485815593589705 0.20170453986661],...
    'NumColumns',2,...
    'FontSize',14,...
    'FontWeight','normal');

%%  Temperature_no_exp
figure('Name','Temperature profiles_no_exp','Color',[1 1 1]);

% Create axes
axes1 = axes('Position',...
    [0.189690721649485 0.172209699442316 0.779381443298969 0.799194056426229]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(timeSave/60,TT(:,1:20:end)-273,'LineWidth',3);
set(plot1(1),'DisplayName','x = 0 mm','Color',[0 0 1]);
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--','Color',[1 0 0]);
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.','Color',[0 1 0]);
set(plot1(4),'DisplayName','x = 12 mm','LineWidth',2,'LineStyle',':');

% Create ylabel
ylabel('Temperature (^oC)','FontWeight','bold','FontName','Times New Roman');

% Create xlabel
xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 25]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 1000]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.0866116583947029 0.598437587714848 0.51987015436978 0.398644056481831],...
    'FontSize',18,...
    'FontWeight','normal');
figure

%%%%%%dehydration
% plot the source term with respect to distance at different time 


% figure
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),:)))
% hold on 
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),:)))
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),:)))
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),:)))
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),:)))
% plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),:)))
% 
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),:))),'DisplayName','t = 2.5 min','LineWidth',2);
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),:))),'DisplayName','t = 5 min','LineWidth',3,'LineStyle','--');
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),:))),'DisplayName','t = 7.5 min','LineWidth',3,'LineStyle',':');
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),:))),'DisplayName','t = 10 min','LineWidth',3,'LineStyle','-.');
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),:))),'DisplayName','t = 12.5 min','MarkerSize',4,'Marker','diamond',...
%     'LineWidth',2,...
%     'LineStyle','--');
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),:))),'DisplayName','t = 15 min','MarkerSize',4,'Marker','+',...
%     'LineWidth',2,...
%     'LineStyle','--');
% 


%%%%%%%%%%

%time_m = [2.5*60/dt 5*60/dt 7.5*60/dt 10*60/dt 12.5*60/dt 15*60/dt];
time_m = [1*60/dt 2*60/dt 3*60/dt 4*60/dt 5*60/dt];
% xmm = x*1000;
% calc_perc1 = (100/H2o_initial)*sum(Sum(1:time_m(1),:));
% calc_perc2 = (100/H2o_initial)*sum(Sum(1:time_m(2),:));
% calc_perc3 = (100/H2o_initial)*sum(Sum(1:time_m(3),:));
% calc_perc4 = (100/H2o_initial)*sum(Sum(1:time_m(4),:));
% calc_perc5 = (100/H2o_initial)*sum(Sum(1:time_m(5),:));
% calc_perc6 = (100/H2o_initial)*sum(Sum(1:time_m(6),:));
% %save('calc_percent5kw.mat','xmm','calc_perc1','calc_perc2','calc_perc3','calc_perc4','calc_perc5','calc_perc6');
% save('calc_percent10kw.mat','xmm','calc_perc1','calc_perc2','calc_perc3','calc_perc4','calc_perc5','calc_perc6');

% Create figure
figure6 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure6,...
    'Position',[0.114825493171472 0.165389714221258 0.62113808801214 0.77453565891307]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
%plot1 = plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m,:)),'Parent',axes1);
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),:))),'DisplayName','t = 1 min','LineWidth',2);
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),:))),'DisplayName','t = 2 min','LineWidth',3,'LineStyle','--');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),:))),'DisplayName','t = 3 min','LineWidth',3,'LineStyle',':');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),:))),'DisplayName','t = 4 min','LineWidth',3,'LineStyle','-.');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),:))),'DisplayName','t = 5 min','MarkerSize',4,'Marker','diamond',...
    'LineWidth',2,...
    'LineStyle','--');
% set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),:))),'DisplayName','t = 18 min','MarkerSize',4,'Marker','+',...
%     'LineWidth',2,...
%     'LineStyle','--');

% Create ylabel
ylabel('Percentage of dehydration (%)','FontWeight','bold',...
    'FontName','Times New Roman');

% Create xlabel
xlabel('Distance (mm)','FontWeight','bold','FontName','Times New Roman');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 100]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.77529302456479 0.402985074626866 0.207890739661609 0.490878938640132],...
    'EdgeColor',[1 1 1]);


%%%%%%%%%%%%%%%combined
time_m = [1*60/savetime 2*60/savetime 3*60/savetime 4*60/savetime 5*60/savetime];
time_m2 = [1*60/dt 2*60/dt 3*60/dt 4*60/dt 5*60/dt];

figure1 = figure;
% plotyy(x*1000,phi_vapor_save(time_m(1),:),x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(1),:)))
% hold on 
% plotyy(x*1000,phi_vapor_save(time_m(2),:),x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(2),:)))
% plotyy(x*1000,phi_vapor_save(time_m(4),:),x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(4),:)))
% plotyy(x*1000,phi_vapor_save(time_m(6),:),x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(6),:)))

axes1 = axes('Parent',figure1,...
    'Position',[0.0938833570412518 0.109207708779443 0.811805265713243 0.819955315096151]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(x*1000,phi_vapor_save(time_m,:),'Parent',axes1,'LineWidth',2);
set(plot1(1),'DisplayName','t = 3 min','Color',[1 0 0]);
set(plot1(2),'DisplayName','t = 6 min','LineStyle','--','Color',[0 0 1]);
set(plot1(3),'DisplayName','t = 12 min','LineStyle',':','Color',[0 1 0]);
set(plot1(4),'DisplayName','t = 18 min','LineStyle','-.',...
    'Color',[1 0.600000023841858 0.7843137383461]);

% Create ylabel
ylabel('Vapor density (kg/m^3)');

% Create xlabel
xlabel('Distance (mm)');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 0.6]);
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',16,'XGrid','on','YColor',...
    [0.494 0.184 0.556],'YTick',[0 0.2 0.4 0.6],'YTickLabel',...
    {'0','0.2','0.4','0.6'});
% Create axes
axes2 = axes('Parent',figure1,...
    'ColorOrder',[0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741]);
hold(axes2,'on');

% Create plot
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(1),:)),'Parent',axes2,'MarkerSize',4,'Marker','o','LineWidth',2,...
    'Color',[1 0 0]);

% Create ylabel
ylabel('Percentage of dehydration (%)');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes2,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes2,[0 100]);
% Set the remaining axes properties
set(axes2,'Color','none','FontName','Times New Roman','FontSize',16,...
    'HitTest','off','YAxisLocation','right','YColor',[0.85 0.325 0.098],'YTick',...
    [0 20 40 60 80 100]);
% Create axes
axes3 = axes('Parent',figure1,...
    'ColorOrder',[0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741]);
axis off
hold(axes3,'on');

% Create plot
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(2),:)),'Parent',axes3,'MarkerSize',4,'Marker','o','LineWidth',2,...
    'LineStyle','--',...
    'Color',[0 0 1]);

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes3,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes3,[0 100]);
% Set the remaining axes properties
set(axes3,'Color','none','HitTest','off','YAxisLocation','right','YColor',...
    [0.85 0.325 0.098],'YTick',zeros(1,0));
% Create axes
axes4 = axes('Parent',figure1,...
    'ColorOrder',[0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741]);
axis off
hold(axes4,'on');

% Create plot
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m2(4),:)),'Parent',axes4,'MarkerSize',4,'Marker','o','LineWidth',2,...
    'LineStyle',':',...
    'Color',[0 1 0]);

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes4,[0 12]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes4,[0 100]);
% Set the remaining axes properties
set(axes4,'Color','none','HitTest','off','YAxisLocation','right','YColor',...
    [0.85 0.325 0.098],'YTick',zeros(1,0));
% Create axes
axes5 = axes('Parent',figure1,...
    'ColorOrder',[0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741],...
    'Position',[0.0938833570412518 0.109207708779443 0.811805265713243 0.819955315096151]);
axis off
hold(axes5,'on');



 figure;
 plot(Tfire(:,1),Tfire(:,2));
%%%%