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
fname2 = 'wakiliValidation.xlsx';
%--------------------------------------------------------------------------
const_1   = 0.0;
const_2   = 0;

%% Discretization scheme for the convection term 
id_scheme = 1;  % 1 use upwind, 2 use centeral differencing scheme
% 3 exponential schem, 4 expo with scource term
%---- option 4 is not yet added.
%% implicitness parameter
alpha    = 1; %0.55; % alpha = 1; fully implicit -- alpha = 0; fully explicit
%% time parameters
dt       =.1;
time_max = 60*18; % 24 min
% time_max = 400; 
savetime = 2; % save data at each 30 seconds
urelax   = .1;
%% Choose solver - matrix solver 
% If ide_solve = 0 Gauss-Seidel, = 1 TDMA
id_solve   = 1;
%urelax_fac = 0.0;
%% molecular weight of water and air, Pressure and Gas constant
MWair = 28.97; % gram/mole
MWh2o = 18.02; % gram/molde
Ru    = 83.144598; % bar.cm3/(mol.K) :
Pamb  = 1.01325;   % bar
%% generate the grid
maxiter = 50;
nx      = 61;
xfirst = 0.;
xtot   = 0.012; % [m] thickness of ympson board
y      = 1.25;  % [m]
z      = 1.05;  % [m]
Volume = xtot*y*z; % [m3]
[x,dxu,dxp] = Mesh_1D(xfirst,xtot,nx);
%% Upload the experimental temperature from the fire side 
load Tfire_exp; % experimentally measured temperature at fire side.
Tfire = Tfire_exp;
%load hicks_c12_t2;
%hicks = hicks_c12_t2;
% load TFire2
% Tfire = TFire2;
%% use the interpolation to get the temperature corresponding to a specific time 
Temp_E = interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273;% Fire side    [K];
%Temp_E = interp1q(hicks(:,1),hicks(:,2),0);% Fire side    [K];
Temp_W = 29+273; %interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273; % Ambient side [K];

%% initial condition
Temp(1:nx) = Temp_W;
Temp_old = Temp;
%% initiliaze the drag velocity
U_D(1:nx) = 0.0;
%% EKU experimental data for model comparison
% % % % fname1  = '90kw_flametrial_7minREDO.xlsx'; 
% % % % [~, ~, raw, dates] = xlsread(fname1,'Untitled','A2:AH19561','',@convertSpreadsheetExcelDates);
% % % % raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
% % % % raw = raw(:,[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]);
% % % % dates = dates(:,1);
% % % % R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
% % % % raw(R) = {NaN}; % Replace non-numeric cells
% % % % R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
% % % % dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN
% % % % data = reshape([raw{:}],size(raw));
% % % % timeR = data(:,1);
% % % % termo_couples = data(:,2:2:26);
% % % % %% ---- Here we should cut the initial data from the experiments. Because in the experiments the heater was not on for initial time of 2 minutes approximately. 
% % % % cutinT = 145; 
% % % % timeRR = timeR(timeR>cutinT); 
% % % % termo_couples = termo_couples(timeR>cutinT,:); 
% % % % timeRR = timeRR-timeRR(1); 
% % % % 
% % % % timeR = timeRR; 
% % % % %
% % % % %
% % % % % Create figure
% % % % figure1 = figure('Color',[1 1 1]);
% % % % axes1 = axes('Parent',figure1);
% % % % hold(axes1,'on');
% % % % plot(timeR,termo_couples);
% % % % ylabel('Temperature (^oC)');
% % % % xlabel('time (s)');
% % % % box(axes1,'on');
% % % % hold(axes1,'off');
% % % % set(axes1,'FontSize',12,'FontWeight','bold','XGrid','on','YGrid','on');
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
hcE        = 25;    % [W/m2.K]
hcW        = 9;    % [W/m2.K]
Pore      = 0.85;  % porosity
K         = 1e-15; % Permiability [m2]
DAB       = 2.56e-5; % [m2/s] Binary diffusion coefficient
Tau       = 1.869; % tortuosity
Deff(1:nx)= Pore*DAB/Tau; % effective diffusion of gases.
e         = 0.85; % emissivity
sigma     = 5.67*10^(-8); % stefan boltzmann constant
Tgas      = 303; % Tgas is the mean temperature at the gas phase on the fire side.
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
% load database_gyp  % --- old data base 
load ChemistryDatabase.mat    % --- New data base 
% load LatentHeat;
% load rho_gyp; % load the density data, rho_gyp [kg/m3], T [C]
% load CpEff_gyp; % load the specific heat data; CpEff_gyp [J/Kg.K], T [C]
% load Heat_Flow; % load the heat flow endothermic reaction of gypson board; HF [mW], T [C]
% load k_gyp; % load the thermal conductivity data, k_gyp [W/m.K], T [C]
% load W_loss; % load weight loss: W_loss = [%loss]
% load W_Loss_rate; % weight loss rate percentage.
% load HF_Lid05mm; % load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry (DSC) of pure CaSO4.2H2O
% load HF_Lid1mm; % load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry (DSC) of pure CaSO4.2H2O
% load HF_NoLid; % load the heat flow endothermic reaction;  Curve of the differential scanning calorimetry (DSC) of pure CaSO4.2H2O
% load Mg; % load data for the dynamic viscosity of air;
% load MLR;
% load kgas;
% % load MLRnew
% load CpVapor;
%% ------------------------------------------------------------------------
HeatingRate = 50; % 20 degree/min is converted to degree/second 
%% ------------------------------------------------------------------------
str1 = mat2str(HeatingRate);
fname = ['HeatingRate', str1,'.xlsx'];
str2 = 'W_loss = W_loss_White_';
str3 = [str2,str1,';'];
eval(str3); 

% W_loss = W_loss_White_20(:,[1,2]); 
% W_Loss_rate = W_loss_White_20(:,[1,3]); 
W_loss(:,2) = 100*W_loss(:,2)/W_loss(1,2);
rho_gyp = W_loss; 
rho_gyp(:,2) = rho_in*W_loss(:,2)/100;
Tempdata = rho_gyp(:,1)+273;
% 
% p = polyfit(Tempdata,rho_gyp(:,2),10); 
% 
% rho_gypFit(:,2) = smooth(Tempdata,rho_gyp(:,2),0.1,'rloess');
% 
% % rho_gypFit(:,2) = polyval(p,Tempdata); 
rhodiff    = -diff(rho_gyp(:,2))./diff(Tempdata)*HeatingRate/60;
% rhodiffFit = -diff(rho_gypFit(:,2))./diff(Tempdata)*HeatingRate/60;
% 
% figure
% plot(Tempdata(1:end-1), rhodiff,'o')
% hold on 
% plot(Tempdata(1:end-1), rhodiffFit,'-')
% figure 
% plot(Tempdata, rho_gyp(:,2),'o')
% hold on 
% plot(Tempdata, rho_gypFit(:,2),'-')

% rho_gyp(:,2) = rho_gypFit(:,2); % uncomment this line when smoothing is needed
%% calculate the effective latent heat of gypson board dehydration.
L1ind = find(CpEff_gyp(:,1) >100 & CpEff_gyp(:,1) < 250);
L1 = 5.6e+5; %

%--- Numerical integration 
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

% plot(LatentHeat(:,1),LatentHeat(:,2))
%% here calculate the MLR
% for i = 2:length(rho_gyp)-1
%     
%     MLRnew(i,2) = -(rho_gyp(i+1,2)-rho_gyp(i-1,2))/(rho_gyp(i+1,1)-rho_gyp(i-1,1))/3;
%     
% end
% MLRnew(1,2) = 0;
% MLRnew(2) = 0;
% MLRnew(length(rho_gyp),2) = MLRnew(length(rho_gyp)-1,2);
% 
% MLRnew(:,1) = rho_gyp(:,1)+273;


% plot(MLRnew(:,1),MLRnew(:,2))


%%

plot(LatentHeat(:,1),LatentHeat(:,2))

% plot(W_Loss_rate(:,1),W_Loss_rate(:,2),'-o')
weight = Volume*interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W);
volume = 10.7e-6/interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W);
% W_loss_rate(:,2) = 0.01*(W_Loss_rate(:,2)).*interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp_W); % this is the source term for vapor


% plot(W_Loss_rate(:,1), W_loss_rate(:,2))
% W_loss_rate(:,2) = 0.01*(W_Loss_rate(:,2)).*63.6173; % this is the source term for vapor
% W_loss_rate(:,1) = W_Loss_rate(:,1)+273;


% plot(CpEff_gyp(:,1),CpEff_gyp(:,2))
% plot(k_gyp(:,1),k_gyp(:,2))
% plot(rho_gyp(:,1),rho_gyp(:,2))
% plot(Heat_Flow(:,1),Heat_Flow(:,2))
% plot(HF_NoLid(:,1)+273,HF_NoLid(:,2))

for i = 1: nx
%     rhogyp(i)    = interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp(i));  % temperature values are converted to Kelvin.
    CpEffgyp(i)  = interp1q(CpEff_gyp(:,1)+273,CpEff_gyp(:,2),Temp(i)); % Temperature values are converted to Kelvin.
    HF(i)        = interp1q(Heat_Flow(:,1)+273,Heat_Flow(:,2),Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W.
    kgyp(i)      = interp1q(k_gyp(:,1)+273,k_gyp(:,2),Temp(i));  % temperature values are converted to Kelvin.
%     Wgyp(i)      = interp1q(W_loss(:,1)+273,W_loss(:,2),Temp(i));  % temperature values are converted to Kelvin.
    HFLid05mm(i) = interp1q(HF_Lid05mm(:,1)+273,HF_Lid05mm(:,2)*1e-3,Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W.
    HFLid1mm(i)  = interp1q(HF_Lid1mm(:,1)+273,HF_Lid1mm(:,2)*1e-3,Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W.
    %     HF(i)   = interp1q(HF_NoLid(:,1)+273,HF_NoLid(:,2),Temp(i)); % Temperature values are converted to Kelvin. and mW is converted to W.
    Mug(i)       = interp1q(Mg(:,1),Mg(:,2)*1e-5,Temp(i));
    
end

CpV(1:nx)   = 1.996*1000; % Specific heat capacity of water steam [J/kg.K].
CpA(1:nx)   = 1*1000;     % Specific heat capcity of air [J/kg.K].
%% end of the thermo-physical data
%% Reaction rate coefficients;
%% 
% range1 = [310,415]'; %% [K] temperature interval where we perform the fitting
% Tind = range1(1)<Tempdata & Tempdata<range1(2);
% [C_low,TableData1] = Reaction_Rate(rho_gyp,rhodiff,Tind,Tempdata,1);
% C_low = C_low';
% 
% range2 = [700,1000]';
% Tind = range2(1)<Tempdata & Tempdata<range2(2);
% [C_up,TableData2] = Reaction_Rate(rho_gyp,rhodiff,Tind,Tempdata,2);
% C_up = C_up';
% str4 = 'coeff_'; 
% str4 = [str4,str1]; 
% table = table(range1,range2,C_low,C_up); 
% 
% % Table = join(TableData1 , TableData2);
% writetable(TableData1,fname,'Sheet','data','Range','A1')
% writetable(TableData2,fname,'Sheet','data','Range','D1')
% writetable(table,fname,'Sheet',1,'Range','G1')
% save(str4, 'C_low', 'C_up');
%% 

HeatingRate = HeatingRate/60; % --- to convert C/min to C/s; 
load polyHR
% --- coefficients are calculated by average rates from each heating Rate 
% C_low = ([meanL1,meanL2]);
% C_up  = ([meanU1,meanU2]);
% --- coefficients are calulated by using polynomial fitting 
% C_low(1) = exp(P1(1)*HeatingRate+P1(2));
% C_low(2) = exp(P2(1)*HeatingRate+P2(2));
% C_up(1)  = exp(P3(1)*HeatingRate+P3(2));
% C_up(2)  = exp(P4(1)*HeatingRate+P4(2));

% load CoeffHeatingRate.mat
% C_low(1) = exp(p1(1)*HeatingRate+p1(2));
% C_low(2) = exp(p2(1)*HeatingRate+p2(2));
% C_up = [coeffACup,coeffBCup];



%--- original data from Wakili et al. 
 C_up  = [0.0873877575974930,336.990116244978];
 C_low = [177262519.522505,10405.0274800779];

% %%% --- trial_20 c/min
% C_up  = [25881.19013,19483.98907];
% C_low = [608206297.2,11568.9044];

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
    incident_HFlux = 60000;
   % heat_flux = e*(incident_HFlux-sigma*Temp(end)^4)+hc*(Tgas-Temp(end)); %W/m^2
%     heat_flux =incident_HFlux-e*sigma*(Temp(end)^4-Tgas^4)-hcE*(Temp(end)-Tgas);
    heat_flux =incident_HFlux-e*sigma*(Temp(end)^4-Tgas^4); %+hcE*(Tgas-Temp(end));
    %Temp_Fire = 800+273; %K
    Temp_E    = interp1q(Tfire(:,1)*60,Tfire(:,2),time)+273;
    %Temp_E = interp1q(hicks(:,1),hicks(:,2),time);% Fire side    [K];
    dTempdx_W = (hcW*(Temp_W-Temp(1))+e*sigma*(Temp_W^4-Temp(1)^4))./kgyp(1);
    % dTempdx_W = hc*(Temp_W-Temp(1))./kgyp(1);
    
    % Temp_W    = 20+273; %interp1q(Tfire(:,1)*60,Tfire(:,2),0)+273; % Ambient side [K];
    %dTempdx_E = (hc*(Temp_Fire-Temp(nx))+e*sigma*(Temp_Fire^4-Temp(nx)^4))./kgyp(nx);
    %dTempdx_E = (hc*(Temp_E-Temp(nx))+e*sigma*(Temp_E^4-Temp(nx)^4))./kgyp(nx);
    dTempdx_E = heat_flux/kgyp(nx);
    BcW       = 0; % BcE = 1, drichlet, BcE = 0, flux
    BcE       = 1; % BcE = 1, drichlet, BcE = 0, flux
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
    BcE       = 1; % BcE = 1, drichlet, BcE = 0, flux
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
%         rhogyp(i)    = interp1q(rho_gyp(:,1)+273,rho_gyp(:,2),Temp(i));
        CpEffgyp(i)  = interp1q(CpEff_gyp(:,1)+273,CpEff_gyp(:,2),Temp(i));
        HF(i)        = max(0,interp1q(Heat_Flow(:,1)+273,Heat_Flow(:,2)-20,Temp(i)));
        %          HF(i)       = interp1q(HF_NoLid(:,1)+273,HF_NoLid(:,2),Temp(i));
        kgyp(i)      = interp1q(k_gyp(:,1)+273,k_gyp(:,2),Temp(i));
%         W_loss_Rate(i) = interp1q(W_loss_rate(:,1),W_loss_rate(:,2),Temp(i));
        LHeat(i)       = interp1q(LatentHeat(:,1),LatentHeat(:,2),Temp(i));
        %       MLRS(i)        = interp1q(MLR(:,1),MLR(:,2),Temp(i));
%         MLRS(i)       = interp1q(MLRnew(:,1),MLRnew(:,2),Temp(i));
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
    
%% --- Here update heating rate to calculate Arhenious coefficients 
    HeatingRate = max(.1,((Temp - Temp_old)/dt));
% %     HeatingRate = 0.333*ones(1,length(Temp)); 
%     C_low1 = exp(p1(1)*HeatingRate+p1(2));
%     C_low2 = exp(p2(1)*HeatingRate+p2(2));
%      C_up = [coeffACup,coeffBCup];
    
% C_low1 = exp(P1(1)*HeatingRate+P1(2));
% C_low2 = exp(P2(1)*HeatingRate+P2(2));
% C_up1  = exp(P3(1)*HeatingRate+P3(2));
% C_up2  = exp(P4(1)*HeatingRate+P4(2));
    
    %%    
%     for jj = 1:length(Temp)
%            mask1 = Temp(jj)<=873;
%             SphiV(jj) = mask1*(C_low1(jj)*exp(-C_low2(jj)./Temp(jj)).*(H2o_st_1-Ssum(jj)))...
%                 +(1-mask1)*(C_up1(1)*exp(-C_up2(2)./(Temp(jj))).*(H2o_initial-Ssum(jj)));           
%     end
    
    for jj = 1:length(Temp)
            mask1 = Temp(jj)<=700;
            SphiV(jj) = mask1*(C_low(1)*exp(-C_low(2)./Temp(jj)).*(H2o_st_1-Ssum(jj)))...
                +(1-mask1)*(C_up(1)*exp(-C_up(2)./(Temp(jj))).*(H2o_initial-Ssum(jj)));
    end
%     
    
    
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
    
    %% Water re-condensation
    % determine the boiling point
    for i = 1:nx
        if Temp(i)<150+273
            Psat(i)    = interp1q(waterTable(:,1)+273,waterTable(:,2)*1000,Temp(i)); %Pa
            BP(i)      = interp1q(waterTable(:,2)/100,waterTable(:,1)+273,P(i));
            Rho_Sat(i) = Psat(i).*MWmix(i)./(Ru.*100*Temp(i)); % kg/m3
        else
            Rho_Sat(i) = phi_vapor(i);
        end
    end
    
    for i = 1:nx
        if Temp(i) < BP(i) && phi_vapor(i) > Rho_Sat(i)   
            hfg(i) = 1000*interp1q(waterTable(:,1)+273,waterTable(:,3),BP(i));
            %SphiV_recond(i) = abs(U_D(i))*(phi_vapor(i) - Rho_Sat(i))/dxu(i);
            SphiV_recond(i) = (phi_vapor(i) - Rho_Sat(i))/dt;
            SphiVnew(i)     = SphiV(i) - SphiV_recond(i);
            Q_recond(i)     = SphiV_recond(i) * hfg(i);   
        else
            Q_recond(i)     = 0;
            SphiV_recond(i) = 0;
            SphiVnew(i)     = SphiV(i);
        end
        SphiT(i) = STP(i)+QT(i)+Q_recond(i)/RhoCpeff(i);
    end
    
    % SphiT = STP+QT;
    
    %% ------------- solve for temperature ------------------------------------
    Temp = sclarTrans_eq(nx,dxp,dxu,alpha,dt,UvelTemp,const_1,const_2,alphaT,...
        id_scheme,Temp_old,BCT,SphiT);
%      Temp = Temp*(1-urelax)+ Temp*urelax;
    %% ----------- solve for species transpor ---------------------------------
    phi_vapor = sclarTrans_eq(nx,dxp,dxu,alpha,dt,U_D/Pore,const_1,const_2,Deff/Pore,...
        id_scheme,phi_vapor_old,BCV,SphiVnew);
    phi_vapor = phi_vapor*urelax+phi_vapor_old*(1-urelax);
    
    phi_air   = sclarTrans_eq(nx,dxp,dxu,alpha,dt,U_D/Pore,const_1,const_2,Deff/Pore,...
        id_scheme,phi_air_old,BCA,SphiA);
    phi_air = phi_air*urelax+phi_air_old*(1-urelax);
    subplot(3,1,1); plot(Temp); grid on 
    subplot(3,1,2); plot(phi_vapor); grid on 
    subplot(3,1,3); plot(phi_air); grid on 
    %% ---------- calculate pressure ------------------------------------------
    Yv     = phi_vapor./(phi_vapor+phi_air);
    Xv     = Yv.*MWair./(Yv*MWair+MWh2o-Yv*MWh2o);
    MWmix  = Xv*MWh2o + (1-Xv)*MWair; % Molecular weight of air-water mixture [gram]
    P      = phi_vapor.*(Ru.*Temp)./(Yv.*MWmix)/1000; % [bar]
    % P      = (phi_vapor+phi_air).*(Ru.*Temp)./(MWmix)/1000; % [bar]
    P(1)   = Pamb;
    P(end) = Pamb;
    
    P = P*urelax+ P_old*(1-urelax);
    
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
% load Tempexp
% figure('Name','Temperature profiles','Color',[1 1 1]);
% 
% % Create axes
% axes1 = axes('Position',...
%     [0.13 0.137543032775649 0.819464012251149 0.799194056426229]);
% hold(axes1,'on');
% 
% % Create multiple lines using matrix input to plot
% plot1 = plot(timeSave,TT(:,[60,40,20,1])-273,'LineWidth',3);
% set(plot1(1),'DisplayName','x = 0 mm, Numerical','Color',[0 0 1]);
% set(plot1(2),'DisplayName','x = 4 mm, Numerical','LineStyle','--',...
%     'Color',[1 0 0]);
% set(plot1(3),'DisplayName','x = 8 mm, Numerical','LineStyle','-.',...
%     'Color',[0 1 0]);
% set(plot1(4),'DisplayName','x = 12 mm, Numerical','LineWidth',2,...
%     'LineStyle',':');
% 
% timeR = timeR - 163;
% plot2 = plot(timeR(timeR>0 & timeR<timeSave(end)),termo_couples(timeR>0 & timeR<timeSave(end),[12,4,8,11]));
% set(plot2(1),'DisplayName','x = 0 mm, experimental','Color',[0 0 1]);
% set(plot2(2),'DisplayName','x = 4 mm, experimental','LineStyle','--',...
%     'Color',[1 0 0]);
% set(plot2(3),'DisplayName','x = 8 mm, experimental','LineStyle','-.',...
%     'Color',[0 1 0]);
% set(plot2(4),'DisplayName','x = 12 mm, experimental','LineWidth',2,...
%     'LineStyle',':');
%% For Report 07/28/2021
% Create figure
% % % % % % figure1 = figure('Name','Temperature profiles','Color',[1 1 1]);
% % % % % % 
% % % % % % % Create axes
% % % % % % axes1 = axes('Parent',figure1,...
% % % % % %     'Position',[0.13 0.137543032775649 0.819464012251149 0.799194056426229]);
% % % % % % hold(axes1,'on');
% % % % % % 
% % % % % % % Create multiple lines using matrix input to plot
% % % % % % plot1 = plot(timeSave,TT(:,[60,40,20,1])-273,'LineWidth',3,'Color',[0 0 0],'Parent',axes1);
% % % % % % set(plot1(1),'DisplayName','x = 0 mm, Numerical');
% % % % % % set(plot1(2),'DisplayName','x = 4 mm, Numerical','LineStyle','--');
% % % % % % set(plot1(3),'DisplayName','x = 8 mm, Numerical','LineStyle','-.');
% % % % % % set(plot1(4),'DisplayName','x = 12 mm, Numerical','LineStyle',':');
% % % % % % 
% % % % % % % Create multiple lines using matrix input to plot
% % % % % % plot2 = plot(timeR(timeR>0 & timeR<timeSave(end)),termo_couples(timeR>0 & timeR<timeSave(end),[13,4,8,11]),'LineWidth',3,'Color',[0 0 1],'Parent',axes1);
% % % % % % set(plot2(1),'DisplayName','x = 0 mm, experimental');
% % % % % % set(plot2(2),'DisplayName','x = 4 mm, experimental','LineStyle','--');
% % % % % % set(plot2(3),'DisplayName','x = 8 mm, experimental','LineStyle','-.');
% % % % % % set(plot2(4),'DisplayName','x = 12 mm, experimental','LineStyle',':');
% % % % % % 
% % % % % % % Create ylabel
% % % % % % ylabel('Temperature (^oC)','FontSize',14);
% % % % % % 
% % % % % % % Create xlabel
% % % % % % xlabel('Time (S)','FontSize',14);
% % % % % % 
% % % % % % box(axes1,'on');
% % % % % % % Set the remaining axes properties
% % % % % % set(axes1,'XGrid','on','XMinorTick','on','XTick',...
% % % % % %     [0:100:2000],...
% % % % % %     'YGrid','on','YMinorTick','on','YTick',...
% % % % % %     [0:50:1000]);
% % % % % % % Create legend
% % % % % % legend1 = legend(axes1,'show');
% % % % % % set(legend1,...
% % % % % %     'Position',[0.157065919642549 0.675890323353678 0.210441585955756 0.237552600020612],...
% % % % % %     'EdgeColor',[1 1 1]);
% % % % % % 
% % % % % % % ----- Save Data to excel for Shijin -------------------------------------
% % % % % % 
% % % % % % fname2 = [fname1(1:4), '_num_exp_data.xlsx'];
% % % % % % Time_num = timeSave; 
% % % % % % num_data = TT(:,[60,40,20,1])-273; 
% % % % % % Time_exp = timeR(timeR>0 & timeR<timeSave(end));
% % % % % % exp_data = termo_couples(timeR>0 & timeR<timeSave(end),[12,4,8,11]);
% % % % % % TableDatanum  = table(Time_num',num_data); 
% % % % % % TableDataexp  = table(Time_exp, exp_data); 
% % % % % % writetable(TableDatanum,fname2,'Sheet','numdata','Range','A1');
% % % % % % writetable(TableDataexp,fname2,'Sheet','expdata','Range','A1');






% % Create plot
% plot(T8mmexp(:,1),T8mmexp(:,2),'DisplayName','x = 8 mm, Experimental [8]',...
%     'MarkerEdgeColor',[0 1 0],...
%     'Marker','o',...
%     'LineWidth',2,...
%     'LineStyle','none',...
%     'Color',[0 1 0]);
% 
% % Create plot
% plot(T4mmexp(:,1),T4mmexp(:,2),'DisplayName','x = 4 mm, Experimental [8]','Marker','^',...
%     'LineWidth',2,...
%     'LineStyle','none',...
%     'Color',[1 0 0]);
% 
% % Create plot
% plot(T0mmexp(:,1),T0mmexp(:,2),'DisplayName','x =0 mm, Experimental [8]','Marker','x',...
%     'LineWidth',2,...
%     'LineStyle','none',...
%     'Color',[0 0 1]);
% 
% % Create ylabel
% ylabel('Temperature (^oC)','FontWeight','bold','FontName','Times New Roman');
% 
% % Create xlabel
% xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');
% 
% % Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 1000]);
% box(axes1,'on');
% % Set the remaining axes properties
% set(axes1,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
% % Create legend
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.133795525966869 0.57979351931209 0.450244687266684 0.355932193287349],...
%     'FontSize',13,...
%     'FontWeight','normal');


%% validation_ppt
load Tempexp
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
%% Time Histories of vapor density: ---> phi_vapor
figure4 = figure('Name','Water vapor profiles','Color',[1 1 1]);
axes1 = axes('Parent',figure4,...
    'Position',[0.157679626287364 0.160934883673525 0.77341578007306 0.781473240117578]);
hold(axes1,'on');

plot1 = plot(timeSave/60,phi_vapor_save(:,1:20:end),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm');
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--');
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.');
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');

ylabel('Vapor Desnity (Kg/m^3)','FontWeight','bold',...
    'FontName','Times New Roman');

xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');

box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.567667057133842 0.664373373014608 0.324646648164352 0.253143139180639],...
    'FontSize',16,...
    'EdgeColor',[1 1 1]);


%% Time Histories of air density: ---> phi_air
figure3 = figure('Name','Air profiles','Color',[1 1 1]);
axes1 = axes('Parent',figure3,...
    'Position',[0.157679626287364 0.151263701362736 0.747320373712636 0.781473242544421]);
hold(axes1,'on');
plot1 = plot(timeSave/60,phi_air_save(:,1:20:end),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','x = 0 mm');
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--');
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.');
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');
ylabel('Air Desnity (Kg/m^3)','FontWeight','bold',...
    'FontName','Times New Roman');
xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.564209207562917 0.652414268289552 0.324646648164352 0.253143139180639],...
    'FontSize',16,...
    'EdgeColor',[1 1 1]);



%% Transients of total pressure ---> PP
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
set(axes1,'FontSize',12,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.133890919319763 0.722842266561374 0.206484638078221 0.181919637801392]);
%% phi_vapor_propagation: plot the vapor density at different time
figure
time_m = [2.5*60/savetime 5*60/savetime 7.5*60/savetime 10*60/savetime 12.5*60/savetime 15*60/savetime];
figure
axes1 = axes('Position',...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182]);
hold(axes1,'on');

plot1 = plot(x*1000,phi_vapor_save(time_m,nx:-1:1),'LineWidth',2);
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
set(plot1(6),'DisplayName','t = 15 min',...
    'MarkerFaceColor',[0.301960796117783 0.745098054409027 0.933333337306976],...
    'MarkerSize',4,...
    'Marker','hexagram',...
    'LineStyle','--',...
    'LineWidth',0.5);

ylabel('Vapor Desnity (Kg/m^3)');

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

%%%%%%%%%%%%
%%%figure

time_m = [3*60/savetime 6*60/savetime 9*60/savetime 12*60/savetime 15*60/savetime 18*60/savetime];
xmm = x*1000;
phivapor1 = phi_vapor_save(time_m(1),:);
phivapor2 = phi_vapor_save(time_m(2),:);
phivapor3 = phi_vapor_save(time_m(3),:);
phivapor4 = phi_vapor_save(time_m(4),:);
phivapor5 = phi_vapor_save(time_m(5),:);
phivapor6 = phi_vapor_save(time_m(6),:);
%save('phivapor5kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6');
save('phivapor10kw.mat','xmm','phivapor1','phivapor2','phivapor3','phivapor4','phivapor5','phivapor6');

figure
plot(x*1000,phi_vapor_save(time_m(1),:))
hold on
plot(x*1000,phi_vapor_save(time_m(2),:))
plot(x*1000,phi_vapor_save(time_m(3),:))
plot(x*1000,phi_vapor_save(time_m(4),:))
plot(x*1000,phi_vapor_save(time_m(5),:))
plot(x*1000,phi_vapor_save(time_m(6),:))

figure
axes1 = axes('Position',...
    [0.13 0.161931818181818 0.806170212765957 0.763068181818182]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(x*1000,phi_vapor_save(time_m,nx:-1:1),'LineWidth',2);
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
set(plot1(6),'DisplayName','t = 18 min',...
    'MarkerFaceColor',[0.301960796117783 0.745098054409027 0.933333337306976],...
    'MarkerSize',4,...
    'Marker','hexagram',...
    'LineStyle','--',...
    'LineWidth',0.5);

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
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.289598115563532 0.718704518549598 0.485815593589705 0.20170453986661],...
    'NumColumns',2,...
    'FontSize',14,...
    'FontWeight','normal');

%  Temperature_no_exp

figure12 = figure('Name','Temperature profiles_no_exp','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure12,...
    'Position',[0.169484471131572 0.151078834588102 0.779548920784069 0.786369298191981]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(timeSave/60,TT(:,1:20:end)-273,'LineWidth',3);
set(plot1(1),'DisplayName','x = 0 mm','Color',[0 0 1]);
set(plot1(2),'DisplayName','x = 4 mm','LineStyle','--','Color',[1 0 0]);
set(plot1(3),'DisplayName','x = 8 mm','LineStyle','-.','Color',[0 1 0]);
set(plot1(4),'DisplayName','x = 12 mm','LineStyle',':');

% Create ylabel
ylabel('Temperature (^oC)','FontWeight','bold','FontName','Times New Roman');

% Create xlabel
xlabel('Time (min)','FontWeight','bold','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 1000]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',18);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.250001893604248 0.777662520627923 0.521968355183023 0.12136929126696],...
    'NumColumns',2,...
    'FontSize',18,...
    'EdgeColor',[0.972549021244049 0.972549021244049 0.972549021244049]);




%%%%%dehydration
% plot the source term with respect to distance at different time


figure
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),:)))
hold on
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),:)))
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),:)))
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),:)))
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),:)))
plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),:)))

set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),:))),'DisplayName','t = 2.5 min','LineWidth',2);
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),:))),'DisplayName','t = 5 min','LineWidth',3,'LineStyle','--');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),:))),'DisplayName','t = 7.5 min','LineWidth',3,'LineStyle',':');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),:))),'DisplayName','t = 10 min','LineWidth',3,'LineStyle','-.');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),:))),'DisplayName','t = 12.5 min','MarkerSize',4,'Marker','diamond',...
    'LineWidth',2,...
    'LineStyle','--');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),:))),'DisplayName','t = 15 min','MarkerSize',4,'Marker','+',...
    'LineWidth',2,...
    'LineStyle','--');



%%%%%%%%%%

%time_m = [2.5*60/dt 5*60/dt 7.5*60/dt 10*60/dt 12.5*60/dt 15*60/dt];
time_m = [3*60/dt 6*60/dt 9*60/dt 12*60/dt 15*60/dt 18*60/dt];
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
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(1),nx:-1:1))),'DisplayName','t = 3 min','LineWidth',2);
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(2),nx:-1:1))),'DisplayName','t = 6 min','LineWidth',3,'LineStyle','--');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(3),nx:-1:1))),'DisplayName','t = 9 min','LineWidth',3,'LineStyle',':');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(4),nx:-1:1))),'DisplayName','t = 12 min','LineWidth',3,'LineStyle','-.');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(5),nx:-1:1))),'DisplayName','t = 15 min','MarkerSize',4,'Marker','diamond',...
    'LineWidth',2,...
    'LineStyle','--');
set(plot(x*1000,(100/H2o_initial)*sum(Sum(1:time_m(6),nx:-1:1))),'DisplayName','t = 18 min','MarkerSize',4,'Marker','+',...
    'LineWidth',2,...
    'LineStyle','--');

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


%% Write data to excel for paper publication 

time_m = [3*60/savetime 6*60/savetime 9*60/savetime 12*60/savetime 15*60/savetime 18*60/savetime];
Time_num = timeSave; 
num_data = TT(:,[60,40,20,1])-273; 
% Time_exp = timeR(timeR>0 & timeR<timeSave(end));
% exp_data = termo_couples(timeR>0 & timeR<timeSave(end),[12,4,8,11]);
TableDatanum  = table(Time_num',num_data); 

TableDataexp  = table(T0mmexp,T4mmexp(1:24,:),T8mmexp(1:24,:)); 

writetable(TableDatanum,fname2,'Sheet','Tempnumdata','Range','A1');
writetable(TableDataexp,fname2,'Sheet','Tempexpdata','Range','A1');

TableVapDens = table(Time_num',phi_vapor_save(:,1:20:end));
TableAirDens = table(Time_num',phi_air_save(:,1:20:end));
TableTotPres = table(Time_num',PP(:,1:20:end));
TableVapDensProf = table([(x*1000)',phi_vapor_save(time_m,nx:-1:1)']);
for i = 1:5
AA(i,:) = (100/H2o_initial)*sum(Sum(1:time_m(i),nx:-1:1));
end

TableDepthDeHyration = table((x*1000)',AA');

writetable(TableVapDens,fname2,'Sheet','VaporDensity','Range','A1');
writetable(TableAirDens,fname2,'Sheet','AirDensity','Range','A1');
writetable(TableTotPres,fname2,'Sheet','totalPressure','Range','A1');
writetable(TableVapDensProf,fname2,'Sheet','VaporDensityProfile','Range','A1');
writetable(TableDepthDeHyration,fname2,'Sheet','Degre_Dehyd','Range','A1');
 














