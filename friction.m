% program to compute friction mean effective pressure
% fmep units in kPa
% inputs
clear;
N = 7718; %engine speed rpm
b = 85.725; % bore (mm)
s = 88.9; % stroke (mm)
nc =3; % # cylinders
pin=101; % intake manifold pressure (kPa)

% everything below this line should be modified
%%
db = 56; % main bearing diameter (mm)
lb = 21; % main bearing length (mm)
niv = 2; % # intake valves/cyl
nev = 2; % # exhaust valves/cyl
div = 1.1875*25.4; % intake valve diameter (mm);
dev= 1.1875*25.4; % exhaust valve diameter (mm);
lv = 0.3995*25.4; % valve lift (mm)
mu = 100.e-3 ; % dynamic viscosity (Pa s)
pa= 101; %atmospheric pressure (kPa)
Up = 2.* N * s/60; % mean piston speed (mm/s)
denom = nc*b^2*s;
nb= nc+1; %# main crankshaft bearings
nv = (niv+nev)*nc; % # valves (total)
% friction coefficients
c_cb=0.0202; % crankshaft bearing
c_cs=93600; % crankshaft seals
c_pb=0.0202; % piston bearings
c_ps=14; % piston seals
c_pr=2707; % pison ringpack
c_vb=6720; % camshaft bearings
c_vh=0.5; % oscillating hydrodynamic
c_vm=10.7; % oscillating mixed
c_vs=1.2; % seals
c_vf= 207; % flat cam follower
c_vr=0.0151;% roller cam follower
c_1o=1.28; c_2o=0.0079; c_3o=-8.4e-7; %oil pump
c_1w=0.13; c_2w=0.002; c_3w=3.e-7; %water pump
c_1f=1.72; c_2f=0.00069; c_3f=1.2e-7; %fuel injection
c_iv=4.12e-3; % inlet valves (kPa sˆ2/mˆ2)
c_ev=c_iv; %exhaust valves (kPa sˆ2/mˆ2)
c_es=0.178; % exhaust system (kPa sˆ2/mˆ2)
% component fmeps
%crankshaft
f_cb=c_cb*nb*N.^(0.6)*db^3*lb/denom;
f_cs=c_cs*db/denom;
f_crank=f_cb+f_cs;
%piston assembly
f_pb=c_pb*nb*N.^(0.6)*db^3*lb/denom; % bearings
f_ps=c_ps*Up.^(0.5)/b; % skirt
f_pr=c_pr*Up.^(0.5)/(b^2); %ringpack
f_piston=f_pb+f_ps+f_pr;
%valvetrain
f_cam=c_vb*nb*N.^(0.6)/denom; %bearings
f_vh=c_vh*nv*lv^(1.5)*N.^(0.5)/denom; %oscill hydro
f_vm=c_vm*(2+10./(5+mu.*N))*lv*nv/(nc*s);
f_vs=c_vs; %seals
%cam followers - choose flat or roller
f_ff=c_vf*(2+10./(5+mu.*N))*nv/(nc*s); % flat
%f_rf=c_nv*N/(nc*s); % or roller
f_valve=f_cam+f_vh+f_vm+f_vs+f_ff;
%auxiliary
f_oil= c_1o + c_2o*N + c_3o*N.^2; %oil pump
f_wat= c_1w + c_2w*N + c_3w*N.^2; %water pump
f_fuel= c_1f + c_2f*N + c_3f*N.^2; %fuel pump
f_aux=f_oil+f_wat+f_fuel;
%pumping
dpis= pa-pin;
dpiv=c_iv*(pin/pa.*Up/1000*b^2/niv/div^2).^2;
dpev=c_ev*(pin/pa.*Up/1000*b^2/nev/dev^2).^2;
dpes=c_es*(pin/pa.*Up/1000).^2;
f_pump=dpis+dpiv+dpev+dpes;
%total
f_tot=f_crank+f_piston+f_valve+f_aux+f_pump;
fprintf(' \n fmep crankshaft (kPa)= %7.1f \n',f_crank);
fprintf(' fmep piston (kPa)= %7.1f \n',f_piston);
fprintf(' fmep valvetrain (kPa)= %7.1f \n',f_valve);
fprintf(' fmep auxiliary (kPa)= %7.1f \n',f_aux);
fprintf(' fmep pumping (kPa)= %7.1f \n',f_pump);
fprintf(' fmep total (kPa)= %7.1f \n',f_tot);