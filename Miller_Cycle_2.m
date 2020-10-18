%% Miller Cycle - 8 states
%clear;
disp('Miller Cycle Calculations with Turbocharger')
disp('Please input your intake Temperature (K), Tin')
disp('Please input your intake Pressure (kPa), Pin')
disp('Here are Your State Values, Work, Power, Efficiency, and other values')

%% Constants
R = 0.287; %Ideal gas constant
cr = 9; %Compression Ratio 9:1
er = 10; %Expansion Raito 10:1
bts = 0.5; %Bore to Stroke 
ce = 0.90909; %Combustion Efficiency
Vd = 0.0015; %Displacement (m^3) V8=V1=V4=V5
Tin = 333; %Intake Temp. (K)
Pin = 101.325; %Intake Press. (kPa)
fa = 14.7; %Fuel to Air Ratio
n = 2; %Number of Cylinders
os = 5000; %Operating Speed in rpm
Cp = 1.005; %Specific Pressure of Air
Cv = 0.718; %Specific Volume of Air
hv = 45000; %Heating Value (kJ/kg)
C = 1.399721448; %Heat Capacity Ratio
Vc = Vd*0.11; %Clearence Volume (m^3) Vc=V2=V3=V6
nc = 0.9; %Combustion Efficiency
V7 = Vc*cr; %(m^3)

%% Temp.(K) & Press.(kPa) Cycle calculations 
%Need to add inefficincies to states 2 & 4
%The turbocharger inefficiencies are essentially the inefficiencies in
%states 5,6,7,8
T6 = 300; %K
P6 = 101.325; %kPa
T7 = Tin;
P7 = Pin;
P8 = P7;
T8 = T7*(V7/Vd)^(C-1); 
P1 = P8*(V7/Vd)^(C); %Vd=V8
T1 = T7*(V7/Vd)^(C-1); %Vd=V8
T2 = T7*(cr)^(C-1);
P2 = P7*(cr)^(C);

% Mass of gas in the cylinder
ma = (P1*Vd)/(R*T1); %mass of air in the cylinder
mf = ma/fa; %mass of fuel 

% Continuing Cycle calcuations 
T1 = T7*(V7/Vd)^(C-1); %Vd=V1
T2 = T7*(cr)^(C-1);
P2 = P7*(cr)^(C);
Qin = mf*hv*nc;
T3 = T2 + (Qin/(ma*Cv));
P3 = P2*(T3/T2);
T4 = T3*(Vc/Vd)^(C-1); %Vc=V3 & Vd=V4
P4 = P3*(1/er)^(C);
V4 = (ma*R*T4)/P4; %double check V4
P5 = 101; %100 kPa
T5 = T4*(P5/P4); 

%% Density Calculations (kg/m^3) Cycle Calculations
D1 = P1/(T1*R);
D2 = P1/(T2*R);
D3 = P3/(T3*R);
D4 = P4/(T4*R);
D5 = P5/(T5*R);
D6 = P6/(T6*R);
D7 = P7/(T7*R);
D8 = P8/(T8*R);

%% Work Cycle Calculations
%Work from State 3 to 4
W34 = (ma*R*(T4-T3))/(1-C);
%Work from State 7 to 2
W72 = (ma*R*(T2-T7))/(1-C);
%Work from State 6 to 7
W67 = P7*(V7-Vc); %Vc=V6
%Work from Sate 5 to 6 
W56 = P5*(Vc-Vd); %Vc=V6 & Vd=V5
%Net Work (kJ)
Wnet = (W34+W72+W67+W56);
%Work per cylinder (kJ/cylinder)
Wcylinder = Wnet/n;
%Work Specific (kJ/kg)
Wspec = Wnet/ma; 

%% Power Cycle Calculations
%Net Power (kW)
Pnet = (Wnet*os)/120;
%Power per cylinder (kW/cylinder)
Pcylinder = Pnet/n;
%Power Specific (kW/kg)
Pspec = Pnet/ma;
%Net Power (hp)
Php = (Pnet*1.34102);

%% Specific Fuel Consumption (SFC) (g/kW*hr)
SFC = (((C/fa)*ma)/Wnet)*3600000;

%% Efficiency (%)
Eff = abs(Wnet/Qin);

%% Indicated Mean Effecitve Pressure
Imep = Wnet/Vd; %(kPa)

%% Exhaust Temperature (K)
Pex = P5; %101 kPa
Tex = T4*(Pex/P4)^((C-1)/C);

%% Printing Values 
% Temperature (K)
fprintf('T1 = %6.3f\r',T1)
fprintf('T2 = %6.3f\r',T2)
fprintf('T3 = %6.3f\r',T3)
fprintf('T4 = %6.3f\r',T4)
fprintf('T5 = %6.3f\r',T5)
fprintf('T6 = %6.3f\r',T6)
fprintf('T7 = %6.3f\r',T7)
fprintf('T8 = %6.3f\r',T8)
% Presure (kPa)
fprintf('P1 = %6.3f\r',P1)
fprintf('P2 = %6.3f\r',P2)
fprintf('P3 = %6.3f\r',P3)
fprintf('P4 = %6.3f\r',P4)
fprintf('P5 = %6.3f\r',P5)
fprintf('P6 = %6.3f\r',P6)
fprintf('P7 = %6.3f\r',P7)
fprintf('P8 = %6.3f\r',P8)
%Density (kg/m^3)
fprintf('D1 = %6.3f\r',D1)
fprintf('D2 = %6.3f\r',D2)
fprintf('D3 = %6.3f\r',D3)
fprintf('D4 = %6.3f\r',D4)
fprintf('D5 = %6.3f\r',D5)
fprintf('D6 = %6.3f\r',D6)
fprintf('D7 = %6.3f\r',D7)
fprintf('D8 = %6.3f\r',D8)
%Qin
fprintf('Qin = %6.3f\r',Qin)
%Work 
fprintf('W34 = %6.3f\r',W34)
fprintf('W72 = %6.3f\r',W72)
fprintf('W67 = %6.3f\r',W67)
fprintf('W56 = %6.3f\r',W56)
fprintf('Wnet = %6.3f\r',Wnet)
fprintf('Wcylinder = %6.3f\r',Wcylinder)
fprintf('Wspec = %6.3f\r',Wspec)
%Power
fprintf('Pnet = %6.3f\r',Pnet)
fprintf('Pcylinder = %6.3f\r',Pcylinder)
fprintf('Pspec = %6.3f\r',Pspec)
fprintf('Php = %6.3f\r',Php)
%SFC
fprintf('SFC = %6.3f\r',SFC)
%Efficiency
fprintf('Eff = %6.3f\r',Eff)
%Indicated Mean Temperature
fprintf('Imep = %6.3f\r',Imep)
%Exhaust Temp.(K) & Press.(kPa)
fprintf('Tex = %6.3f\r',Tex)
fprintf('Pex = %6.3f\r',Pex)

