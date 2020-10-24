function [alpha,beta,gamma,delta,h,s,cp,mw,Fs,q ] = fuel( id, T )
% [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T )
%
% Parameters
% id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
% T - Temperature (K) at which to eval 300<T<1000 K
% Outputs
% alpha - # carbon
% beta - # hydrogen
% gamma - # oxygen
% delta - # nitrogen
% h - specific enthalpy (kJ/kg)
% s - specific entropy (kJ/kgK)
% cp - specific heat (kJ/kgK)
% mw - molecular weight (kg/kmol)
% Fs - stoichiometric fuel-air ratio
% q - heat of combustion (kJ/kg)
% Curve fit coefficients for thermodynamic properties of selected fuels
% a1 a2 a3 a6 a7
FuelProps = [ [ 1.971324, 7.871586e-3, -1.048592e-06, -9.930422e+3,8.873728 ]; ... % Methane
[ 4.0652, 6.0977e-2, -1.8801e-05, -3.588e+4,1.545e+1 ]; ... % Gasoline
[ 7.971, 1.1954e-01, -3.6858e-05, -1.9385e+4,-1.7879 ]; ... % Diesel
[ 1.779819, 1.262503e-02, -3.624890e-6,-2.525420e+4,1.50884e+1 ]; ... % Methanol
[ 1.412633, 2.0871e-02, -8.14213e-6, -1.02635e+4,1.917126e+1 ] ]; % Nitromethane
% Fuel chemical formula
% C H O N
% alpha beta gamma delta
FuelInfo = [[ 1 4 0 0 ]; ... % Methane
[ 7 17 0 0 ]; ... % Gasoline
[ 14.4 24.9 0 0 ]; ... % Diesel
[ 1 4 1 0 ]; ... % Methanol
[ 1 3 2 1 ] ]; % Nitromethane
% stoichiometric fuel-air ratio
FSv = [ 0.0584 0.06548 0.06907 0.1555 0.5924 ];
% available energy of combustion ac
ac = [ 52420 47870 45730 22680 12430 ];
% stoichiometric fuel-air ratio
Fs = FSv(id);
% available energy
q = ac(id);
% Get fuel composition
alpha = FuelInfo(id, 1);
beta = FuelInfo(id, 2);
gamma = FuelInfo(id, 3);
delta = FuelInfo(id, 4);
% compute fuel properties
ao = FuelProps(id, 1);
bo = FuelProps(id, 2);
co = FuelProps(id, 3);
do = FuelProps(id, 4);
eo = FuelProps(id, 5);
% compute thermodynamic properties
h = ao + bo/2*T +co/3*T^2 +do/T;
s = ao*log(T) + bo*T +co/2*T^2 + eo;
cp = ao + bo*T + co*T^2;
% Calculate molecular weight of fuel
mw = 12.01*alpha + 1.008*beta + 16.00*gamma + 14.01*delta;
fprintf('Fuel Composition \n')
fprintf('Carbon = \t %6.1f \n', alpha)
fprintf('Hydrogen = \t %6.1f \n', beta)
fprintf('Oxygen = \t %6.1f \n', gamma)
fprintf('Nitrogen = \t %6.1f \n', delta)
end