%Input-Output program for running farg.m
clear;
T = 300; % enter temperature (K) input
P = 101.; % enter pressure (kPa) input
phi = 1; % enter equivalence ratio input
f=0.1111; %residual fraction input
fuel_id = 2;
% fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
% call farg function
[Y,h,u,s,v,R,Cp,MW,dvdT,dvdP] = farg(T,P,phi,f,fuel_id);

%[Y,h,u,s,v,R,Cp,MW,dvdT,dvdP,dMWdT,dMWdP]
%echo input
fprintf(' \n Fuel Air Residual Gas \n' );
fprintf(' Pressure (kPa) = %6.1f \n', P );
fprintf(' Temperature (K) = %6.1f \n', T); ...
fprintf(' Fuel Air Equivalence ratio = % 3.1f \n', phi);
fprintf(' Residual Fraction = \t% 3.1f \n ', f);
%print output mole fractions and properties
fprintf(' \n Mole Fractions \n' );
fprintf(' CO2 = \t %6.4f \n', Y(1) );
fprintf(' H2O = \t %6.4f \n', Y(2) );
fprintf(' N2 = \t %6.4f \n', Y(3) );
fprintf(' O2 = \t %6.4f \n', Y(4) );
fprintf(' CO = \t %6.4f \n', Y(5) );
fprintf(' H2 = \t %6.4f \n', Y(6) );
%fprintf(' H = \t %6.4f \n', Y(7) );
%fprintf(' O = \t %6.4f \n', Y(8) );
%fprintf(' OH = \t %6.4f \n', Y(9) );
%fprintf(' NO = \t %6.4f \n', Y(10) );
fprintf(' \n Mixture Properties \n' );
fprintf(' h(kJ/kg) = \t %6.1f \n', h );
fprintf(' u(kJ/kg) = \t %6.1f \n', u );
fprintf(' s (kJ/Kg K) = \t %6.3f \n', s );
fprintf(' v (m3/kg) = \t %6.3f \n', v );
fprintf(' Cp (kJ/Kg K) =\t %6.3f \n', Cp );
fprintf(' Molecular Mass = %5.2f \n', MW );
fprintf(' dvdt = %8.2e \n', dvdT );
fprintf(' dvdp = %8.2e \n', dvdP );