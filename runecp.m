%Input-Output program for ecp.m clear;
phi = 0.8; % enter equivalence ratio input
T = 3000; % enter temperature (K) input
P = 5000.; % enter pressure (kPa) input
fuel_id = 2;
% fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
% call ecp function
[ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T, P, phi, fuel_id );
%echo input
fprintf(' \n Equilibrium Combustion Solver \n' );
fprintf(' Pressure (kPa) = \t \t %6.1f \n', P );
fprintf(' Temperature (K) = \t \t %6.1f \n', T); ...
fprintf(' Fuel Air Equivalence ratio = \t% 3.1f \n ', phi);
%print output mole fractions and properties fprintf(' \n Mole Fractions \n' );
fprintf(' CO2 = \t %6.4f \n', Y(1) ); 
fprintf(' H2O = \t %6.4f \n', Y(2) ); 
fprintf(' N2 = \t %6.4f \n', Y(3) );
fprintf(' O2 = \t %6.4f \n', Y(4) ); 
fprintf(' CO = \t %6.4f \n', Y(5) ); 
fprintf(' H2 = \t %6.4f \n', Y(6) ); 
fprintf(' H = \t %6.4f \n', Y(7) ); 
fprintf(' O = \t %6.4f \n', Y(8) ); 
fprintf(' OH = \t %6.4f \n', Y(9) ); 
fprintf(' NO = \t %6.4f \n', Y(10) ); 
fprintf(' \n Mixture Properties \n' ); 
fprintf(' h(kJ/kg) = \t %6.1f \n', h ); 
fprintf(' u(kJ/kg) = \t %6.1f \n', u ); 
fprintf(' s (kJ/Kg K) = \t %6.3f \n', s ); 
fprintf(' v (m3/kg) = \t %6.3f \n', v ); 
fprintf(' cp (kJ/Kg K) =\t %6.3f \n', Cp ); 
fprintf(' Molecular Mass = %5.2f \n', MW ); 
fprintf(' dvdt = %8.2e \n', dvdT ); 
fprintf(' dvdp = %8.2e \n', dvdP );