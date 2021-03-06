%% lets try to run both runecp.m and runfarg.m from the same script!
% this is the run file for ecp, farg, fuel, and crankanglecalcs_new
close all
clear
clc

phi = 0.765; % enter equivalence ratio input
T = 3000; % enter temperature (K) input
P = 5000.; % enter pressure (kPa) input
thetas=-6; % start of heat release (deg)
thetad=39.4; % duration of heat release (deg)
rc=9; %compression ratio
re=10; %expansion ratio
a= 5; %weibe parameter a
n= 3; %weibe exponent n
b=.085725; %bore (m)
stroke=.0889; %stroke (m)
st=8.89; %stroke (cm)
len= 13.35; %connecting rod length (cm)

fuel_id = 2;
if T<1000
    %Input-Output program for running farg.m clear;
    f=0.1111; %residual fraction input
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
    [T]=CrankAngleCalcs(rc,re,thetas,thetad,a,n,b,stroke,len,Cp,MW,h,R,v,P);

elseif (T>600 || T<3500) || (P>20 || P<30000)
    
    % fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
    % call ecp function
    [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP,dMdT,dMdP] = ecp( T, P, phi, fuel_id );
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
    [T]=CrankAngleCalcs(rc,re,thetas,thetad,a,n,b,stroke,len,Cp,MW,dMdT,h,R,dMdP,v,P);
else
    fprintf('Welp try a diff temp or pressure value \n')
end
%% Gamma as a function of crank angle
%[T]=CrankAngleCalcs_jack(rc,re,thetas,thetad,a,n,b,stroke,len,Cp,MW,hate,h,R);
