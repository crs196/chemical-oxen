function [Y,h,u,s,v,R,Cp,MW,dvdT,dvdP]=farg(T,P,phi,f,fuel_id)
% Subroutine for Fuel Air Residual Gas
%
% inputs:
% T - temperature (K) [ 300 --> 1000 K ]
% P - pressure (kPa)
% phi - equivalence ratio = actual fuel-to-air ratio/stoichimetric
% fuel-air ratio
% stoich fuel-air ratio: 14.7
% f - residual fraction
% fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
%
% outputs:
% y - mole fraction of constituents
% y(1) : CO2
% y(2) : H2O
% y(3) : N2
% y(4) : O2
% y(5) : CO
% y(6) : H2
% h - specific enthalpy of mixture, kJ/kg
% u - specific internal energy of mixture, kJ/kg
% s - specific entropy of mixture, kJ/kgK
% v - specific volume of mixture, m3/kg
% r - specific ideal gas constant, kJ/kgK
% cp - specific heat at constant pressure, kJ/kgK
% mw - molecular weight of mixture, kg/kmol
% dvdt - (dv/dT) at const P, m3/kg per K
% dvdp - (dv/dP) at const T, m3/kg per kPa
% Get fuel composition information
[alpha,beta,gamma,delta,h_fuel,so_fuel,cp_fuel,m_fuel ]=fuel( fuel_id, T );
% Curve fit coefficients for thermodynamic properties
% 300 < T < 1000 K
% Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
% h/RT = a1 + a2/2*T + a3/3*T^2 + a4/4*T^3 + a5/5*T^4 + a6/T
% so/R = a1*ln(T) + a2*T + a3/2*T^2 + a4/3*T^3 + a5/4*T^4 + a7
A = [ [ 0.24007797e+1, 0.87350957e-2, -0.66070878e-5, 0.20021861e-8,0.63274039e-15, -0.48377527e+5, 0.96951457e+1 ]; ... % CO2
    [ 0.40701275e+1, -0.11084499e-2, 0.41521180e-5, -0.29637404e-8,0.80702103e-12, -0.30279722e+5, -0.32270046 ]; ... % H2O
    [ 0.36748261e+1, -0.12081500e-2, 0.23240102e-5, -0.63217559e-9, -0.2257725e-12, -0.10611588e+4, 0.23580424e+1 ]; ... % N2
    [ 0.36255985e+1, -0.18782184e-2, 0.70554544e-5, -0.67635137e-8, 0.21555993e-11, -0.10475226e+4, 0.43052778e+1 ]; ... % O2
    [ 0.37100928e+1, -0.16190964e-2, 0.36923594e-5, -0.20319674e-8, 0.23953344e-12, -0.14356310e+5, 0.2955535e+1 ]; ... % CO
    [ 0.30574451e+1, 0.26765200e-2, -0.58099162e-5, 0.55210391e-8, -0.1812273e-11, -0.98890474e+3, -0.22997056e+1 ] ]; % H2
% molar mass of constituents
% CO2 H2O N2 O2 CO H2
Mi = [ 44.01, 18.02, 28.013, 32.00, 28.01, 2.016 ];
% Calculate stoichiometric molar air-fuel ratio
a_s = alpha + beta/4 - gamma/2
% mole fraction of fuel, O2, N2
y_1 = 1 / (1 + 4.76*a_s/phi); % mole fraction for one mole of reactant
y_fuel = y_1; % assuming 1 mole fuel
y_O2 = a_s/phi * y_1; % a_s/phi moles O2
y_N2 = a_s/phi*3.76 * y_1; % a_s/phi * 3.76 moles N2
% mass of fuel air mixture (M’’)
m_fa = y_fuel*m_fuel + y_O2*32.00 + y_N2*28.013;
% default case: no residual gas
Y = zeros(6,1);
m_r = 0; % mass of residual gas
y_r = 0; % mole fraction of residual gas in mixture
n = zeros(6,1);
dcdt = 0;
if ( phi <= 1 )
    % lean combustion
    n(1) = alpha;
    n(2) = beta/2;
    n(3) = delta/2 + 3.76*a_s/phi;
    n(4) = a_s*(1/phi - 1);
else
    % rich combustion
    d1 = 2*a_s*(1-1/phi);
    z = T/1000;
    K = exp( 2.743 - 1.761/z - 1.611/z^2 + 0.2803/z^3 );
    a1 = 1-K;
    b1 = beta/2 + alpha*K - d1*(1-K);
    c1 = -alpha*d1*K;
    n(5) = (-b1 + sqrt(b1^2 - 4*a1*c1))/(2*a1);
    % Required derivatives for Cp calculation of mixture
    % calculate dcdt = dn5/dK * dK/dT
    dkdt = -K*(-1.761+z*(-3.222+z*.8409))/1000;
    dn5dk = -((alpha - n(5))*(n(5) + 2*a_s*(1/phi - 1)))/(beta/2 + n(5)+ 2*a_s*(1/phi - 1));
    dcdt = dn5dk * dkdt;
    n(1) = alpha - n(5);
    n(2) = beta/2 - d1 + n(5);
    n(3) = delta/2 + 3.76*a_s/phi;
    n(6) = d1 - n(5);
end
% total moles
N = sum(n);
% calculate mole fractions and mass of residual gas
m_r = 0;
for i=1:6
    Y(i) = n(i)/N;
    m_r = m_r + Y(i)*Mi(i);
end
% compute residual mole fraction
y_r = 1/(1 + m_r/m_fa * (1/f-1));
% compute total mole fractions in mixture
for i=1:6
    Y(i) = Y(i)*y_r;
end
% fuel mole fraction based on all moles
y_fuel = y_fuel*(1 - y_r);
% include intake N2 and O2
Y(3) = Y(3) + y_N2*(1 - y_r);
Y(4) = Y(4) + y_O2*(1 - y_r);
% compute properties of mixture
h = h_fuel*y_fuel;
s = (so_fuel-log(max(y_fuel,1e-15)))*y_fuel;
Cp = cp_fuel*y_fuel;
MW = m_fuel*y_fuel;
% compute component properties according to curve fits
cpo = zeros(6,1);
ho = zeros(6,1);
so = zeros(6,1);
for i=1:6
    cpo(i) = A(i,1) + A(i,2)*T + A(i,3)*T^2 + A(i,4)*T^3 + A(i,5)*T^4;
    ho(i) = A(i,1) + A(i,2)/2*T + A(i,3)/3*T^2 + A(i,4)/4*T^3+A(i,5)/5*T^4 + A(i,6)/T;
    so(i) = A(i,1)*log(T) + A(i,2)*T + A(i,3)/2*T^2 + A(i,4)/3*T^3+A(i,5)/4*T^4 +A(i,7);
end
table = [-1,1,0,0,1,-1];
for i=1:6
    if(Y(i)>1.e-25)
        h = h + ho(i)*Y(i);
        s = s + Y(i)*(so(i)-log(Y(i)));
        Cp = Cp+cpo(i)*Y(i)+ho(i)*T*table(i)*dcdt*y_r/N;
        MW = MW + Y(i)*Mi(i);
    end
end
% compute thermodynamic properties
R = 8.31434/MW; % compute mixture gas constant
h = R*T*h; % curve fit for h is h/rt
u = h-R*T;
v = R*T/P;
s = R*(-log(P/101.325)+s);
Cp = R*Cp; % curve fit for cp is cp/r
dvdT = v/T; % derivative of volume wrt temp
dvdP = -v/P; % derivative of volume wrt pres