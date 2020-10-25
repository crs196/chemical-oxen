function [ierr,Y,h,u,s,v,R,Cp,MW,dvdT,dvdP]=ecp(T,P,phi,ifuel ) 
% Subroutine for Equilibrium Combustion Products
%
% inputs:
% T - temperature (K) [ 600 --> 3500 ]
% P - pressure (kPa) [ 20 --> 30000 ]
% phi - equivalence ratio [ 0.01 --> ~3 ]
% ifuel - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane 
%
% outputs:
%   ierr - Error codes:
%   0 = success
%   1 = singular matrix
%   2 = maximal pivot error in gaussian elimination
%   3 = no solution in maximum number of iterations
%   4 = result failed consistency check sum(Y)=1
%   5 = failure to obtain initial guess for oxygen concentration
%   6 = negative oxygen concentration in initial guess calculation
%   7 = maximum iterations reached in initial guess solution
%   8 = temperature out of range
%   9 = pressure out of range
%   10 = equivalence ratio too lean
%   11 = eqivalence ratio too rich, solid will be formed
%
%   y - mole fraction of constituents
%   y(1)    : CO2
%   y(2)    : H2O
%   y(3)    : N2
%   y(4)    : O2
%   y(5)    : CO
%   y(6)    : H2
%   y(7)    : H
%   y(8)    : O
%   y(9)    : OH
%   y(10)   : NO
%   h - specific enthalpy of mixture, kJ/kg
%   u - internal energy of mixture, kJ/kg
%   s - specific entropy of mixture, kJ/kgK
%   v - specific volume of mixture, m3/kg
%   R - specific ideal gas constant, kJ/kgK
%   Cp - specific heat at constant, kJ/kgK
%   MW - molecular weight of mixture, kg/kmol
%   dvdt - (dv/dT) at const P, m3/kg per K
%   dvdp - (dv/dP) at const T, m3/kg per kPa

% initialize outputs
Y = zeros(10,1);
h = 0;
u = 0;
s = 0;
v = 0;
R = 0;
Cp = 0;
MW = 0;
dvdT = 0;
dvdP = 0;
% solution parameters 
prec = 1e-3;
MaxIter = 20;
% square root of pressure (used many times below) 
PATM = P/101.325;
sqp = sqrt(PATM);
if ( T < 600 || T > 3500 )
    ierr = 8;
    return; 
end
if ( P < 20 || P > 30000 ) 
    ierr = 9;
    return; 
end
if ( phi < 0.01 ) 
    ierr = 10;
    return; 
end
% Get fuel composition information
[ alpha, beta, gamma, delta ] = fuel( ifuel, T );
% Equilibrium constant curve fit coefficients.
% Valid in range: 600 K < T < 4000 K
%          Ai            Bi           Ci           Di             Ei
Kp = [ [ 0.432168, -0.112464e+5, 0.267269e+1, -0.745744e-4, 0.242484e-8 ];...
    [ 0.310805, -0.129540e+5, 0.321779e+1, -0.738336e-4, 0.344645e-8 ]; ...
    [ -0.141784, -0.213308e+4, 0.853461, 0.355015e-4, -0.310227e-8 ]; ...
    [ 0.150879e-1, -0.470959e+4, 0.646096, 0.272805e-5, -0.154444e-8 ]; ...
    [ -0.752364, 0.124210e+5, -0.260286e+1, 0.259556e-3, -0.162687e-7 ]; ...
    [ -0.415302e-2, 0.148627e+5, -0.475746e+1, 0.124699e-3, -0.900227e-8 ] ];

    
K = zeros(6,1);
for i=1:6
log10ki = Kp(i,1)*log(T/1000) + Kp(i,2)/T + Kp(i,3) + Kp(i,4)*T + ... 
Kp(i,5)*T*T;
K(i) = 10^log10ki;
end
c1 = K(1)/sqp;
c2 = K(2)/sqp;
c3 = K(3);
c4 = K(4);
c5 = K(5)*sqp;
c6 = K(6)*sqp;
[ierr,y3,y4,y5,y6 ] = guess( T, phi, alpha, beta, gamma, delta, c5, c6 ); 
if ( ierr  ~= 0 )
    return; 
end
a_s = alpha + beta/4 - gamma/2;
D1 = beta/alpha;
D2 = gamma/alpha + 2*a_s/(alpha*phi);
D3 = delta/alpha + 2*3.7619047619*a_s/(alpha*phi); 
A = zeros(4,4);
final = 0;
for jj=1:MaxIter
    sqy6 = sqrt(y6); 
    sqy4 = sqrt(y4); 
    sqy3 = sqrt(y3); 
    y7= c1*sqy6;
    y8= c2*sqy4;
    y9= c3*sqy4*sqy6; 
    y10= c4*sqy4*sqy3; 
    y2= c5*sqy4*y6; 
    y1= c6*sqy4*y5;
    d76 = 0.5*c1/sqy6;
    d84 = 0.5*c2/sqy4;
    d94 = 0.5*c3*sqy6/sqy4; 
    d96 = 0.5*c3*sqy4/sqy6; 
    d103 = 0.5*c4*sqy4/sqy3; 
    d104 = 0.5*c4*sqy3/sqy4; 
    d24 = 0.5*c5*y6/sqy4; 
    d26 = c5*sqy4;
    d14 = 0.5*c6*y5/sqy4; 
    d15 = c6*sqy4;
    % form the Jacobian matrix
A = [ [ 1+d103, d14+d24+1+d84+d104+d94, d15+1, d26+1+d76+d96 ]; ... 
    [ 0, 2.*d24+d94-D1*d14, -D1*d15-D1, 2*d26+2+d76+d96; ]; ...
    [ d103, 2*d14+d24+2+d84+d94+d104-D2*d14,2*d15+1-D2*d15-D2, d26+d96 ]; ...
    [ 2+d103, d104-D3*d14, -D3*d15-D3,0 ] ];
if ( final )
    break; 
end
B = [ -(y1+y2+y3+y4+y5+y6+y7+y8+y9+y10-1); ...
    -(2.*y2 + 2.*y6 + y7 + y9 -D1*y1 -D1*y5); ...
    -(2.*y1 + y2 +2.*y4 + y5 + y8 + y9 + y10 -D2*y1 -D2*y5); ... 
    -(2.*y3 + y10 -D3*y1 -D3*y5) ];
[ B, ierr ] = gauss( A, B ); 
if ( ierr  ~= 0 )
    return;
end
y3 = y3 + B(1); 
y4 = y4 + B(2); 
y5 = y5 + B(3); 
y6 = y6 + B(4); 
nck = 0;
if ( abs(B(1)/y3) > prec ) 
    nck = nck+1;
end
if ( abs(B(2)/y4) > prec )
    nck = nck+1; 
end
if ( abs(B(3)/y5) > prec ) 
    nck = nck+1;
end
if ( abs(B(4)/y6) > prec )
    nck = nck+1; 
end
if( nck == 0 )
    % perform top half of loop to update remaining mole fractions 
    % and Jacobian matrix
    final = 1;
    continue;
end
end
if (jj>=MaxIter)
    ierr = 3;
    return; 
end
Y = [ y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 ];
% consistency check
if( abs( sum(Y)-1 ) > 0.0000001 )
    ierr = 4;
    return; 
end
% constants for partial derivatives of properties
dkdt = zeros(6,1);
for i=1:6
    dkdt(i)=2.302585*K(i)*( Kp(i,1)/T - Kp(i,2)/(T*T)+ Kp(i,4)+2*Kp(i,5)*T ); 
end
dcdt = zeros(6,1);
dcdt(1) = dkdt(1)/sqp;
dcdt(2) = dkdt(2)/sqp;
dcdt(3) = dkdt(3);
dcdt(4) = dkdt(4);
dcdt(5) = dkdt(5)*sqp;
dcdt(6) = dkdt(6)*sqp;
dcdp = zeros(6,1);
dcdp(1) = -0.5*c1/P;
dcdp(2) = -0.5*c2/P;
dcdp(5) = 0.5*c5/P;
dcdp(6) = 0.5*c6/P;
x1 = Y(1)/c6;
x2 = Y(2)/c5;
x7 = Y(7)/c1;
x8 = Y(8)/c2;
x9 = Y(9)/c3;
x10 = Y(10)/c4;
dfdt(1) = dcdt(6)*x1 + dcdt(5)*x2 + dcdt(1)*x7 + dcdt(2)*x8 + dcdt(3)*x9 + dcdt(4)*x10;
dfdt(2) = 2.*dcdt(5)*x2 + dcdt(1)*x7 + dcdt(3)*x9 -D1*dcdt(6)*x1;
dfdt(3) = 2.*dcdt(6)*x1+dcdt(5)*x2+dcdt(2)*x8+dcdt(3)*x9+dcdt(4)*x10 - D2*dcdt(6)*x1;
dfdt(4) = dcdt(4)*x10 -D3*dcdt(6)*x1;
dfdp(1) = dcdp(6)*x1 + dcdp(5)*x2 + dcdp(1)*x7 +dcdp(2)*x8;
dfdp(2) = 2.*dcdp(5)*x2 + dcdp(1)*x7 -D1*dcdp(6)*x1;
dfdp(3) = 2.*dcdp(6)*x1 + dcdp(5)*x2 + dcdp(2)*x8 - D2*dcdp(6)*x1;
dfdp(4) = -D3*dcdp(6)*x1;
dfdphi(1) = 0;
dfdphi(2) = 0;
dfdphi(3) = 2*a_s/(alpha*phi*phi)*(Y(1)+Y(5));
dfdphi(4) = 2*3.7619047619*a_s/(alpha*phi*phi)*(Y(1)+Y(5));
% solve matrix equations for independent temperature derivatives 
b = -1.0 .* dfdt'; %element by element mult.
[b, ierr] = gauss(A,b);% solve for new b with t derivatives
if ( ierr  ~= 0 )
    return; 
end
dydt(3) = b(1);
dydt(4) = b(2);
dydt(5) = b(3);
dydt(6) = b(4);
dydt(1) = sqrt(Y(4))*Y(5)*dcdt(6) + d14*dydt(4) + d15*dydt(5); dydt(2) = sqrt(Y(4))*Y(6)*dcdt(5) + d24*dydt(4) + d26*dydt(6); dydt(7) = sqrt(Y(6))*dcdt(1) + d76*dydt(6);
dydt(8) = sqrt(Y(4))*dcdt(2) + d84*dydt(4);
dydt(9) = sqrt(Y(4)*Y(6))*dcdt(3) + d94*dydt(4) + d96*dydt(6); dydt(10) = sqrt(Y(4)*Y(3))*dcdt(4) + d104*dydt(4) + d103*dydt(3); % solve matrix equations for independent pressure derivatives
b = -1.0 .* dfdp'; %element by element mult.
[b,ierr] = gauss(A,b); % solve for new b with p derivatives
if ( ierr ~=0 )
    return; 
end
dydp(3) = b(1);
dydp(4) = b(2);
dydp(5) = b(3);
dydp(6) = b(4);
dydp(1) = sqrt(Y(4))*Y(5)*dcdp(6) + d14*dydp(4) + d15*dydp(5); dydp(2) = sqrt(Y(4))*Y(6)*dcdp(5) + d24*dydp(4) + d26*dydp(6); dydp(7) = sqrt(Y(6))*dcdp(1) + d76*dydp(6);
dydp(8) = sqrt(Y(4))*dcdp(2) + d84*dydp(4);
dydp(9) = d94*dydp(4) + d96*dydp(6);
dydp(10)= d104*dydp(4) + d103*dydp(3);
% molecular weights of constituents (g/mol)
%       CO2    H2O     N2      O2     CO     H2     H     O     OH      NO
Mi = [ 44.01, 18.02, 28.013, 32.00, 28.01, 2.016, 1.009, 16., 17.009, 30.004];
if ( T > 1000 )
    % high temp curve fit coefficients for thermodynamic properties ... 1000 < T < 3000 K
    AAC = [ ...
        [.446080e+1,.309817e-2,-.123925e-5,.227413e-9, -.155259e-13, ...
            -.489614e+5,-.986359 ];
        [.271676e+1,.294513e-2,-.802243e-6,.102266e-9, -.484721e-14, ... 
            -.299058e+5,.663056e+1 ];
        [.289631e+1,.151548e-2,-.572352e-6,.998073e-10,-.652235e-14, ... 
            -.905861e+3,.616151e+1 ];
        [.362195e+1,.736182e-3,-.196522e-6,.362015e-10,-.289456e-14,... 
            -.120198e+4,.361509e+1 ];
        [.298406e+1,.148913e-2,-.578996e-6,.103645e-9, -.693535e-14,... 
            -.142452e+5,.634791e+1 ];
        [.310019e+1,.511194e-3, .526442e-7,-.349099e-10,.369453e-14,... 
            -.877380e+3,-.196294e+1 ];
        [.25e+1,0,0,0,0,.254716e+5,-.460117 ]; 
        [.254205e+1,-.275506e-4,-.310280e-8,.455106e-11,-.436805e-15,...
            .292308e+5,.492030e+1 ]; 
        [.291064e+1,.959316e-3,-.194417e-6,.137566e-10,.142245e-15,...
            .393538e+4,.544234e+1 ];
        [.3189e+1 ,.133822e-2,-.528993e-6,.959193e-10,-.648479e-14,...
            .982832e+4,.674581e+1 ]; ];
else
    % low temp curve fit coefficients for thermodynamic properties, 300 < T <= 1000 K 
    AAC = [ ...
        [ 0.24007797e+1, 0.87350957e-2, -0.66070878e-5, 0.20021861e-8, ... 
            0.63274039e-15, -0.48377527e+5, 0.96951457e+1 ]; % CO2
        [ 0.40701275e+1, -0.11084499e-2, 0.41521180e-5, -0.29637404e-8, ...
            0.80702103e-12, -0.30279722e+5, -0.32270046 ]; %H2O
        [ 0.36748261e+1, -0.12081500e-2, 0.23240102e-5, -0.63217559e-9, ...
            -0.22577253e-12, -0.10611588e+4, 0.23580424e+1 ]; % N2
        [ 0.36255985e+1, -0.18782184e-2, 0.70554544e-5, -0.67635137e-8, ... 
            0.21555993e-11, -0.10475226e+4, 0.43052778e+1 ]; % O2
        [ 0.37100928e+1, -0.16190964e-2, 0.36923594e-5, -0.20319674e-8, ... 
            0.23953344e-12, -0.14356310e+5, 0.2955535e+1 ]; % CO
        [ 0.30574451e+1, 0.26765200e-2, -0.58099162e-5, 0.55210391e-8, ... 
            -0.18122739e-11, -0.98890474e+3, -0.22997056e+1 ]; % H2
        [ 0.25000000e+1, 0, 0, 0, 0, 0.25471627e+5, ... 
            -0.46011762e+0 ]; % H
        [ 0.29464287e+1, -0.16381665e-2, 0.24210316e-5, -0.16028432e-8, ... 
            0.38906964e-12, 0.29147644e+5, 0.29639949e+1 ]; % O
        [ 0.38375943e+1, -0.10778858e-2, 0.96830378e-6, 0.18713972e-9, .... 
            -0.22571094e-12, 0.36412823e+4, 0.49370009e+0 ]; % OH
        [ 0.40459521e+1, -0.34181783e-2, 0.79819190e-5, -0.61139316e-8, ... 
            0.15919076e-11, 0.97453934e+4, 0.29974988e+1 ]; % H2
        ];
end

% Compute cp,h,s
% initialize h, etc to zero 
MW = 0;
Cp = 0;
h = 0;
s = 0; 
dMWdT = 0; 
dMWdP = 0; 
for i=1:10
    cpo = AAC(i,1) + AAC(i,2)*T + AAC(i,3)*T^2 + AAC(i,4)*T^3 + AAC(i,5)*T^4;
    ho = AAC(i,1) + AAC(i,2)/2*T + AAC(i,3)/3*T^2 + AAC(i,4)/4*T^3 + ... 
        AAC(i,5)/5*T^4 + AAC(i,6)/T;
    so = AAC(i,1)*log(T) + AAC(i,2)*T + AAC(i,3)/2*T^2 + AAC(i,4)/3*T^3 + ...
        AAC(i,5)/4*T^4 +AAC(i,7);
        h = h + ho*Y(i); % h is h/rt here
        MW = MW + Mi(i)*Y(i);
        dMWdT = dMWdT + Mi(i)*dydt(i); 
        dMWdP = dMWdP + Mi(i)*dydp(i);
        Cp = Cp+Y(i)*cpo + ho*T*dydt(i);
        if (Y(i)> 1.0e-37)
            s = s + Y(i)*(so - log(Y(i))); 
        end
end
        
R = 8.31434/MW;
v = R*T/P;
Cp = R*(Cp - h*T*dMWdT/MW);
h = h*R*T;
s = R*(-log(PATM) + s);
u=h-R*T;
dvdT = v/T*(1 - T*dMWdT/MW);
dvdP = v/P*(-1 + P*dMWdP/MW);
ierr = 0;
return;

function [ierr,y3,y4,y5,y6] = guess(T,phi,alpha,beta,gamma,delta,c5,c6) 
        ierr = 0;
        y3 = 0;
        y4 = 0;
        y5 = 0;
        y6 = 0;
        % estimate number of total moles produced, N 
        n = zeros(6,1);
        % Calculate stoichiometric molar air-fuel ratio 
        a_s = alpha + beta/4 - gamma/2;
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
            KK = exp( 2.743 - 1.761/z - 1.611/z^2 + 0.2803/z^3 ); 
            aa = 1-KK;
            bb = beta/2 + alpha*KK - d1*(1-KK);
            cc = -alpha*d1*KK;
            n(5) = (-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa);
            n(1) = alpha - n(5);
            n(2) = beta/2 - d1 + n(5);
            n(3) = delta/2 + 3.76*a_s/phi;
            n(6) = d1 - n(5);
        end
        % total product moles per 1 mole fuel
        N = sum(n);
        % try to get close to a reasonable value of ox mole fraction 
        % by finding zero crossing of ’f’ function
        ox = 1;
        nIterMax=40;
        for ii=1:nIterMax
            f = 2*N*ox - gamma - (2*a_s)/phi + (alpha*(2*c6*ox^(1/2) + 1))/ ... 
                (c6*ox^(1/2) + 1) + (beta*c5*ox^(1/2))/(2*c5*ox^(1/2) + 2);
            if ( f < 0 )
                break;
            else
                ox = ox*0.1;
                if ( ox < 1e-37 ) 
                    ierr = 5;
                    return; 
                end
            end
        end
% now zero in on the ox mole fraction using Newton-Raphson iteration 
for ii=1:nIterMax
    f = 2*N*ox - gamma - (2*a_s)/phi + (alpha*(2*c6*ox^(1/2) + 1))/ ... 
        (c6*ox^(1/2) + 1) + (beta*c5*ox^(1/2))/(2*c5*ox^(1/2) + 2);
    df = 2*N - (beta*c5^2)/(2*c5*ox^(1/2) + 2)^2 + (alpha*c6)/ ... 
        (ox^(1/2)*(c6*ox^(1/2) + 1)) +(beta*c5)/ ... 
        (2*ox^(1/2)*(2*c5*ox^(1/2) + 2)) - (alpha*c6*(2*c6*ox^(1/2) + 1))/ ... 
        (2*ox^(1/2)*(c6*ox^(1/2) + 1)^2);
            dox = f/df;
            ox = ox - dox; 
            if ( ox < 0.0 )
                ierr = 6;
                return; 
            end
            if ( abs(dox/ox) < 0.001 ) 
                break;
            end
end
        if( ii == nIterMax ) 
            ierr = 7;
            return; 
        end
        y3 = 0.5*(delta + a_s/phi*2*3.76)/N; 
        y4 = ox;
        y5 = alpha/N/(1+c6*sqrt(ox));
        y6 = beta/2/N/(1+c5*sqrt(ox));
end % guess

function [B, IERQ] = gauss( A, B )
% maximum pivot gaussian elimination routine adapted
% from FORTRAN in Olikara &RBorman, SAE 750468, 1975
% not using built-in MATLAB routines because they issue % lots of warnings for close to singular matrices
% that haven’t seemed to cause problems in this application % routine below does check however for true singularity
    IERQ = 0; 
    for N=1:3
        NP1=N+1;
        BIG = abs( A(N,N) ); 
        if ( BIG < 1.0e-05)
            IBIG=N;
            for I=NP1:4
                if( abs(A(I,N)) <= BIG ) 
                    continue;
                end
                BIG = abs(A(I,N)); 
                IBIG = I;
            end
            if(BIG <= 0.)
                IERQ=2;
                return; 
            end
            if( IBIG  ~= N) 
                for J=N:4
                    TERM = A(N,J);
                    A(N,J) = A(IBIG,J);
                    A(IBIG,J) = TERM;
                end
                TERM = B(N); 
                B(N) = B(IBIG); 
                B(IBIG) = TERM;
            end
        end
        for I=NP1:4
            TERM = A(I,N)/A(N,N); 
            for J=NP1:4
                A(I,J) = A(I,J)-A(N,J)*TERM; 
            end
            B(I) = B(I)-B(N)*TERM; 
        end
    end
    if( abs(A(4,4)) > 0.0 )
        B(4) = B(4)/A(4,4);
        B(3) = (B(3)-A(3,4)*B(4))/A(3,3);
        B(2) = (B(2)-A(2,3)*B(3)-A(2,4)*B(4))/A(2,2);
        B(1) = (B(1)-A(1,2)*B(2)-A(1,3)*B(3)-A(1,4)*B(4))/A(1,1);
    else
        IERQ=1; % singular matrix 
            return;
    end
    end % gauss()
end % ecp()

