%program four stroke ottofuel - computes const vol fuel air cycle
%establish initial conditions for intake stroke
clear;
Ti = 300; % intake temperature (K)
Pi = 101; % intake pressure (kPa)
Pe = 101; % exhaust pressure (kPa)
phi = .8; % equivalence ratio
rc = 9.; % compression ratio
%get specific heat input qin and stoichiometric fuel-air ratio FS
fuel_id=2 ; %id: 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
[ ~, ~, ~, ~, ~, ~, ~, ~, FS, qin ] = fuel( fuel_id, Ti );
% find enthalpy hi of intake fuel-air mixture
ff=0; % no residual fraction in intake air
[yi,hi,ui, si, vi, r, cpi, mw, dvdT, dvdP] = farg( Ti, Pi, phi, ff, fuel_id );
maxits = 100;
tol =0.0001;
%initial estimates of residual fraction and initial temp
f= .1; % residual fraction,
T1=350.; % initial temperature
P1=Pi; % no pressure drop in intake
% main iteration loop around cycle to get converged values of f and T1
for imain = 1:maxits
    % isentropic compression from v1 to known v2
    % call farg to get properties at 1
    [y1,h1,u1, s1, v1, r, cp1, mw, dvdT, dvdP] = farg( T1, P1, phi, f, fuel_id );
    %initial estimates of T2,P2
    v2=v1/rc;
    s2=s1;
    cv1=cp1+ T1*(dvdT^2)/dvdP;
    gam= cp1/cv1;
    T2=T1*(v1/v2)^(gam-1.);
    P2=P1*(v1/v2)^gam;
    %do the iteration to get T2 and P2 at end of compression
    for i2 = 1:maxits
        [y2, h2,u2, s2, v2, r, cp2, mw, dvdT, dvdP]=farg(T2,P2,phi,f,fuel_id);
        f1=s1-s2;
        f2=v1/rc - v2;
        det= cp2/T2*dvdP + dvdT^2;
        dt=(dvdP*f1 + dvdT*f2)/det;
        dp= (-dvdT*f1 + cp2/T2*f2)/det;
        %update T2 and P2
        T2=T2 + dt;
        P2=P2 + dp;
        %check for convergence
        %check for convergence
        if ( abs(dt)/T2 < tol && abs(dp)/P2 < tol )
            break;
        end
    end
    w12=-(u2-u1);% compression work
    % combustion from 2-3 with v and u constant
    %initial estimates of T3,P3 at state 3
    T3=3000;%Kelvin
    P3=7000; % kPa
    %do the iteration to get T3 and P3
    for i3 = 1:maxits
        [ierr,y3,h3,u3,s3,v3,r,cp3,mw,dvdT,dvdP] =ecp(T3,P3,phi,fuel_id);
        f1= u2-u3;
        f2= v2-v3;
        det= cp3*dvdP + T3*dvdT^2;
        dt= (-f1*dvdP - f2*(dvdT+dvdP))/det;
        dp= ((P3*dvdT-cp3)*f2 + f1*dvdT)/det;
        %update T3 and P3
        T3=T3 - dt;
        P3=P3 - dp;
        %check for convergence
        %check for convergence
        if ( abs(dt)/T3 < tol && abs(dp)/P3 < tol )
            break;
        end
    end
    % isentropic expansion of combustion products from v3 to known v4
    % initial estimates of T4,P4
    v4=v1;
    cv3=cp3+ P3*v3/T3*(dvdT^2)/dvdP;
    gam= cp3/cv3;
    T4=T3*(v3/v4)^(gam-1.);
    P4=P3*(v3/v4)^gam;
    %do the iteration to get T4 and P4
    for i4 = 1:maxits
        [ierr,y4,h4,u4,s4,v4,r,cp4,mw,dvdT,dvdP]=ecp(T4,P4,phi,fuel_id);
        f1=s3-s4;
        f2=rc*v3 - v4;
        det= cp4*dvdP/T4 + dvdT^2;
        dt=(dvdP*f1 + dvdT*f2)/det;
        dp= (-dvdT*f1 + cp4/T4*f2)/det;
        %update T4 and P4
        T4=T4 + dt;
        P4=P4 + dp;
        %check for convergence
        if ( abs(dt)/T4 < tol && abs(dp)/P4 < tol )
            break;
        end
    end
    % isentropic blowdown of control mass to exhaust pressure
    P5=Pe;
    s5=s4;
    %initial estimates of T5
    cv4=cp4+ P4*v4/T4*(dvdT^2)/dvdP;
    gam= cp4/cv4;
    T5=T4*(P5/P4)^((gam-1.)/gam);
    % do iteration for T5
    for i5 = 1:maxits
        [ierr,y5,h5,u5,s5,v5,r,cp5,mw,dvdT,dvdP]=ecp(T5,P5,phi,fuel_id);
        f1=s4-s5;
        dt=T5*f1/cp5;
        T5=T5 + dt;
        %check for convergence
        if ( abs(dt)/T5 < tol )
            break;
        end
    end
    % recompute residual fraction
    fold = f;
    v6=v5;
    f=v4/v6/rc;
    % constant pressure exhaust stroke
    h6=h5;
    %recompute h1 with fuel-air mixture and new residual fraction
    h1= f*(h6 + (Pi-Pe)*v6) + (1-f)*hi;
    T1old=T1;
    %recompute T1 with latest f and h1
    for i6 = 1:maxits
        [y1,h1new,u1,s1,v1,r,cp1,mw,dvdT,dvdP]=farg(T1,Pi,phi,f,fuel_id);
        g=h1new-h1;
        dt=-g/cp1;
        T1=T1+dt;
        %check for convergence
        if ( abs(dt)/T1 < tol )
            break;
        end
    end
    %check for convergence of main iteration loop
    dt=T1old - T1;
    df=fold - f;
    if ( abs(dt)/T1 < tol && abs(df)/f < tol )
        break;
    end
end %end of main iteration loop
%compute cycle parameters
w = u1-u4;% net work
imep = w/(v1-v2); %imep
eta = w*(1+phi*FS)/phi/FS/(1.-f)/qin; % thermal efficiency
pmep = Pe -Pi; %pmep
etanet = eta*(1.- pmep/imep); %net thermal efficiency
ev = rc*(1.-f)*vi/(rc-1.)/v1;
%output state and cycle parameters
fprintf(' \n Ottofuel inlet: Temp (K)= %5.1f Pressure (kPa)= %5.1f phi= %6.2f fuel= %3d \n', Ti, Pi, phi, fuel_id );
fprintf(' State \t\t 1 \t 2 \t \t 3\t \t 4 \n')
fprintf('Pressure (kPa)= %7.1f \t %7.1f \t %7.1f \t %7.1f \n',P1,P2,P3,P4);
fprintf('Temperature (K)=%7.1f \t %7.1f \t %7.1f \t %7.1f \n',T1,T2,T3,T4);
fprintf('Enthalpy(kJ/kgK)=%7.1f \t %7.1f \t %7.1f \t %7.1f \n',h1,h2,h3,h4);
fprintf('Int. Energy(kJ/kg)=%6.1f \t %7.1f \t %7.1f \t %7.1f \n',u1,u2,u3,u4);
fprintf('Volume (m^3/kg)= %7.3f \t %7.3f \t %7.3f \t %7.3f \n', ....
    v1,v2,v3,v4);
fprintf('Entropy(kJ/kgK)= %6.3f \t %7.3f \t %7.3f \t %7.3f \n', ...
    s1,s2,s3,s4);
fprintf('Cp (kJ/kg K) = %7.3f \t %7.3f \t %7.3f \t %7.3f \n \n', ...
    cp1,cp2,cp3,cp4);
fprintf('Work (kJ/kg)= %7.1f \t \t Volumetric Efficiency= %7.4f \n', w, ev);
fprintf('Ideal Thermal Efficiency= %7.3f \t Net Thermal Efficiency= %7.4f \n', eta, etanet);
fprintf('Imep (kPa)= %7.1f \t \t \t Pmep (kPa)= %7.1f \n', imep, pmep);
fprintf('Exhaust Temperature (K)= %7.1f \t Residual Mass Fraction f =%7.4f \n', T5,f);
fprintf('Iterations = \t %4d\t %4d \t %4d \t %4d \n', imain,i2,i3,i4);