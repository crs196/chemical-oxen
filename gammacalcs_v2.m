%% lets try to run both runecp.m and runfarg.m from the same script!
% this is the run file for ecp, farg, fuel, and crankanglecalcs_new

phi = 0.765; % enter equivalence ratio input
[Temp,pressure]=FiniteHeatRelease();
%Temp = [463.686,652.816,2920.290,1413.622,209.891,300,300,298.797,463.686]; % enter temperature (K) input
%Press = [141.870,3141.366,14052.509,680.238,101,101,160,160,141.870]; % enter pressure (kPa) input
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
gamv=zeros(9,1);

for i=1:length(Temp)
    T=Temp(i);
    P=pressure(i);
    
    fuel_id = 2;
    if T<1000
        %Input-Output program for running farg.m clear;
        f=0.1111; %residual fraction input
        % fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
        % call farg function
        [Y,h,u,s,v,R,Cp,MW,dvdT,dvdP,gam_test] = farg(T,P,phi,f,fuel_id);
        
    elseif (T>600 || T<3500) || (P>20 || P<30000)
        
        % fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
        % call ecp function
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP,dMdT,dMdP,gam_test] = ecp( T, P, phi, fuel_id );
        
    else
        fprintf('Welp try a diff temp or pressure value \n')
    end
    gamv(i)=gam_test;
    
end

figure()
plot(-180:180,gamv)
xlabel('Theta(deg)')
ylabel('Gamma')
[Temp,Gamma]=CrankAngleCalcs_nam(rc,re,thetas,thetad,a,n,b,stroke,len,Cp,MW,h,R,v,P,dvdT,dvdP);