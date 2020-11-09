function [ ]=FunctionofCrankAngle()
%% This is a run file for FiniteHeatRelease.m, Volume.m, and others
% The graphs we are looking for here are:
% Temperature v. Crank angle
% Pressure v. Crank angle
% Volume v. Crank angle
% Gamma v. Crank angle

% We are then looking to find these graphs with a Gamma that varies with
% temperature

% First lets clear our Workspace
close all
clear

% In order to calculate gamma we need:
% Cp = specific heat capacity [kJ/kg K]
% T = temperature [K]
% dvdT = derivative of volume with respect to temperature [m^3/kg K]
% dvdP = derivative of volume with respect to pressure

% It is cleanest to calculate gamma using ecp.m or farg.m depending on
% the pressure and temperature at that point in the cycle

% The inputs for ecp.m and farg.m are:

% ecp inputs:
% T - temperature (K) [ 600 --> 3500 ]
% P - pressure (kPa) [ 20 --> 30000 ]
% phi - equivalence ratio [ 0.01 --> ~3 ]
% fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane

% farg inputs:
% T - temperature (K) [ 300 --> 1000 K ]
% P - pressure (kPa)
% phi - equiv ratio = actual fuel-to-air ratio/stoich fuel-air ratio
% (stoich fuel-air ratio: 14.7)
% f - residual fraction
% fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane

% From these inputs the Temperature and the Pressure will change as a
% function of crank angle. The other values will not.

% Let's start by getting Temperature and Pressure from FiniteHeatRelease.m

[Temp,press,theta_test]=FiniteHeatRelease();

% The pressure value is correct however the temperature values are
% experimental

% Now let's set Temp and pressure to be T and P so that we can use them in
% later functions.

Temperature=Temp(:,1); % Temperature was is now 1 column and is in kelvin
Pressure=press(:,1); % Pressure is now 1 column and is in kPa

% Now we'll initialize the rest of our values and go into a for loop for
% to get Cp values for all temperatures
gamv=zeros(length(Temperature),1);
cp_gamma=zeros(length(Temperature),1);
dvdT_gamma=zeros(length(Temperature),1);
dvdP_gamma=zeros(length(Temperature),1);

for i=1:length(Temperature)
    
    % Choosing the temperature and pressure values for the iteration and
    % setting the fuel_id and equivalence ratio for our fuel/engine
    
    T=Temperature(i);
    P=Pressure(i);
    fuel_id = 2;
    phi = 0.765;
    
    if T<1000
        % disp('running farg.m')
        f=0.1111; % residual fraction input
        % call farg function
        [Y,h,u,s,v,R,Cp,MW,dvdT,dvdP] = farg(T,P,phi,f,fuel_id);
        
    elseif T<3500 && P<30000
        % disp('running ecp.m')
        % call ecp function
        [ierr,Y,h,u,s,v,R,Cp,MW,dvdT,dvdP] = ecp(T,P,phi,fuel_id );
        
    else
        disp('something went wrong')
    end
    
    % Now that we have calculated these values we can calculate gamma.
    % Gamma will be calculated in the variable gam_test
    
    gam_test=Cp/(Cp+T*(dvdT^2)/dvdP);
    
    % Now that we have calculated gamma we can store gamma, Cp, dvdT, and
    % dvdP in the arrays that we initialized before going into the for loop
    
    gamv(i)=gam_test;
    cp_gamma(i)=Cp;
    dvdT_gamma(i)=dvdT;
    dvdP_gamma(i)=dvdP;
end

%% Plot Gamma v. Crank Angle
%theta_plot=linspace(-180,180,length(Temperature));
figure()
plot(theta_test,gamv)
xlabel('Crank Angle (deg)')
ylabel('Gamma')

%
% %% Initializing Variables
%     % r=9; %compression ratio
%     % re=10; %expansio ratio
%     b=.085725; %bore (m)
%     stroke=.0889; %stroke (m)
%     Vd=(pi/4)*(b^2)*stroke;
%     % V1=Vd/(1-(1/r));
%
% %% Volume
%     [dim_vol]=Volume( );
%     v_theta=dim_vol*Vd;
%
%     figure()
%     plot(-180:180,v_theta)
%     xlabel('Crank Angle (deg)')
%     ylabel('Cylinder Volume (m^3)')
%
% %% Temperature
% %     N=74.25; %Number of moles in octane. from ecp and farg. this does not change
% %     v_spec=1.423*10^-3; %ideal gas specific volume of octane
% %     R=8.314; %universal gas constant
% %     T=(save.press*100).*dim_vol(1:end-1)'/(N*R);
% %
% %     figure()
% %     plot(save.theta,T)
% %     xlabel('Crank Angle (deg)')
% %     ylabel('Temperature (K)')
% %
end % FunctionofCrankAngle.m ends here