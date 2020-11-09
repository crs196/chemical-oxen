function [ ]=BurnFraction( )
% this program computes and plots the cumulative burn fraction
% and the instantanous burnrate
clear();
close all
a = 5; % Weibe efficiency factor
n = 3; % Weibe form factor
thetas = -6; % start of combustion
thetad = 39.4; % duration of combustion
theta=linspace(thetas,thetas+thetad,100); %crankangle theta vector
dum=(theta-thetas)/thetad; % theta diference vector
temp=-a*dum.^n;
xb=1.-exp(temp); %burn fraction
dxb=n*a*(1-xb).*dum.^(n-1); %element by element vector multiplication
%plot results
plot(theta,xb,'b','linewidth',1);
set(gca, 'fontsize', 18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Cumulative Burn Fraction','fontsize', 18);
figure();
plot(theta,dxb,'b','linewidth',1);
set(gca, 'fontsize', 18,'linewidth',1);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Burn Rate (1/deg)','fontsize', 18);
end