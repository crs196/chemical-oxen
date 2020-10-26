function [ ] = WoschniHeatTransfer( )
% Gas cycle heat release code with Woschni Heat Transfer
clear();
thetas = -15; %start of heat release (deg)
thetad = 40; %duration of heat release (deg)
r = 9; %compression ratio
R = 0.287; 
gamma =1.31509; %gas constant
beta = 1.5; %dimensionless volume
a = 5; %weibe parameter a
n = 3; %weibe exponent n
omega = 200.; %engine speed rad/s
c = 0; %mass loss coeff
s = 0.0889; %stroke (m)
b = 0.085725; %bore (m)
Vd=(pi/4)*(b^2)*s;
V1=Vd/(1-(1/r));
Q= 4.8756/(148.249*V1); % dimensionless total heat release Qin/P1V1
% from initial code Q = 20; %dimensionless total heat release
T_bdc = 300; %temp at bdc (K)
tw = 1.2; %dimensionless cylinder wall temp
P_bdc = 100; %pressure at bdc (kPa)
Up = s*omega/pi; %mean piston speed (m/s)
%dens = R/(T_bdc*P_bdc);
%m_bdc = dens*Up*(V1/s);
step=1; %crankangle interval for calculation/plot
NN = 360/step; %number of data points
theta = -120; %initial crankangle
thetae = theta + step; %final crankangle in step
%initalize results data structure
save.theta=zeros(NN,1);
save.vol=zeros(NN,1); %volume
save.press=zeros(NN,1); %pressure
save.work=zeros(NN,1); %work
save.heatloss=zeros(NN,1); %heat loss
save.mass=zeros(NN,1); %mass left
save.htcoeff=zeros(NN,1); %heat transfer coeff
save.heatflux=zeros(NN,1); %heat flux (W/m^2)
save.massflow=zeros(NN,1); %massflow (kg-m/s) - added
fy = zeros(4,1); %vector for calculated pressure, work, heat and mass loss
fy(1) = 1; %initial pressure (P/P_bdc)
fy(4) = 1; %initial mass (-)
%for loop for pressure and work calculations
for i = 1:NN
    [fy, vol, ht,hflux] = integrate_ht(theta,thetae,fy);
    %print values
    %fprintf('%7.1f %7.2f %7.2f %7.2\n', theta,vol,fy(1),fy(2),fy(3));
    
    %reset to next interval 
    theta=thetae;
    thetae=theta+step;
    save.theta(i)=theta; %put results in output vectors
    save.vol(i)=vol;
    save.press(i)=fy(1);
    save.work(i)=fy(2);
    save.heatloss(i)=fy(3);
    save.mass(i)=fy(4);
    save.htcoeff(i)=ht;
    save.hflux(i)=hflux;
    %save.massflow(i)=m; %added
end %end of pressure and work of loop
[pmax, id_max] = max(save.press(:,1)); %find max pressure
thmax=save.theta(id_max); %and crank angle
ptdc=save.press(NN/2)/pmax; %pressure top-dead-center
w=save.work(NN,1); %w is cumulative work
massloss=1- save.mass(NN,1);
eta=w/Q; %thermal efficiency
imep=eta*Q*(r(r-1)); %imep/P1V1
eta=w/Q; %thermal efficiency
eta_rat=eta/(1-r^(1-gamma));
%output overall results
fprintf('Weibe Heat Release with Heat and Mass Loss \n');
fprintf('Theta_start = %5.2f \n', thetas);
fprintf('Theta_dur = %5.2f\n', thetad);
fprintf('P_max/P1 = %5.2f\n', pmax);
fprintf('Theta @ P_max = %7.1f\n',thmax);
fprintf('Net Work/P1V1 = %7.2f]n',w);
fprintf('Heat Loss/P1V1 = %7.2f\n',save.heatloss(NN,1));
fprintf('Mass Loss/P1V1 = %7.2f\n', massloss);
fprintf('Efficiency= %5.2f\n',eta);
fprintf('Eff./Eff. Otto = %5.2f\n',eta_rat);
fprintf('Imep/P1= %5.2f\n', imep);
%plot results
figure(1);
plot(save.theta,save.work,'-',save.theta,save.heatloss,'--','linewidth',2);
set(gc,'Xlim',[-180 180], 'fontsize',18,'inewidth',1.5);
hleg1=legend('Work','Heat Loss','Location','NorthWest');
set(hlegl,'Box','off')
xlabel('Crank Angle \ theta (deg)','fontsize', 18);
ylabel('Cumulative Work and Heat Loss','fontsize',18);
plot(save.theta,save.press,'-','linewidth',2);
set(gca,'fontsize',19,'inewidth',1.5,'Xlim', [-180 180]);
xlabel('Crank Angle (deg)','fontsize',18);
ylabel('Pressure (bar)','fontsize',18);
figure(2);
plot(save.theta,save.htcoeff,'-','linewidth',2);
set(gca,'fontsize',18,'linewidth',1.5,'Xlim', [-180 180]);
xlabel('Crank Angle \theta (deg)','fontsize', 18);
ylabel('Heat Transfer Coefficient (h) (W/m^s-K)','fontsize',18);
figure(3);
plot(save.theta,save.hflux,'-','linewidth',2);
set(gca,'fontsize',18,'linewidth',1.5,'Xlim',[-180 180]);
xlabel('Crank Angle \theta (deg)','fontsize',18);
ylabel('Heat flux q{"} (MW/m^2)','fontsize',18);

function[fy,vol,ht,hflux] = integrate_ht(theta,thetae,fy)
%ode23 integration of the pressure differential equation
%from theta to thetae with current values of fy as initial conditions
[tt,yy] = ode23(@rates, [theta thetae], fy);
%put last element of yy into fy vector
for j=1:4
    fy(j) = yy(length(tt,j));
end
%pressure differential equation 
function [yprime] = rates(theta,fy)
vol=(1.+ (r-1)/2.*(1-cosd(theta)))/r;
dvol=(r-1)/2.*sind(theta)/r*pi/180.; %dvol/dtheta
dx=0;
if (theta>thetas) %heat release > 0 
    dum1= (theta -thetas)/thetad;
    x=1-exp(-(a*dum1^n));
    dx=(1-x)*a*n*dum1^(n-1)/thetad; %dx/dtheta
end
P=P_bdc*fy(1); %P(theta) (kPa)
T=T_bdc*fy(1)*vol; %T(theta) (K)
term4=T_bdc*(r-1)*(fy(1)-vol^(-gamma))/r; %comb. vel. increase
U=2.28*Up + 0.00324*term4; %Woschni vel (m/s)
ht = 3.26*(P^(0.8))*U*(0.8)*b^(-0.2)*T*(-0.55); %Woschni ht coeff
hflux=ht*T_bdc*(fy(1)*vol/fy(4) - tw)/10^6; %heat flux MW/m^2
h = ht*T_bdc*4/(1000*P_bdc*omega*beta*b); %dimensionless ht coeff
term1 = -gamma*fy(1)*dvol/vol;
term3 = h*(1. +beta*vol)*(fy(1)*vol/fy(4) - tw)*pi/180.;
term2 = (gamma-1)/vol*(Q*dx - term3);
yprime(1,1)=term1 + term2 - gamma*c/omega*fy(1)*pi/180;
yprime(2,1)=fy(1)*dvol;
yprime(3,1)=term3;
yprime(4,1)=-c*fy(4)/omega*pi/180;
%mass flow 
d = R/(P*T);
Area = (save.vol)/s;
m = d*U*Area; %mass flow after combusion for engine, will need compression scale for after turbo
m_turbo = m*7.55986;

end %end of function rates

end %end of function integrate_ht

end %end of function HeatReleaseHeatTransfer


