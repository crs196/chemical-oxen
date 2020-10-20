function [ ]=CrankAngleCalcs( )
% r = 9; % compression ratio
% s = 8.89; % stroke (cm)
% len= 13.335; %connecting rod length (cm)
% ep=s/(2*len);
% theta=-180:1:180; %crankangle theta vector
% ys1=(1-cosd(theta))/2; %approx y/s
% ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
% vol1 = 1+(r-1)*ys1; %approx volume
% vol2= 1+(r-1)*ys2; % exact volume
% %plot results
% figure(3)
% plot(theta,vol1,'--',theta,vol2,'-','linewidth',2);
% set(gca,'Xlim',[-180 180],'Ylim',[0 r],'fontsize',18,'linewidth',2);
% xlabel('Crank Angle (deg)','fontsize', 18);
% ylabel('Dim. Cylinder Volume','fontsize', 18);
% legend('Approx. Volume', 'Exact Volume','Location', 'North');

%% Pressure as a function of crank angle
% Gas cycle heat release code
% engine parameters
%clear();
thetas=-15; % start of heat release (deg)
thetad=40; % duration of heat release (deg)
r=9; %compression ratio
gamma= 1.31509; %gas const
b=.085725; %bore (m)
s=.0889; %stroke (m)
Vd=(pi/4)*(b^2)*s;
V1=Vd/(1-(1/r));
q= 4.8756/(148.249*V1); % dimensionless total heat release Qin/P1V1
a= 5; %weibe parameter a
n= 3; %weibe exponent n
step=1; % crankangle interval for calculation/plot
NN=360/step; % number of data points

% initialize the results data structure
save.theta=zeros(NN,1); % crankangle
save.vol=zeros(NN,1); % volume
save.press=zeros(NN,1); % pressure
save.work=zeros(NN,1); % work
pinit(1) = 1; %initial dimensionless pressure P/P1

j=1;
theta = -120; %initial crankangle
thetae = theta + step; %final crankangle in step
fy(1) = pinit(j); % assign initial pressure to working vector
fy(2) =0.; % reset work vector
% for loop for pressure and work calculation
for i=1:NN
    [fy, vol] = integrate(theta,thetae,fy);
    % reset to next interval
    theta = thetae;
    thetae = theta+step;
    % copy results to output vectors
    save.theta(i)=theta;
    save.vol(i)=vol;
    save.press(i,j)=fy(1);
    save.work(i,j)=fy(2);
end %end of pressure and work iteration loop
[pmax1, id_max1] = max(save.press(:,1)); %max pressure
thmax1=save.theta(id_max1);%crank angle
w1=save.work(NN,1);
eta1= w1/q; % thermal efficiency
imep1 = eta1*q*(r/(r -1)); %imep
eta_rat1 = eta1/(1-r^(1-gamma));
% output overall results
fprintf(' Inline-3 \n');
fprintf(' Theta_start %5.2f \n', thetas(1,1));
fprintf(' Theta_dur %5.2f \n', thetad(1,1));
fprintf(' P_max/P_1 %5.2f \n', pmax1);
fprintf(' Theta_max %7.1f \n',thmax1);
fprintf(' Net Work/P1V1 %7.2f \n', w1);
fprintf(' Efficiency %5.3f \n', eta1);
fprintf(' Eff. Ratio %5.3f \n', eta_rat1);
fprintf(' Imep/P1 %5.2f \n', imep1);

%plot results
figure(1)
plot(save.theta,save.press(:,1),'-','linewidth',1 )
set(gca, 'fontsize', 12,'linewidth',1);
legend('Inline-3', 'Location','NorthWest')
xlabel('Crank Angle (deg)','fontsize', 18)
ylabel('Pressure (bar)','fontsize', 18)
print -deps2 heatrelpressure

% figure(2);
% plot(save.theta,save.work(:,1),'-','linewidth',1)
% set(gca, 'fontsize', 12,'linewidth',1);
% legend('Inline-3','Location','NorthWest')
% xlabel('Theta (deg)','fontsize', 18)
% ylabel('Work','fontsize', 18)

    function[fy,vol] = integrate(theta,thetae,fy)
        %ode23 integration of the pressure differential equation
        %from theta to thetae with current values of fy as initial conditions
        [tt, yy] = ode23(@rates, [theta thetae], fy);
        for k=1:2
            fy(k) = yy(length(tt),k); %put last element of yy into fy vector
        end
        %pressure differential equation
        function [yprime] = rates(theta,fy)
            vol=(1.+ (r -1)/2.*(1-cosd(theta)))/r;
            dvol=(r - 1)/2.*sind(theta)/r*pi/180.; %dvol/dtheta
            dx=0.; %set heat release to zero
            if(theta > thetas) % then heat release dx > 0
                dum1=(theta -thetas)/thetad;
                x=1.- exp(-(a*dum1^n));
                dx=(1-x)*a*n*dum1^(n-1)/thetad; %dx/dthetha
            end
            term1= -gamma*fy(1)*dvol/vol;
            term2= (gamma-1)*q*dx/vol;
            yprime(1,1)= term1 + term2;
            yprime(2,1)= fy(1)*dvol;
        end %end of function rates
    end %end of function integrate2

%% Burn Fraction
theta=linspace(thetas,thetas+thetad,100); %crankangle theta vector
dum=(theta-thetas)/thetad; % theta diference vector
temp=-a*dum.^n;
xb=1.-exp(temp); %burn fraction
dxb=n*a*(1-xb).*dum.^(n-1); %element by element vector multiplication
%plot results
figure ()
plot(theta,xb,'b','linewidth',1);
set(gca, 'fontsize', 18,'linewidth',1);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Cumulative Burn Fraction','fontsize', 18);
figure();
plot(theta,dxb,'b','linewidth',1);
set(gca, 'fontsize', 18,'linewidth',1);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Burn Rate (1/deg)','fontsize', 18);

%% Cylinder Volume
len= 13.35; %connecting rod length (cm)
ep=(s*100)/(2*len); %epsilon
theta=-180:1:180; %crankangle theta vector
ys1=(1-cosd(theta))/2; %approx y/s
ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
vol2= 1+(r-1)*ys2; % exact volume
%plot results
plot(theta,vol2,'-','linewidth',1);
set(gca,'Xlim',[-180 180],'Ylim',[0 r],'fontsize',18,'linewidth',1);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Dim. Cylinder Volume (cc)','fontsize', 18);
legend('Exact Volume','Location', 'North');

%% Temperature
cv=0.718;
T=((save.press*100).*save.vol)/(cv*(gamma-1));

figure()
plot(save.theta,T)
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Temperature (K)','fontsize',18)

%% Instantaneous piston speed/mean piston speed
crankrad=4.445; %crank radius (cm) 4.445
len=120;
R=len/crankrad; %ratio of connecting rod length to the crank radius
sp_msp=(pi/2)*sin(theta)*(1+((cos(theta)/sqrt((R.^2)-(sin(theta).^2))))); %instantaneoud piston speed/mean piston speed as a function of crank angle as per Heywood p.45

figure()
plot(theta,sp_msp)
set(gca,'Xlim',[0 180],'fontsize',18,'linewidth',1);
xlabel('Crank Angle (deg)','fontsize', 18);

end % heat_release_weibe2