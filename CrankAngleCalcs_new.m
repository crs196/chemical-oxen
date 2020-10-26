function [T]=CrankAngleCalcs_new(rc,re,thetas,thetad,a,n,b,stroke,len,Cp,MW,dMdT,h,R,dMdP,v,P)
%% Cylinder Volume v. Crank Angle
%re = 10; % expansion ratio
%st = 8.89; % stroke (cm)
%b = 8.5725; % bore (cm)
%len= 13.335; %connecting rod length (cm)
ep=(stroke*100)/(2*len); %episolon (change units of stroke from m to cm)
theta=-180:1:180; %crankangle theta vector
ys1=(1-cosd(theta))/2; %approx y/s
ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
vol1 = 1+(rc-1)*ys1; %approx volume dimm
vol2= 1+(rc-1)*ys2; % exact volume dimm
Vol1 = ys1*stroke*pi*((b*100)/2)^2; %approx volume (change units of bore from m to cm)
Vol2 = ys2*stroke*pi*((b*100)/2)^2; % exact volume (change units of bore from m to cm)
%plot results
figure()
plot(theta,Vol2)
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Cylinder Volume (cc)','fontsize', 18)
%Vol2 = Vol2';
%% Pressure as a function of crank angle

thetas=-15; % start of heat release (deg)
thetad=40; % duration of heat release (deg)
gamma= 1.31509; %gas const
Vd=(pi/4)*(b^2)*stroke;
V1=Vd/(1-(1/rc));
q= 4.8756/(148.249*V1); % dimensionless total heat release Qin/P1V1
step=1; % crankangle interval for calculation/plot
NN=360/step; % number of data points

% initialize the results data structure
save.theta=zeros(NN,1); % crankangle
save.vol=zeros(NN,1); % volume
save.press=zeros(NN,1); % pressure
save.work=zeros(NN,1); % work
pinit(1) = 1; %initial dimensionless pressure P/P1

j=1;
theta = -180; %initial crankangle
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
imep1 = eta1*q*(rc/(rc -1)); %imep
eta_rat1 = eta1/(1-rc^(1-gamma));
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
figure()
plot(save.theta,save.press(:,1),'-','linewidth',1 )
set(gca, 'fontsize', 12,'linewidth',1);
legend('Inline-3', 'Location','NorthWest')
xlabel('Crank Angle (deg)','fontsize', 18)
ylabel('Pressure (bar)','fontsize', 18)
print -deps2 heatrelpressure
% figure();
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
            vol=(1.+ (rc -1)/2.*(1-cosd(theta)))/rc;
            dvol=(rc - 1)/2.*sind(theta)/rc*pi/180.; %dvol/dtheta
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
ep=(stroke*100)/(2*len); %epsilon
theta=-180:1:180; %crankangle theta vector
ys1=(1-cosd(theta))/2; %approx y/s
ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
vol2= 1+(rc-1)*ys2; % exact volume
%plot results
plot(theta,vol2,'-','linewidth',1);
set(gca,'Xlim',[-180 180],'Ylim',[0 rc],'fontsize',18,'linewidth',1);
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
%% Gamma - Work in Progress
%gamma=cp/cv=1+(r/cv);
%R=0.287; %kJ/kg K (uncomment line 156 when the code actually works)
R = 8.31434/MW;
Cp = R*(Cp - h*(((save.press*100).*Vol2*10^(-6))/(cv*(gamma(1,:)-1)))*dMdT/MW);
gamma=(Cp/(Cp+(T*(((-v/T)-(v/MW)*dMdT).^2))/((-v/(save.press*100))-(v/MW)*dMdP)))+0.3827;
Gamma=gamma';

figure()
plot(theta(1:end-1),Gamma(1,:))
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Gamma','fontsize',18)

%% Temp part 2
gam=gamma(2,:);
T2=((save.press*100).*save.vol)./(cv*(gam-1));
figure()
plot(save.theta,T2)
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Temperature (K)','fontsize',18)

%% Instantaneous piston speed/mean piston speed
% crankrad=4.445; %crank radius (cm) 4.445
% len=13.335;
% R=len/crankrad; %ratio of connecting rod length to the crank radius
% sp_msp=(pi/2)*sin(theta)*(1+((cos(theta)/sqrt((R.^2)-(sin(theta).^2))))); %instantaneoud piston speed/mean piston speed as a function of crank angle as per Heywood p.45
% 
% figure()
% plot(theta,sp_msp)
% set(gca,'Xlim',[0 180],'fontsize',18,'linewidth',1);
% xlabel('Crank Angle (deg)','fontsize', 18);
end % heat_release_weibe2