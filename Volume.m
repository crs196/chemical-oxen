function [ ]=Volume( )
clear();
r = 9; % compression ratio
s = 8.901796; % stroke (cm)
len= 120; %connecting rod length (cm)
ep=s/(2*len);
theta=-180:1:180; %crankangle theta vector
ys1=(1-cosd(theta))/2; %approx y/s
ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
vol1 = 1+(r-1)*ys1; %approx volume
vol2= 1+(r-1)*ys2; % exact volume
%plot results
plot(theta,vol1,'--',theta,vol2,'-','linewidth',2);
set(gca,'Xlim',[-180 180],'Ylim',[0 r],'fontsize',18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Dim. Cylinder Volume','fontsize', 18);
legend('Approx. Volume', 'Exact Volume','Location', 'North');
end
