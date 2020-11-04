function [LD,CD_e]=cd_chart_exhaust()
% tl=[152,57]
% bl=[152,503]
% br=[873,503]
% tr=[873,57]
clean

tl=[0,503-57];
bl=[0,0];
tr=[873-152,503-57];
br=[873-152,0];
bound=[tl;tr;br;bl];
plot(bound(:,1)/(873-152)*0.4,bound(:,2)/(503-57))
cd_location=[193,182;239,188;286,183;333,185;372,178;419,162;464,161;559,175;650,199;739,235;829,266];
cdy=(1-(cd_location(:,2)-57)/(503-57));
cdx=(cd_location(:,1)-152)/(873-152)*0.4;
hold on
plot(cdx,cdy,'.')

dx=0.01*10^-3;
LD=round(cdx(1):dx:cdx(end),5);
CD_e=interp1(cdx,cdy,LD,'pchip');
plot(LD,CD_e)
end