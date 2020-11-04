% heywood cd graph
function [LD,CD]=cd_chart_intake()
% clean
ay=[15.3/16*.2+.4,.6,19/16*.2+.4,20.3/16*.2+.4];
ax=[8.5,9.5,13.4,14.3]*.1/15.8;

by=[13,15,19.3]/16*.2+.4;
bx=[1.5,5,9.5]*0.1/15.5+.1;

cy=[18.2,12,8.8,6.8]/16*0.2+0.4;
cx=[-3 1.8 6.2 9.4]/15.5*0.1+0.2;

cdx=[ax,bx,cx];
cdy=[ay,by,cy]; 

plot(cdx,cdy,'.','markersize',20)
xlim([0,.3])
ylim([0.4,.8])
%%
i=3;
maxnondimlift=0.5;
dx=0.01*10^-i;
aintpl=ax(1):dx:ax(end);
a=interp1(ax,ay,aintpl,'pchip');
% b=interp1(bx,by,linspace(bx(1),bx(end)),'pchip','extrap');
bintpl=bx(1):dx:0.2;
b=interp1(bx,by,bintpl,'pchip','extrap');
% c=interp1(cx,cy,linspace(cx(1),cx(end)),'pchip');
cintpl=0.15:dx:maxnondimlift;
c=interp1(cx,cy,cintpl,'pchip','extrap');

axbx=ax(end):dx:bx(1);
% ab=interp1([ax(end),bx(1)],[ay(end),by(1)],axbx,'linear');
ab=round(interp1([ax(end),bx(1)],[ay(end),by(1)],axbx,'linear'),4);
%%

figure(1)
hold on
plot(aintpl,a)
plot(axbx,ab)

% plot(linspace(bx(1),bx(end)),b)
% plot(linspace(cx(1),cx(end)),c)
plot(bintpl(1:6026),b(1:6026))
plot(cintpl(1998:end),c(1998:end))
%%
LD=round([aintpl axbx bintpl(1:6026) cintpl(1998:end)],5);
CD=[a,ab,b(1:6026),c(1998:end)];
figure(2)
plot(LD,CD)
end
