
%%coolant: 0.5'
clear
close all
clc

oil_table
vpiston %fps

% ALWYAS CHOOSE OPTION 2 BY DEFAULT
option=input('option');
theta=0:359; %degree
d_crank=[2.5 2.025 0.75]*0.0254;
% d_pin=2.025*0.0254;
d_cam=*0.0254;

bore=3.375*0.0254; %in
stroke=3.5*0.0254; %in
engine_rpm=[7718,7718,7718/2];

% sig_n= ;%fluid pressure as a function of theta
skirt=(0.8-0.65)*bore;
r_oil=0.25*0.0254;
% db_rod
%% crankshaft & connecting rod & camshaft

% W bearing load
% N angular velocity
N=engine_rpm/60;
lb_crank=0.35*bore;
for k=1:3
    db_crank=d_crank(k)+0.00001*0.0254; %bearing diameter
    rb_crank=db_crank/2;
    c_crank=db_crank-d_crank(k); %annular clearance of bearing, in

    N_crank=N(k);
    for j=[1:4,6]
        for i=1:12
            mu=oiltable(i,j);
            switch option
                case 1
                    ff.Ff_petrov{k}(i,j)=pi^2*db_crank^2*lb_crank*mu*N_crank/c_crank;
                case 2
                    eccen=zeros(1,length(W));
                    for p=1:length(W)
                        syms eccen_crank
                        eccen(p)=solve((W(p)/lb_crank)/(rb_crank*mu*N_crank)*(c_crank/rb_crank)^2*(db_crank/lb_crank)^2==...
                            pi*eccen_crank/(1-eccen_crank^2)^2*sqrt(0.62*eccen_crank^2+1),eccen_crank);
                    end
                    phi=atan(pi*sqrt(1-eccen.^2)./(4*eccen));
                    ff.film{k}(i,j)=min(c_crank*(1-eccen));
                    ff.Ff_eccen{k}(i,j)=2*pi*mu*(rb_crank)^2*N_crank*lb_crank/(c_crank*sqrt(1-eccen^2))+e*W/lb_crank*sin(phi);
                case 3
                    stribeck=mu*N_crank/sig_n;
                    f=stribeck*pi*db_crank(k)/c_crank;
            end
        end
    end
    oil_per_bearing(k)=pi*(db_crank^2-d_crank(k)^2)*lb_crank/4;
end

oil_flow_rate=4*sum(oil_per_bearing([1,3]))*N(1)/(2*pi)...
    +3*oil_per_bearing(2)*N(1)/(2*pi);

oil_channel_area=r_oil^2*pi;
oil_velocity=oil_vol/oil_channel_area;

Re=rho*oil_velocity*l_channel/mu;
%%
% phi=atan(pi*sqrt(1-eccen_crank^2)/(4*eccen_crank));
% ff.eccen=2*pi*mu*(rb_crank)^2*N_crank*lb_crank/(c_crank*sqrt(1-eccen_crank^2))+e*W/lb_crank*sin(phi);
% 
% stribeck=mu*N_crank/sig_n;



%% Piston assembly friction
% stribeck=mu*abs(v_piston)/