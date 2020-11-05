
%%coolant: 0.5'
clear
close all
clc
%%
Loading
oil_table %m2/s
vpiston %fps

% ALWYAS CHOOSE OPTION 2 BY DEFAULT
% option=input('option');
option=2;
theta=0:359; %degree
d_crank=[2.5 2.025 0.75]*0.0254; %m
% d_pin=2.025*0.0254;
d_cam=0.75*0.0254; %m

bore=3.375*0.0254; %m
stroke=3.5*0.0254; %m
engine_rpm=[7718,7718,7718/2]; %rpm
N=engine_rpm/60; %rps
lb_crank=0.35*bore; %m
r_oil=0.25*0.0254; %m
oil_channel_area=r_oil^2*pi; %m2
% sig_n= ;%fluid pressure as a function of theta
skirt=(0.8-0.65)*bore; %m
%% crankshaft & connecting rod & camshaft

% First loop iterates through rod selection
k=1;% for k=1:3 
    % Bearing diameter selection and annular clearing calculation
    db_crank=d_crank(k)+0.00001*0.0254; %bearing diameter
    rb_crank=db_crank/2;
    c_crank=db_crank-d_crank(k); %annular clearance of bearing, m
    
    % Angular velocity selection
    N_crank=N(k);
    
    % Second and third loop going through oil viscosity table
    for j=[1:4,6]
        for i=1:12
            mu=oiltable(i,j);
            
            % Selecting friction estimation method
            %%%% option 1: Petrov model
            %%%% option 2: Petrov model with eccentricity and bearing load
            %%%% Option 3: Coefficient of friction based on sommerfield 
                                                                %%%% number
            switch option
                case 1
                    ff.Ff_petrov{k}(i,j)=pi^2*db_crank^2*lb_crank*mu*N_crank/c_crank;
                case 2
                    % Initialization of epsilon(theta)
                    eccen=zeros(1,length(W));
                    % Solving for eccentricity as a function of theta
                    for p=1:length(W)
                        syms eccen_crank
                        eccentricity=vpasolve(-(W(p)/lb_crank)/(rb_crank*mu*N_crank)*...
                            (c_crank/rb_crank)^2*(db_crank/lb_crank)^2+...
                            pi*eccen_crank/(1-eccen_crank^2)^2*sqrt(0.62*eccen_crank^2+1)==0,eccen_crank);
                        
                        % Data correction
                        if isempty(eccentricity)||eccentricity<0
                            eccen(p)=0;
                        else
                            eccen(p)=eccentricity;
                        end
                    end
                    % Solving for attitude angle
                    phi=atan(pi*sqrt(1-eccen.^2)./(4*eccen));
                    
                    % Solving for minimum film thickness
                    ff.film{k}(i,j)=min(c_crank*(1-eccen));
                    e=c_crank*eccen;
                    
                    % Calculating friction force
                    ff.Ff_eccen{k}{i,j}=2*pi*mu*(rb_crank)^2*N_crank*lb_crank./(c_crank*sqrt(1-eccen.^2))+...
                        e*W/lb_crank*sin(phi);
                case 3
                    stribeck=mu*N_crank/sig_n;
                    f=stribeck*pi*db_crank(k)/c_crank;
            end
        end
    end
    oil_per_bearing(k)=pi*(db_crank^2-d_crank(k)^2)*lb_crank/4;
    %%
% end
oil_flow_rate=4*sum(oil_per_bearing([1,3]))*N(1)/(2*pi)...
    +3*oil_per_bearing(2)*N(1)/(2*pi);
oil_velocity=oil_flow_rate/oil_channel_area;
for i=1:12
    for j=[1:4,6]
        Re=rho*oil_velocity*l_channel/mu(i,j);
    end
end

%%
% phi=atan(pi*sqrt(1-eccen_crank^2)/(4*eccen_crank));
% ff.eccen=2*pi*mu*(rb_crank)^2*N_crank*lb_crank/(c_crank*sqrt(1-eccen_crank^2))+e*W/lb_crank*sin(phi);
%
% stribeck=mu*N_crank/sig_n;



%% Piston assembly friction
% stribeck=mu*abs(v_piston)/