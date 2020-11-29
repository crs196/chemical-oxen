% Run Homogenous.m
% Outputs plots of eta vs rpm, T vs rpm, NOx vs rpm, grams of Nox per... 
%   km vs rpm at 60mph

% takes a while to run btw
RPM = 1250;
nox = zeros(25,1);
eta = zeros(25,1);
index = 1;
while RPM <= 7500
    RPM = RPM + 250;
    [ ETA, ~, NOX_ppm,~,~] = Homogeneous_range(RPM);
    nox(index) = NOX_ppm;
    eta(index) = ETA;
    index = index + 1;
end

nox_begin = zeros(3,1);
eta_begin = zeros(3,1);

for i = 1:3
    [ ETA, ~, NOX_ppm,~,~] = Homogeneous_range(800);
    nox_begin(i) = NOX_ppm;
    eta_begin(i) = ETA;
    [ ETA, ~, NOX_ppm,~,~] = Homogeneous_range(1000);
    nox_begin(i) = NOX_ppm;
    eta_begin(i) = ETA;
    [ ETA, ~, NOX_ppm,~,~] = Homogeneous_range(1124);
    nox_begin(i) = NOX_ppm;
    eta_begin(i) = ETA;
end

nox_end = zeros(1,1);
eta_end = zeros(1,1);

for j = 1
    [ ETA, ~, NOX_ppm,~,~] = Homogeneous_range(7714);
    nox_end(j) = NOX_ppm;
    eta_end(j) = ETA;
end

nox = [nox_begin ; nox ; nox_end];
eta = [eta_begin ; eta ; eta_end];

rpm = [1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,...
    4250,4500,4750,5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500]';

rpm = [800;1000;1124;rpm;7714];

figure()
plot(rpm,nox)
ylabel('Concentration of NOx (ppm)')
xlabel('RPM')

figure()
plot(rpm,eta)
ylabel('Thermal Efficiency')
xlabel('RPM')
    
efficiency = [21.08,16.86,15.00,13.49,11.24,9.63,8.43,7.49,6.74,6.13,...
    5.62,5.19,4.82,4.5,4.22,3.97,3.75,3.55,3.37,3.21,3.07,2.93,2.81,...
    2.7,2.59,2.5,2.41,2.33,2.25,2.19]'; % (km/L) for 60MPH

emission_nox = (nox./1000).*efficiency;

figure()
plot(rpm,emission_nox)
ylabel('Grams of NOx per km (60MPH)')
xlabel('RPM')

figure()
plot(rpm,emission_nox.*0.05)
ylabel('Grams of NOx per km (60MPH) (with Catalytic Converter)')
xlabel('RPM')
