% Follower Friction and Lubrication

% flat follower
cff = 133; %flat follower coefficient, kPa-mm
nv = 12; %total number of valves
s = 88.9; %piston stroke, mm
b = 85.725; %cylinder bore, mm
N = 7714.3; %engine speed/relative rotational speed between surfaces
nc = 3; %number of cylinders
T = 2920.290; %max temperature of cycles in K
T = T - 273;
C1 = [1.09e-4, 9.38-5, 9.73e-5, 8.35e-5, 1.17e-4, 1.29e-4];
C2 = [1157.5, 1271.6, 1360.0, 1474.4, 1509.6, 1564.0];
mu = C1.*exp(C2./(1.8.*T + 127)); %lubricant dynamic viscosity, N-s/m^2
fmepff = cff.*(2+10./(5+mu.*N)).*(nv./(nc.*s)); %mean effective pressure
fmep_ff = min(fmepff);
index = find(fmep_ff);
mu_desired = mu(index);
fprintf('The mean effective pressure of the flat follower: %2.2f kPa\n', fmep_ff)
fprintf('The desired lubrucating oil is SAE 20 with a dynamic viscocity of %f N-s/m2\n', mu_desired)

% % roller follower
% crf = 0.0050; %roller follower coefficient, kPa-mm-min/rev
% fmep_rf = crf.*(nv.*N)./(nc.*s);
% fprintf('The mean effective pressure of roller follower: %f\n', fmep_rf)
