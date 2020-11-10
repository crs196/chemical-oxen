function [ ETA, IMEP, NOX_ppm,TU,TB] = Homogeneous(varargin)
close all
% Two Zone Arbitrary Heat Release (Fuel Inducted Engine)
%
% Find:
% 1. Indicated thermal efficiency - ETA
% 2. Indicated mean effective pressure - IMEP (kPa)
%
% Note:
% 1. Cosine burning law is employed
% 2. Fuel is gasoline, C7H17
R = 9; % Compression ratio - R
B = .085725; % Bore - B (m)
S = .0889; % Stroke - S (m)
EPS = .5*S/.1335; % Half stroke to rod ratio - EPS
RPM = 7714; % Engine speed - RPM
HEAT = 500; % Heat transfer coefficient
BLOWBY = 0.8; % Blowby coefficient
THETAS = -6; % Start of heat release (deg ATDC)
THETAB = 39; % Burn angle (deg)
PHI = 0.765; % Equivalence ratio - PHI
F = 0.1111; % Residual fraction - F
TW = 420; % Wall temperature - TW
fuel_type = 2; % gasoline
FS = 0.06548; % stoichiometric fuel-air ratio for gasoline
A0 = 47870; %
T1 = 463.686; %K
P1 = 141.870; % kPa
if ( nargin == 3 )
    PHI = varargin{1};
    F = varargin{2};
    RPM = varargin{3};
end
OMEGA = RPM*pi/30;
to_ppm = 10^6; % convert from mass fraction to ppm
MW_NO = 30; % molecular weight of NO, g/mol
THETA = -180;
DTHETA = 1;
THETAE = THETA+DTHETA;
[ VOL, X, EM ] = auxiliary( THETA );
NNOX = THETAB/DTHETA;
NY = 6+NNOX;
Y = zeros(NY,1);
Y(1) = P1;
Y(2) = nan;
Y(3) = T1;
[~, ~, ~, ~, vU, ~, ~, ~, ~, ~] = farg( Y(3), Y(1), PHI, F, fuel_type );
MNOT = VOL/vU;
M = EM*MNOT;
NN = 36*10;
SAVE.THETA = zeros( NN, 1 );
SAVE.VOL = zeros( NN, 1 );
SAVE.T = zeros(NN, 1 );
SAVE.P = zeros( NN, 1 );
SAVE.MDOTFI = zeros( NN, 1 );
SAVE.NOx = zeros(NN,5);
fprintf( 'THETA VOL BURN FRAC PRESS BURN TEMP UNBURNED T WORK HEAT LOSS MASS H-LEAK NOx\n' );
fprintf( ' deg cm^3 -- kPa K K J J g J ppm\n' );
fprintf('%7.1f %6.1f %3.3f %6.1f %6.1f %6.1f %5.0f %5.0f %5.3f %5.2f %6.1f\n', THETA, VOL*1000000, X, Y(1), Y(2), Y(3), Y(4)*1000, ...
    Y(5)*1000, M*1000, Y(6)*1000, 0.0 );
II = 1;
for III=1:36
    for JJJ=1:10
        [ Y ] = integrate( THETA, THETAE, Y );
        [ VOL, X, EM ] = auxiliary( THETA );
        M = EM*MNOT;
        THETA=THETAE;
        THETAE=THETA+DTHETA;
        % save data for plotting later
        SAVE.THETA(II) = THETA;
        SAVE.VOL(II) = VOL;
        SAVE.P(II) = Y(1);
        SAVE.TB(II) = Y(2);
        SAVE.TU(II) = Y(3);
        SAVE.X(II) = X;
        SAVE.NOX(II,:) = [ Y(6+1), Y(round(6+0.25*NNOX)),Y(round(6+0.5*NNOX)), Y(round(6+0.75*NNOX)), Y(6+NNOX) ]*to_ppm;
        II=II+1;
        if ( THETAS >= THETA && THETAS < THETAE )
            Y(2) = tinitial( Y(1), Y(3), PHI, F );
        end
        if ( THETA > THETAS + THETAB )
            Y(3) = nan;
        end
    end
    fprintf('%7.1f %6.1f %3.3f %6.1f %6.1f %6.1f %5.0f %5.0f %5.3f %5.2f %6.1f\n', THETA, VOL*1000000, X, Y(1), Y(2), Y(3), Y(4)*1000, Y(5)*1000, M*1000, Y(6)*1000, Y(7)*to_ppm );
end
% integrate total NOx value
NOX_ppm = 0;
for nn=1:NNOX
    THETA = THETAS + (nn-1)/(NNOX-1)*THETAB;
    dxbdtheta = 0.5*sin(pi*(THETA-THETAS)/THETAB)*pi/THETAB;
    dxb = dxbdtheta*DTHETA;
    NOX_ppm = NOX_ppm + Y(6+nn)*dxb*to_ppm;
end
ETA = Y(4)/MNOT*(1+PHI*FS*(1-F))/PHI/FS/(1-F)/A0;
IMEP = Y(4)/(pi/4*B^2*S);
fprintf('ETA=%1.4f IMEP=%7.3f kPa NOx = %6.1f ppm\n',ETA,IMEP,NOX_ppm );
if ( nargin == 0 )
    % if not called externally with custom PHI, F, and RPM parameters,
    % generate theta v. NOx, theta v. temp, theta v. pressure, and theta v.
    % burn fraction
    sTitle = sprintf('Homogenous 2 zone, gasoline, PHI=%.2f F=%.2f RPM=%.1f\nETA=%.3f IMEP=%.2f kPa NOx=%.1f ppm ', PHI, F, RPM, ETA, IMEP, NOX_ppm );
%     figure()
%     plot( SAVE.THETA, SAVE.X, 'linewidth',2 );
%     set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]);
%     xlabel( '\theta','fontsize',18);
%     ylabel('burn fraction','fontsize',18);
    
%     figure()
%     plot( SAVE.THETA, SAVE.P,'linewidth',2 );
%     set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]);
%     xlabel( '\theta','fontsize',18);
%     ylabel('pressure (kPa)','fontsize',18);
%     
    figure()
    plot( SAVE.THETA, SAVE.TU, '-',SAVE.THETA, SAVE.TB,'--','linewidth',2 );
    set(gca,'fontsize',18,'linewidth',2,'Xlim',[-180 180]);
    xlabel( '\theta','fontsize',18);
    ylabel( 'temperature (K)', 'fontsize',18);
    legend('Unburned','Burned', 'Location', 'SouthEast');
    
%     figure()
%     plot( SAVE.THETA, SAVE.NOX,'linewidth',2 );
%     set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]);
%     xlabel('\theta','fontsize',18);
%     ylabel('NOx (ppm)','fontsize',18);
%     axis( [ THETAS, 110, 0, max(max(SAVE.NOX)*1.1) ] );
%     legend( 'X=0', 'X=0.25', 'X=0.5', 'X=0.75', 'X=1', 'Location','SouthEast' );
%     title( sTitle );
end
    function [ TB ] = tinitial( P, TU, PHI, F )
        TB = 2000;
        [~, HU,~, ~, ~, ~, ~, ~, ~, ~] = farg( TU, P, PHI, F, fuel_type );
        for ITER=1:50
            [ierr, ~, HB,~, ~, ~, ~, CP, ~, ~, ~] = ecp( TB, P, PHI, fuel_type );
            if ( ierr ~= 0 )
                fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr )
            end
            DELT = +(HU-HB)/CP;
            TB = TB + DELT;
            if ( abs(DELT/TB) < 0.001 )
                break;
            end
        end
    end
    function [ VOL, X, EM ] = auxiliary( THETA )
        VTDC = pi/4*B^2*S/(R-1); % m3
        VOL = VTDC*(1 + (R-1)/2*(1-cosd(THETA) + 1/EPS*(1-sqrt(1- ...
            (EPS*sind(THETA))^2))));
        X = 0.5*(1-cos(pi*(THETA-THETAS)/THETAB));
        if ( THETA <= THETAS )
            X = 0.;
        end
        if ( THETA >= THETAS+THETAB )
            X = 1.;
        end
        EM = exp(-BLOWBY*(THETA*pi/180 + pi)/OMEGA);
    end
    function [Y] = integrate( THETA, THETAE, Y )
        [TT, YY ] = ode23( @rates, [ THETA, THETAE ], Y );
        for J=1:NY
            Y(J) = YY(length(TT),J);
        end
        function [ YPRIME ] = rates( THETA, Y )
            YPRIME = zeros(NY,1);
            [ VOL, X, EM ] = auxiliary( THETA );
            M = EM*MNOT;
            DUMB = sqrt(1-(EPS*sind(THETA))^2);
            DV = pi/8*B^2*S*sind(THETA)*(1+EPS*cosd(THETA)/DUMB);
            AA = (DV + VOL*BLOWBY/OMEGA)/M;
            C1 = HEAT*(pi*B^2/2 + 4*VOL/B)/OMEGA/M/1000;
            C0 = sqrt(X);
            P = Y(1);
            TB = Y(2);
            TU = Y(3);
            % three different computations are required depending upon the size
            % of the mass fraction burned
            if ( X > 0.999 )
                % EXPANSION
                [ierr,YB,HL, ~, ~,VB, ~,CP, ~,DVDT,DVDP]=ecp(TB,P,PHI,fuel_type);
                if ( ierr ~= 0 )
                    fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr );
                end
                BB = C1/CP*DVDT*TB*(1-TW/TB);
                CC = 0;
                DD = 1/CP*TB*DVDT^2 + DVDP;
                EE = 0;
                YPRIME(1) = (AA + BB + CC)/(DD + EE);
                YPRIME(2) = -C1/CP*(TB-TW) + 1/CP*DVDT*TB*YPRIME(1);
                YPRIME(3) = 0;
            elseif ( X > 0.001 )
                % COMBUSTION
                [~,HU, ~, ~,VU, ~,CPU, ~,DVDTU,DVDPU]=farg(TU,P,PHI,F,fuel_type);
                [ierr,YB,HB, ~, ~,VB, ~,CPB, ~,DVDTB,DVDPB]=ecp(TB,P,PHI,fuel_type);
                if ( ierr ~= 0 )
                    fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr );
                end
                BB = C1*(1/CPB*TB*DVDTB*C0*(1-TW/TB) + 1/CPU*TU*DVDTU*(1-C0)* ...
                    (1-TW/TU));
                DX = 0.5*sin( pi*(THETA-THETAS)/THETAB )*180/THETAB;
                CC = -(VB-VU)*DX - DVDTB*(HU-HB)/CPB*(DX-(X-X^2)*BLOWBY/OMEGA);
                DD = X*(VB^2/CPB/TB*(TB/VB*DVDTB)^2 + DVDPB);
                EE = (1-X)*(1/CPU*TU*DVDTU^2 + DVDPU);
                HL = (1-X^2)*HU + X^2*HB;
                YPRIME(1) = (AA + BB + CC)/(DD + EE);
                YPRIME(2) = -C1/CPB/C0*(TB-TW) + 1/CPB*TB*DVDTB*YPRIME(1) + ...
                    (HU-HB)/CPB*(DX/X - (1-X)*BLOWBY/OMEGA);
                YPRIME(3) = -C1/CPU/(1+C0)*(TU-TW) + 1/CPU*TU*DVDTU*YPRIME(1);
            else
                % COMPRESSION
                [~, HL, ~, ~, ~, ~,CP, ~,DVDT,DVDP]=farg(TU,P,PHI,F,fuel_type);
                BB = C1*1/CP*TU*DVDT*(1-TW/Y(3));
                CC = 0;
                DD = 0;
                EE = 1/CP*TU*DVDT^2 + DVDP;
                YPRIME(1) = ( AA + BB + CC )/(DD + EE);
                YPRIME(2) = 0;
                YPRIME(3) = -C1/CP*(Y(3)-TW) + 1/CP*Y(3)*DVDT*YPRIME(1);
            end
            % common to all cases
            YPRIME(4) = Y(1)*DV;
            YPRIME(5) = 0;
            if ( ~isnan(TB) )
                YPRIME(5) = YPRIME(5) + C1*M*C0*(TB-TW);
            end
            if ( ~isnan(TU) )
                YPRIME(5) = YPRIME(5) + C1*M*(1-C0)*(TU-TW);
            end
            YPRIME(6) = BLOWBY*M/OMEGA*HL;
            % perform NOx integration for each element burned
            if ( X > 0.001 )
                % COMBUSTION OR EXPANSION
                for k=1:NNOX
                    if ( THETA >= THETAS + (k-1)/(NNOX-1)*THETAB )
                        % convert Y(6+k) to [NO] mol/cm^3 from mass fraction
                        % and then back
                        YPRIME(6+k) = zeldovich( TB, P/100, YB, Y(6+k)/(MW_NO*VB*1000) ) ...
                            *MW_NO*VB*1000/OMEGA;
                    end
                end
            end
            % 1/omega is s/rad, so convert to s/deg
            for JJ=1:NY
                YPRIME(JJ) = YPRIME(JJ)*pi/180;
            end
        end
    end
    function [ dNOdt ] = zeldovich( T, P, y, NO )
        % calculate rate of NO formation d[NO]/dt given
        % inputs:
        % T [K] : gas mixture temperature, kelvin
        % P [bar] : cylinder pressure, bar
        % y [...] : equilibrium mole fraction of constituents
        % NO [mol/cm^3] : current NOx concentration
        % outputs:
        % dNOdt [ (mol/cm^3) / sec ] : rate of NO formation
        % extended zeldovich rate constants from Heywood Table 11.1 (cm^3/mol-s)
        k1 = 7.6*10^13*exp(-38000/T);
        k2r = 1.5*10^9*T*exp(-19500/T);
        k3r = 2*10^14*exp(-23650/T);
        % calculate molar concentration [mol/cm^3]
        N_V = (100000*P)/(8.314*T)*(1/100)^3;
        N2e = y(3)*N_V;
        He = y(7)*N_V;
        Oe = y(8)*N_V;
        NOe = y(10)*N_V;
        R1 = k1*Oe*N2e;
        R2 = k2r*NOe*Oe;
        R3 = k3r*NOe*He;
        alpha = NO/NOe;
        dNOdt = 2*R1*(1-alpha*alpha)/(1+alpha*R1/(R2+R3));
    end
TB=SAVE.TB;
TU=SAVE.TU;
T=[TU(1:183),TB(196:360)];
T_interpolate=interp1([-180:2,16:180],T,-180:180,'pchip');
T_smooth=smooth(T_interpolate);


figure()
plot(-180:180,T_interpolate)
xlabel('Theta')
ylabel('T_interpolate')

figure()
plot(-180:180,T_smooth)
xlabel('Theta')
ylabel('T_smooth')

end