clear
% INPUT PARAMETERS *****************************************************
c = 299792.458;                      	%speed of ligth nm/ps

% Input Field Paramenters
N2 = 1^2;                             	% Soliton Order
tfwhm = 170e-3;                     	% ps
lamda_pulse = 1550;                  	% pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

% smf1 module parameters
smf1.Aeff = 84.95;                       % Effective mode area (um^2)
smf1.n2 = 30;                            % Kerr coefficient (10^-16*cm^2/W)
smf1.gamma = 2*pi*smf1.n2/lamda_pulse/smf1.Aeff*1e4;	% W^-1 * km^-1
smf1.alpha = 0;                        	% atenuation coef. (km^-1)
smf1.L = 0.001;                       	% fiber length (km)
smf1.betaw = [0 0 -20.62 120.7e-6];      % beta coefficients (ps^n/nm)

smf1.raman = 0;                          % exclude raman effect
smf1.ssp = 0;                            % disable self sharpen effect

% smf2 module parameters
smf2 = smf1;
smf2.L = 0.002-smf1.L;

% NRPS parameters
PhaseShift = pi/8;

% Numerical Parameters
nt = 2^11;                              % number of spectral points
dt = 3.5451e-3;                         % ps
time = dt*nt;                        	% ps
t = -time/2:dt:(time/2-dt);             % ps

df=1/(nt*dt);                           % frequencies separation (Thz)
f=-(nt/2)*df:df:(nt/2-1)*df;            % frequencies vector (en THz)
lambda = c./(f + c/lamda_pulse);    	% lambdas vector (nm)
w = 2*pi*f;                             % angular frequencies vector (en THz)

dz = 0.00001;                        	% longitudinal step (km)
tol = 1e-5; % photon error number, or local error, depending on the used method. 

% INPUT FIELD ************************************************************
P_peak = 2*N2*abs(smf1.betaw(3))/smf1.gamma/tfwhm^2;    % Peak power of the initial pulse (W)   	
u0 =sqrt(P_peak)*sech(t/tfwhm);  %initial field shape in W^0.5
PeakPower = max(u0.^2);

fprintf(1,'\n----------------------------------------------\n');
fprintf('Input Peak Power (W) = %5.2f\n',PeakPower);

b = dt*sum(abs(u0.^2));
fprintf('Input Pulse Energy in pJ = %5.2f\n', b );

% uncomment next codes to exclude gain
% if isfield(mod,'gssdB')
%     mod = rmfield(mod,'gssdB');
% end

% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf('\nInteraction Picture Method started\n');
tic

% uncomment one of the following four lines according to the method youwant
% to use
rho = 0.30;
[uf,ub] = coupler(u0,0,rho);
% forward light cp1o-smf1-NRPS-smf2-cp2o
[ufo,nf,Plotdata_f] = IP_CQEM_FD(uf,dt,dz,smf1,fo,tol,1);
ufo = ufo.*exp(i*PhaseShift);
[ufo,nf,Plotdata_f] = IP_CQEM_FD(ufo,dt,dz,smf2,fo,tol,1);
% forward light cp1o-smf1-NRPS-smf2-cp2o
[ubo,nf,Plotdata_b] = IP_CQEM_FD(ub,dt,dz,smf2,fo,tol,1);
ubo = ubo.*exp(-i*PhaseShift);
[ubo,nf,Plotdata_b] = IP_CQEM_FD(ubo,dt,dz,smf1,fo,tol,1);

[ur,ut] = coupler(ubo,ufo,rho);
u = ut;

% uncomment this line to apply bandpass filter to exclude the SBS component
% u = filter_gauss(u,10,c/928.7,1,fo,df); 

tx = toc;

fprintf('\nSimulation lasted (s) = %5.2f%\n', tx );

% PLOT RESULTS ************************************************************
% plotsimp
plotSA
% Plotdata = Plotdata_f;
% plotall