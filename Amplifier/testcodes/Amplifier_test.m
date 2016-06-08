clear
% INPUT PARAMETERS *****************************************************
c = 299792.458;                             %speed of ligth nm/ps

% Input Field Paramenters
A = 2.5;                                % Soliton Order
tfwhm = 35*1e-3;                        % ps
ni = 1/tfwhm;                         	% ps^-1
duration = 1/ni;                        % ps
lamda_pulse = 830;                      % pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

% Active fiber module parameters
mod.gamma = 78;                     	% W^-1 * km^-1
mod.alpha = 0;                        	% atenuation coef. (km^-1)
mod.L = 0.00075;                       	% fiber length (km)
mod.betaw = [0 0 -12.405 7.3363e-2];  	% beta coefficients (ps^n/nm)
% mod.gssdB = 30;                         % small signal gain coefficient(dB)
% mod.PsatdBm = 40;                     	% saturation input power(dBm)
% mod.lamda_gain = lamda_pulse;         	% central wavelength of gain (nm)
% mod.landa_bw = 150;                   	% gain bandwidth (FWHM, nm)
% mod.fc = c/mod.lamda_gain;              % central frequency of gain (THz)
% mod.fbw = c/(mod.lamda_gain)^2*mod.landa_bw;   % gain bandwidth (THz)

% Numerical Parameters
nt = 2^12;                              % number of spectral points
time = 10;                            	% ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps

df=1/(nt*dt);                           % frequencies separation (Thz)
f=-(nt/2)*df:df:(nt/2-1)*df;            % frequencies vector (en THz)
lambda = c./(f + c/lamda_pulse);    	% lambdas vector (nm)
w = 2*pi*f;                             % angular frequencies vector (en THz)

dz = 0.0000005;                        	% longitudinal step (km)

% INPUT FIELD ************************************************************
u0 = 1*(A*ni*sqrt(abs(mod.betaw(3))/78)*sech(ni*t)); %initial field shape in W^0.5
% randn('state',0);
% u0 = wgn(nt,1,10)';
PeakPower = max(u0.^2);

fprintf(1, '\n\n\nInput Peak Power (W) = ');
fprintf(1, '%5.2f%', PeakPower );
fprintf(1, '\n\n');

b = dt*sum(abs(u0.^2));
fprintf(1, '\n\n\nInput Pulse Energy in pJ = ');
fprintf(1, '%5.2f%', b );
fprintf(1, '\n\n');

% uncomment next line to exclude gain
% mod = rmfield(mod,'gssdB');

% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf(1,'\n\nInteraction Picture Method started');
tol = 1e-3; % photon error number, or local error, depending on the used method. 
tic

% uncomment one of the following four lines according to the method youwant
% to use
[u,nf,Plotdata] = IP_CQEM_FD(u0,dt,dz,mod,fo,tol,1);
% [u,nf,Plotdata] = SSFM_with_Raman(u0,dt,dz,mod,fo,1);

% uncomment this line to apply bandpass filter to exclude the SBS component
% u = filter_gauss(u,10,c/928.7,1,fo,df); 

tx = toc;

fprintf(1, '\n\nSimulation lasted (s) = ');
fprintf(1, '%5.2f%', tx );

% PLOT RESULTS ************************************************************
plotsimp
% plotall
