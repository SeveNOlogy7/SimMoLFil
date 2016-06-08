clear
% INPUT PARAMETERS *****************************************************
c = 299792.458;                      	%speed of ligth nm/ps

% Input Field Paramenters
N2 = 1^2;                             	% Soliton Order
tfwhm = 200e-3;                     	% ps
lamda_pulse = 1550;                  	% pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

% fiber module parameters
mod.Aeff = 84.95;                       % Effective mode area (um^2)
mod.n2 = 30;                            % Kerr coefficient (10^-16*cm^2/W)
mod.gamma = 2*pi*mod.n2/lamda_pulse/mod.Aeff*1e4;	% W^-1 * km^-1
mod.alpha = 0;                        	% atenuation coef. (km^-1)
mod.L = 0.001;                       	% fiber length (km)
mod.betaw = [0 0 -20.62 120.7e-6];      % beta coefficients (ps^n/nm)
% mod.gssdB = 30;                         % small signal gain coefficient(dB)
% mod.PsatdBm = inf;                     	% saturation input power(dBm)
% mod.lamda_gain = lamda_pulse;         	% central wavelength of gain (nm)
% mod.landa_bw = 100;                   	% gain bandwidth (FWHM, nm)
% mod.fc = c/mod.lamda_gain;              % central frequency of gain (THz)
% mod.fbw = c/(mod.lamda_gain)^2*mod.landa_bw;    % gain bandwidth (THz)

mod.raman = 0;                          % exclude raman effect
mod.ssp = 0;                            % disable self sharpen effect

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
% P_peak = 1.1*N2*abs(mod.betaw(3))/mod.gamma/tfwhm^2;
P_peak = N2*abs(mod.betaw(3))/mod.gamma/tfwhm^2;    % Peak power of the initial pulse (W)   	
u0 =sqrt(P_peak)*sech(t/tfwhm);  %initial field shape in W^0.5
PeakPower = max(u0.^2);

fprintf(1, '\n\n\nInput Peak Power (W) = ');
fprintf(1, '%5.2f%', PeakPower );
fprintf(1, '\n\n');

b = dt*sum(abs(u0.^2));
fprintf(1, '\n\n\nInput Pulse Energy in pJ = ');
fprintf(1, '%5.2f%', b );
fprintf(1, '\n\n');

% uncomment next codes to exclude gain
% if isfield(mod,'gssdB')
%     mod = rmfield(mod,'gssdB');
% end

% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf(1,'\n\nInteraction Picture Method started');
tic

% uncomment one of the following four lines according to the method youwant
% to use
[u,nf,Plotdata] = IP_CQEM_FD(u0,dt,dz,mod,fo,tol,1);
% [u,nf,Plotdata] = SSFM_with_Raman(u0,dt,dz,mod,fo,1);
% u = UPM_SSFM(u0,dt,dz,mod.L/dz,mod.alpha,mod.betaw,mod.gamma,tol);

% uncomment this line to apply bandpass filter to exclude the SBS component
% u = filter_gauss(u,10,c/928.7,1,fo,df); 

tx = toc;

fprintf(1, '\n\nSimulation lasted (s) = ');
fprintf(1, '%5.2f%', tx );

% PLOT RESULTS ************************************************************
figure(79),surf(t,Plotdata.z,abs(Plotdata.u).^2);
colorbar;axis tight;shading interp;
ylabel ('x (km)');
xlabel ('t (ps)');
zlabel ('|u(z,t)|^2 (W)');
% axis([-0.5,0.5,0,max(Plotdata.z),0,max(max(abs(Plotdata.u).^2))])
view(0,90);
 
spec = abs(Plotdata.ufft').^2;
specnorm = spec ./(lambda'*ones(1,size(spec,2))).^2;
% specnorm = specnorm./(ones(size(spec,1),1)*max(specnorm));
specnorm = specnorm/max(specnorm(:));
% figure(80),waterfall(c./(f + fo),Plotdata.z,specnorm');
figure(80),surf(c./(f + fo),Plotdata.z,specnorm');
colorbar;axis tight;
shading interp;
ylabel ('x (km)');
xlabel ('lambda (nm)');
zlabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');
axis([1500,1600,0,max(Plotdata.z),0,1])
view(0,90);