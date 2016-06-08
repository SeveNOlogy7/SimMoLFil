clear all
format long e
fprintf(1,'\n\n\n----------------------------------------------');
fprintf(1,'\nSimulating Supercontinuum Generation in PCF');


% INPUT PARAMETERS *****************************************************

c = 299792.458;                     %speed of ligth nm/ps

% Input Field Paramenters
tfwhm = 28.4e-3;                    % ps
ni = 1/tfwhm;                       % ps^-1
lamda_central = 835;
fo=c/lamda_central;                 % central pulse frequency (Thz)

% Fiber Parameters
gamma = 110;                        % W^-1 * km^-1
alpha = 0;                          % atenuation coef. (km^-1)
L = 0.0001;                         % fiber length (km)
betaw = [0 0 -11.830 8.1038e-2 -9.5205e-5 2.0737e-7 -5.3943e-10 1.3486e-12 -2.5495e-15 3.0524e-18 -1.714e-21]; % beta coefficients (ps^n/ nm)


% Numerical Parameters
nt = 2^15;                              % number of spectral points`
time = 32;                              % ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps
dz = 1e-7;                              % initial longitudinal step (km)
v = [(0:nt/2-1),(-nt/2:-1)]/(dt*nt);    % frequencies frequencies (THz)

% INPUT FIELD ***********************************************************
PeakPower = 10000; % W, change here the Input Power!
u0 = sqrt(PeakPower)*sech(ni*t); %initial field shape in W^0.5


% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf(1,'\n\nInteraction Picture Method started');
tol = 1e-2; % photon error number, or local error, depending on the used method. 
tic

% uncomment one of the following four lines according to the method youwant
% to use
[u,nf] = IP_CQEM_FD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol);
% [u,nf] = IP_CQEM_TD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol);
% [u,nf] = IP_LEM_FD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol);
% [u,nf] = IP_LEM_TD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol);


tx = toc;

fprintf(1, '\n\nSimulation lasted (s) = ');
fprintf(1, '%5.2f%', tx );

%  ---------plot output spectrum------------------
specnorm = fftshift(abs(fft(u)).^2);
specnorm = specnorm/max(specnorm);
figure(2)
hold on
plot(c./(fftshift(v) + fo),10*log10(specnorm),'b.-');
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');
axis([480 1550 -70 1])

fprintf(1,'\n----------------------------------------------\n\n');
