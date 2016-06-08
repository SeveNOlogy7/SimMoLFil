clear all
format long e
fprintf(1,'\n\n\n----------------------------------------------');
fprintf(1,'\nSimulating Soliton Propagation in PCF'); 

% INPUT PARAMETERS *****************************************************

c = 299792.458;                      %speed of ligthnm/ps

% Input Field Paramenters
A = 2.5;                             % Soliton Order
tfwhm = 35*1e-3;                     % ps
ni = 1/tfwhm;                        % ps^-1
duration = 1/ni;                     % ps
lamda_pulse = 830;                   % pulse central lambda 
fo=c/lamda_pulse;                    % central pulse frequency (Thz)

% Fiber Parameters
d=33.92;                             % D parameter ps/nmn/km
beta2=-lamda_pulse^2*d/(2*pi*c);     % ps^2/km
beta3 = 7.336283185840707e-02;       % medido en lamda_pulse nm [ps3/km]
gamma = 78;                          % W^-1 * km^-1
alpha = 0;                           % atenuation coef. (km^-1)
L = 0.00075;                         % fiber length (km)
betaw = [0 0 beta2 beta3];          % vector with beta0, beta1, beta2,...

% Numerical Parameters
nt = 2^12;                          % number of spectral points
time = 8;                           % ps
dt = time/nt;                       % ps
t = -time/2:dt:(time/2-dt);         % ps

df=1/(nt*dt);                       % frequencies separation (Thz)
f=-(nt/2)*df:df:(nt/2-1)*df;        % frequencies vector (en THz)
lambda = c./(f + c/lamda_pulse); 	% lambdas vector (nm)
w = 2*pi*f;                         % angular frequencies vector (en THz)

dz = 0.000001;                      % longitudinal step (km)
nz =  L/dz;                         % number of longitudinal steps

% INPUT FIELD ************************************************************
u0 = (A*ni*sqrt(abs(beta2)/gamma)*sech(ni*t)); %initial field shape in W^0.5
PeakPower = max(u0.^2);

fprintf(1, '\n\n\nInput Peak Power (W) = ');
fprintf(1, '%5.2f%', PeakPower );
fprintf(1, '\n\n');
% ************************************************************************

% PR0PAGATE finding numerical solution **********************************
fprintf(1,'\n\nProapagation using the SSFM started');
tx=cputime;
% u = UPM_SSFM(u0,dt,dz,nz,alpha,betaw,gamma,1e-5); % to see the result
% without Raman uncomment here and comment the following line
u = SSFM_with_Raman(u0,dt,dz,nz,alpha,betaw,gamma,fo);

tol = 1e-4; % photon error number, or local error, depending on the used method. 
% [u,nf,Plotdata] = IP_CQEM_FD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol,1);

fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);
fprintf(1,'----------------------------------------------');
fprintf(1,'\n');

% PLOT RESULTS ************************************************************
fprintf(1,'\n\nPloting Results');

%---------plot temporal input and output pulses
figure(1);
plot (t,abs(u0).^2,'*-',t,abs(u).^2,'o-');axis tight;
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Initial (blue) and Final (green) Pulse Shapes');

%---------plot output spectrum------------------

spec = fftshift(abs(fft(u)).^2);
specnorm = spec ./lambda.^2;
specnorm = specnorm/max(specnorm);

figure(2)
plot(c./(f + fo),specnorm,'r.-');axis tight;
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');

fprintf(1,'\n\n----------------------------------------------');

%---------plot input spectrum------------------

spec = fftshift(abs(fft(u0)).^2);
specnorm = spec ./lambda.^2;
specnorm = specnorm/max(specnorm);

figure(3)
plot(c./(f + fo),specnorm,'r.-');axis tight;
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Input Spectrum');

fprintf(1,'\n\n----------------------------------------------');

% *********  Filtered Pulse  ******************************

% filtering windows to see output soliton
df=1/(nt*dt);
f=-(nt/2)*df:df:(nt/2-1)*df; %vector frecuencia (en THz)
lambda = c./(f + c/830);

ventw = ones(1,nt);
filter_cut_at = 870; % wavelengths shorter than this will be filtered
iv = find(lambda<filter_cut_at); 
ventw(iv) = 0;
ventw = fftshift(ventw)';

ufft = fft(u).*ventw';

% plot filtered spectrum
figure(4)
aux = fftshift( abs(ufft).^2)./lambda.^2;
aux = aux / max(aux);
plot(lambda,aux);axis tight;
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Filtered Spectrum');

u = ifft(ufft);

b = dt*sum(abs(u).^2);

fprintf(1, '\n\n\nPulse Energy in pJ = ');
fprintf(1, '%5.2f%', b );
fprintf(1, '\n\n');

% plot ouput soliton after filering spectrum
figure(5)
plot(t,fftshift(abs(u).^2/max(abs(u).^2)));axis tight;
grid on;
xlabel ('t (ps)');
ylabel ('Normalized Pulse Power');
title ('Ouput Soliton');