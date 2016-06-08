clear
% INPUT PARAMETERS *****************************************************
c = 299792.458;                      %speed of ligthnm/ps
% Input Field Paramenters
A = 2.35;                          	 % Soliton Order
tfwhm = 35*1e-3;                     % ps
ni = 1/tfwhm;                        % ps^-1
lamda_pulse = 830;                   % pulse central lambda 
fo=c/lamda_pulse;                    % central pulse frequency (Thz)

% Fiber Parameters                  
d=33.92;                             % D parameter ps/nmn/km
beta2=-lamda_pulse^2*d/(2*pi*c);     % ps^2/km
gamma = 78;                          % W^-1 * km^-1
alpha = 0;                           % atenuation coef. (km^-1)
L = 0.00075;                         % fiber length (km)

% Numerical Parameters
nt = 2^11;                          % number of spectral points
time = 8;                           % ps
dt = time/nt;                       % ps
t = -time/2:dt:(time/2-dt);         % ps

df = 1/(nt*dt);                     % frequencies separation (THz)
f = -(nt/2)*df:df:(nt/2-1)*df;      % frequencies vector (THz)
lambda = c./(f + fo);               % lambdas vector (nm)
w = 2*pi*f;                         % angular frequencies vector (THz)

dz = 0.000001;                      % longitudinal step (km)
nz =  L/dz;                         % number of longitudinal steps

% INPUT FIELD ************************************************************
chirp = 0;
u0 = (A*ni*sqrt(abs(beta2)/gamma)*sech(ni*t)).*exp(1i*chirp*t.^2); %initial field shape in W^0.5
PeakPower = max(u0.^2);

fprintf(1, '\n\n\nInput Peak Power (W) = ');
fprintf(1, '%5.2f%', PeakPower );
fprintf(1, '\n\n');

% CALCULATING beta(w) FROM SAVED beta2(w) ********************************
load Frequencies_for_beta2
load beta2w_Thz;
Frequencies_for_beta2 = fliplr(Frequencies_for_beta2);
beta2w_Thz = fliplr(beta2w_Thz);

% Uncomment the following lines to visualize beta2(w) and beta2(lambda)
% PLOT saved beta2(w) and beta2(lambda)
figure(5)
plot(2*pi*c./Frequencies_for_beta2,beta2w_Thz)
grid on;
xlabel ('lambda (nm)');
ylabel ('beta2 (Trad/s)');
figure(6)
plot(Frequencies_for_beta2,beta2w_Thz)
grid on;
xlabel ('angular frequencies (Thz)');
ylabel ('beta2 (Trad/s)');

betawaux = 0;
beta2waux = 0;

% beta2(w) polinomial fit
poly = polyfit(Frequencies_for_beta2 - 2*pi*fo,beta2w_Thz,9);
poly = fliplr(poly);
for mm=1:1:length(poly),
    beta2waux = beta2waux + poly(mm)*(Frequencies_for_beta2 - 2*pi*fo).^(mm-1);
end

% finding beta(w)
for mm=1:1:length(poly),
    betawaux = betawaux + (poly(mm)/((mm+1)*mm))*(Frequencies_for_beta2 - 2*pi*fo).^(mm+1);
end
betaw = interp1(Frequencies_for_beta2 - 2*pi*fo,betawaux,w,'spline','extrap');
% betaw = betaw; % this is beta(w) used in the propagation
% ************************************************************************

% PR0PAGATE finding numerical solution **********************************
fprintf(1,'\n\nProapagation using the SSFM started');
tx=cputime;
%u = UPM_SSFM(u0,dt,dz,nz,alpha,betaw,gamma,1e-4); % to see the result
% to see the result without the Raman effect uncomment the above line and
% comment the following one
% u = SSFM_with_Raman(u0,dt,dz,nz,alpha,betaw,gamma,fo);

tol = 1e-4; % photon error number, or local error, depending on the used method. 
[u,nf] = IP_CQEM_FD(u0,dt,L,dz,alpha,betaw,gamma,fo,tol);

fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);
fprintf(1,'----------------------------------------------');
fprintf(1,'\n');


% PLOT RESULTS ************************************************************
fprintf(1,'\n\nPloting Results');

%---------plot temporal input and output pulses
figure(1);
plot (t,abs(u0).^2,'*-',t,abs(u).^2,'o-');
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Initial (blue) and Final (green) Pulse Shapes');


%---------plot output spectrum------------------

spec = fftshift(abs(fft(u)).^2);
specnorm = spec ./lambda.^2;
specnorm = specnorm/max(specnorm);

figure(2)
plot(c./(f + fo),specnorm,'r.-')
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
plot(c./(f + fo),specnorm,'r.-')
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
plot(lambda,aux)
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
figure(10)
plot(t,fftshift(abs(u).^2/max(abs(u).^2)))
grid on;
xlabel ('t (ps)');
ylabel ('Normalized Pulse Power');
title ('Ouput Soliton');