clear
% INPUT PARAMETERS *****************************************************
c = 299792.458;                      	%speed of ligth nm/ps

% Input Field Paramenters
N2 = 1^2;                             	% Soliton Order
tfwhm = 50;                             % ps
lamda_pulse = 1030;                  	% pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

% smf1 module parameters
smf1.Aeff = 42.6488;                   	% Effective mode area (um^2)
smf1.n2 = 30;                        	% Kerr coefficient (10^-16*cm^2/W)
smf1.gamma = 2*pi*smf1.n2/lamda_pulse/smf1.Aeff*1e4;	% W^-1 * km^-1
smf1.alpha = 0;                        	% atenuation coef. (km^-1)
smf1.L = 0.0016;                       	% fiber length (km)
smf1.betaw = [0 0 17.8777 39.4876e-6]; 	% beta coefficients (ps^n/nm)

smf1.raman = 0;                      	% exclude raman effect
smf1.ssp = 0;                        	% disable self sharpen effect

% amf module parameters
amf1 = smf1;
amf1.L = 0.00015;
amf1.gssdB = 30;                            % small signal gain coefficient(dB)
amf1.PsatdBm = 30;                          % saturation input power(dBm)
amf1.lamda_gain = lamda_pulse;              % central wavelength of gain (nm)
amf1.landa_bw = 100;                        % gain bandwidth (FWHM, nm)
amf1.fc = c/amf1.lamda_gain;                % central frequency of gain (THz)
amf1.fbw = c/(amf1.lamda_gain)^2*amf1.landa_bw;     % gain bandwidth (THz)

% smf2 module parameters
smf2 = smf1;
smf2.L = 0.002-smf1.L-amf1.L;

% smf3 module parameters
smf3 = smf1;
smf3.L = 0.0005;

% smf4 module parameters
smf4 = smf1;

% smf5 module parameters
smf5 = smf1;
smf5.L = 0.0005;

% amf2 module parameters
amf2 = amf1;
amf2.gssdB = 36;                         % small signal gain coefficient(dB)
amf2.L = 0.001;

% filter parameters
filter.lamda_c = lamda_pulse;
filter.landa_bw = 40;
filter.fc = c/filter.lamda_c;
filter.f3dB = c/(filter.lamda_c)^2*filter.landa_bw;
filter.n = 1;

% coupler parameters
rho = 0.40;             % NALM
rho_out = 0.95;         % Output coupler

% Numerical Parameters
nt = 2^12;                              % number of spectral points
time = 200;                             % ps
dt = time/nt;                           % ps
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
randn('state',0);
u0 = wgn(nt,1,20)'.*u0;
PeakPower = max(u0.^2);

fprintf(1,'\n----------------------------------------------\n');
fprintf('Input Peak Power (W) = %5.2f\n',PeakPower);

b = dt*sum(abs(u0.^2));
fprintf('Input Pulse Energy in pJ = %5.2f\n', b );

% uncomment next codes to exclude gain
% if isfield(amf,'gssdB')
%     amf = rmfield(amf,'gssdB');
% end

% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf('\nInteraction Picture Method started\n');
tic

spec_z = [];
u_z = [];

u = u0;
h1 = waitbar( 0,'The program is running...');
N_trip = 120;
for ii = 1:N_trip
    waitbar( (ii-1)/N_trip,h1);
    [u,nf,Plotdata] = IP_CQEM_FD(u,dt,dz,smf4,fo,tol,1,1);
    
    [u,nf,Plotdata] = IP_CQEM_FD(u,dt,dz,amf2,fo,tol,1,1);
    
    [u,nf,Plotdata] = IP_CQEM_FD(u,dt,dz,smf5,fo,tol,1,1);
    
    [uf,ub] = coupler(u,0,rho);
    % forward light cp1o-smf1-amf-smf2-cp2o
    [ufo,nf,Plotdata] = IP_CQEM_FD(uf,dt,dz,smf1,fo,tol,1,1);
    [ufo,nf,Plotdata] = IP_CQEM_FD(ufo,dt,dz,amf1,fo,tol,1,1);
    [ufo,nf,Plotdata] = IP_CQEM_FD(ufo,dt,dz,smf2,fo,tol,1,1);
    % backward light cp2o-smf2-amf-smf1-cp1o
    [ubo,nf,Plotdata] = IP_CQEM_FD(ub,dt,dz,smf2,fo,tol,1,1);
    [ubo,nf,Plotdata] = IP_CQEM_FD(ubo,dt,dz,amf1,fo,tol,1,1);
    [ubo,nf,Plotdata] = IP_CQEM_FD(ubo,dt,dz,smf1,fo,tol,1,1);
    
    [ur,ut] = coupler(ubo,ufo,rho);
    u = ut;
    
    % ut->smf3->filter->output coupler->smf4->amf2->smf5->NALM
    [u,nf,Plotdata] = IP_CQEM_FD(u,dt,dz,smf3,fo,tol,1,1);
    
    %     u = filter_gauss(u,filter.f3dB,filter.fc,filter.n,fo,df);
    
    [u,uout] = coupler(u,0,rho_out);
    
    %---------plot temporal input and output pulses
    figure(1);
    plot (t,abs(u0).^2,'b.-',t,abs(uout).^2,'r.-');axis tight;
    grid on;
    xlabel ('t (ps)');
    ylabel ('|u(z,t)|^2 (W)');
    title ('Initial (blue) and Final (green) Pulse Shapes');
    %---------plot output spectrum------------------
    spec = fftshift(abs(fft(uout)).^2);
    specnorm = spec ./lambda.^2;
    specnorm = specnorm/max(specnorm);
    figure(3)
    plot(c./(f + fo),specnorm,'r.-');axis tight;
    grid on;
    xlabel ('lambda (nm)');
    ylabel ('Normalized Spectrum (a.u.)');
    title ('Output Spectrum');
    
    spec_z = [spec_z;specnorm];
    u_z = [u_z;uout];
    
end
close(h1)
tx = toc;

fprintf('\nSimulation lasted (s) = %5.2f%\n', tx );

% PLOT RESULTS ************************************************************
figure(4),surf(t,1:N_trip,abs(u_z).^2);
colorbar;axis tight;shading interp;
ylabel ('run trip');
xlabel ('t (ps)');
zlabel ('|u(z,t)|^2 (W)');
title ('Output Optical Field Evolution');
% axis([-0.5,0.5,0,max(Plotdata.z),0,max(max(abs(Plotdata.u).^2))])
view(0,90);

figure(5),surf(c./(f + fo),1:N_trip,spec_z);
colorbar;axis tight;
shading interp;
ylabel ('run trip');
xlabel ('lambda (nm)');
zlabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum Evolution');
% axis([1500,1600,0,max(Plotdata.z),0,1])
view(0,90);

phase_out = phase(uout);
delta_w = diff(phase_out);
Eout = abs(uout).^2;
[uout_max,I_uout_max] = max(Eout);
[width,I_l,I_r] = fwhm(Eout);
range = I_uout_max-4*width:I_uout_max+4*width;
figure(6);[ax,p1,p2] = plotyy(t,Eout,t(range),delta_w(range),'plot','plot');
grid on; axis tight;
xlabel ('t (ps)');
ylabel (ax(1),'|u(z,t)|^2 (W)');
ylabel(ax(2),'Chirp(THz)') % label right y-axis