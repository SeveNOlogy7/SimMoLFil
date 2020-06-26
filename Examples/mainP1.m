% This example reproduces Point 1 in Figure 2 of the following paper:
% Cai, J.-H., Chen, H., Chen, S.-P. & Hou, J. Opt. Express 25, 4414 (2017).

%% Structural Parameters to investigate
lambda_bw = 11; 
smf5_L = 0.025;
NOLM_L = 0.0015;
PsatdBm = 23.5;                         % use cal_AmpSimpPara.m to calculate from pump power
smf4_L = 0.002;

%% Simulation Parameter 
c = 299792.458;                      	% speed of light nm/ps
lambda_pulse = 1030;                  	% pulse central lambda (nm)
fo=c/lambda_pulse;                       % central pulse frequency (THz)

% Numerical Parameters
nt = 2^12;                              % number of spectral points
time = 120;                            	% ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps

df=1/(nt*dt);                           % frequencies separation (Thz)
f=-(nt/2)*df:df:(nt/2-1)*df;            % frequencies vector (en THz)
lambda = c./(f + c/lambda_pulse);    	% lambdas vector (nm)
w = 2*pi*f;                             % angular frequencies vector (in THz)
dz = 0.00001;                        	% Initial longitudinal step (km)
tol = 1e-7;
N_trip = 300;

% Some helpful settings
% if set to 1, all figures will be plotted after simulation hence the 
% computation time of the last trip would be considerable
plot_all_fig = 1;
reinject_flag = 0;

%% Set Component Parameters
% Predefine of fiber parameters
PM980.Amod = pi*3.3^2;                  % Mode Field Area (um^2) (um^2)
PM980.n2 = 2.6;                      	% Kerr coefficient (10^-16*cm^2/W)
PM980.gamma = 2*pi*PM980.n2/lambda_pulse/PM980.Amod*1e4;	% W^-1 * km^-1
PM980.alpha = 0;                       	% atenuation coef. (km^-1)
PM980.betaw = [0 0 24.5864 26.1949e-3];	% beta coefficients (ps^n/nm)

PM980.raman = 0;                      	% 0 = exclude raman effect
PM980.ssp = 0;                        	% 0 = disable self sharpen effect

PM1025 = PM980;
PM1025.Amod = pi*5^2;
PM1025.gamma = 2*pi*PM1025.n2/lambda_pulse/PM1025.Amod*1e4;	% W^-1 * km^-1
PM1025.betaw = [0 0 20.8634 34.9621e-3];

F10125 = PM980;
F10125.Amod = pi*5.25^2;
F10125.gamma = 2*pi*F10125.n2/lambda_pulse/F10125.Amod*1e4;	% W^-1 * km^-1
F10125.betaw = [0 0 20.1814 36.8057e-3];

% smf1 module parameters
NOLM.smf1 = F10125;
NOLM.smf1.L = 0.0007;                       	% fiber length (km)

% smf2 module parameters
NOLM.smf2 = F10125;
NOLM.smf2.L = NOLM_L-NOLM.smf1.L;

% smf3 module parameters
smf3 = F10125;
smf3.L = 0.001;

% smf4 module parameters
smf4 = F10125;
smf4.L = smf4_L; 

% smf5 module parameters
smf5 = F10125;
smf5.L = smf5_L; 

% smf5 module parameters
smf6 = F10125;
smf6.L = 0.001;

% amf1 module parameters
amf1 = F10125;
amf1.L = 0.0008;
amf1.gssdB = 37.5;          % small signal gain coefficient(dB)
amf1.PsatdBm = PsatdBm; 	% saturation input power(dBm)

amf1.fbw = c/(1030)^2*80; 
amf1.fc = c/1030;

% filter parameters
filter.lambda_c = 1030;
filter.lambda_bw = lambda_bw;
filter.fc = c/filter.lambda_c;
filter.f3dB = c/(filter.lambda_c)^2*filter.lambda_bw;
filter.n = 7;

% coupler parameters
NOLM.rho = 0.30;          
rho_out = 0.40;         % Output ratio = 1-rho_out
PhaseShift = 0; % Nonreciprocal Phase Shift for NOLM

Length_cavity = smf4.L+smf5.L+amf1.L+NOLM.smf1.L+NOLM.smf2.L;   % km
RepeatFre = c/2/Length_cavity;                               	% Hz
amf1.RepeatFre = RepeatFre;

%% INPUT FIELD
u0 = rand_sech(nt,time);

if ~reinject_flag
    u = u0;
end

fprintf(1,'\n----------------------------------------------\n');
fprintf('Input Peak Power (W) = %5.2f\n',max(u0.*conj(u0)));
fprintf('Input Pulse Energy in pJ = %5.2f\n', dt*sum(u0.*conj(u0)));

%% PROPAGATE finding numerical solution
fprintf('\nLaser Simulation started\n');
tic

spec_z = zeros(N_trip,nt);
u_z = spec_z;

h1 = waitbar( 0,'The program is running...');

dplot = 0;
dstop = 0;
dconv = 0;

for ii = 1:N_trip
    if plot_all_fig == 1 && ii == N_trip
        dplot = 1;
    end
    %u->BPF->SMF4->AMF1->SMF5->NOLM->OC->u
    u_f = filter_gauss(u,filter.f3dB,filter.fc,filter.n,fo,df);
    u = u_f;
    
    [u,~,Plotdata_smf4] = IP_CQEM_FD(u,dt,dz,smf4,fo,tol,dplot,1);
    [u,~,Plotdata_amf1] = IP_CQEM_FD(u,dt,dz,amf1,fo,tol,dplot,0); 
    [u,~,Plotdata_smf5] = IP_CQEM_FD(u,dt,dz,smf5,fo,tol,dplot,0);
    
    [uf,ub] = coupler(u,0,NOLM.rho);
    % forward light cp1o-smf1-NRPS-smf2-cp2o
    [ufo,~,Plotdata_smf1_f] = IP_CQEM_FD(uf,dt,dz,NOLM.smf1,fo,tol,dplot,1);
%     ufo = ufo.*exp(1i*PhaseShift);
    [ufo,~,Plotdata_smf2_f] = IP_CQEM_FD(ufo,dt,dz,NOLM.smf2,fo,tol,dplot,1);
    % forward light cp1o-smf1-NRPS-smf2-cp2o
    [ubo,~,Plotdata_smf1] = IP_CQEM_FD(ub,dt,dz,NOLM.smf2,fo,tol,dplot,1);
%     ubo = ubo.*exp(-1i*PhaseShift);
    [ubo,~,Plotdata_smf2] = IP_CQEM_FD(ubo,dt,dz,NOLM.smf1,fo,tol,dplot,1);
    
    [ur,ut] = coupler(ubo,ufo,NOLM.rho);

    u = ut;  
    [u,uout] = coupler(u,0,rho_out);    
    
    %---------plot temporal input and output pulses
    figure(1);
    plot (t,abs(uout).^2,'r.-');axis tight;
    grid on;
    xlabel ('t (ps)');
    ylabel ('|u(z,t)|^2 (W)');
    title ('Pulse Shape');
    %---------plot output spectrum------------------
    spec = fftshift(fft(uout));
    specnorm = spec.*conj(spec);
    specnorm = specnorm/max(specnorm);
    figure(3)
    plot(c./(f + fo),specnorm,'r.-');axis tight;
    grid on;
    xlabel ('lambda (nm)');
    ylabel ('Normalized Spectrum (a.u.)');
    title ('Output Spectrum');
    
    spec_z(ii,:) = specnorm;
    u_z(ii,:) = uout;
    
    if dstop
        break;
    end
   
    % Close the waitbar manually, the simulation will end
    try
        waitbar( ii/N_trip,h1,sprintf('The program is running...%d/%d',ii,N_trip));
    catch err
        if strcmp(err.identifier,'MATLAB:waitbar:InvalidSecondInput')
            fprintf('\nSimulation is going to be terminated\n');
        else
            fprintf(err.message);
        end       
        if plot_all_fig ==1 % if figures are required 
            dstop = 1;
            dplot = 1;
        else
            break;
        end
    end 
    
    % Check terminating condition
    if ii>5
        if max(xcorr_normal(spec_z(ii,:),spec_z(ii-1,:),nt/8))>0.99999
            dconv = dconv + 1 ;
            if dconv >= 5 
                fprintf('\nSimulation is going to terminate\n');
                if plot_all_fig ==1 % if figures are required 
                    dstop = 1;
                    dplot = 1;
                else
                    break;
                end
            end
        else
            dconv = 0;
        end
    end
    
end
try
    close(h1)
catch err
end

tx = toc;
fprintf('\nSimulation lasted (s) = %5.2f%\n', tx );

fprintf('\nOutput Peak Power (W) = %5.2f\n',max(uout.*conj(uout)));

PulseEnergy = dt*sum(uout.*conj(uout));
fprintf('Output Pulse Energy in pJ = %5.2f\n', PulseEnergy);


%% Plot results

phase_out = phase(uout);
delta_w = diff(phase_out);
Eout = abs(uout).^2;
[uout_max,I_uout_max] = max(Eout);
[width,~,~] = fwhm(Eout);
[~,I_l,I_r] = fwzm(Eout,1e-4);
range = I_l:I_r;
figure(6);[ax,p1,p2] = plotyy(t,Eout,t(range),delta_w(range),'plot','plot');
grid on; axis(ax(1),'tight');
axis(ax(2),[min(t),max(t),-2,2]);
xlabel (['t (ps)',sprintf(', %.1fps',width*dt),...
    sprintf(', %.1fW(Peak)',PulseEnergy/width/dt),...
    sprintf(', %.1fMHz',RepeatFre/1e6)]);
ylabel (ax(1),'|u(z,t)|^2 (W)');
ylabel(ax(2),'Chirp(THz)') % label right y-axis

%% More plots
if dplot
    Plot_ML_result
end
