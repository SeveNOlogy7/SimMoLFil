% Calculate the small signal gain coefficient and saturation power
% Modify pump power here
pump_power = 8;           % pump power (W)
% Modify the active fiber parameter in the YDF10125 section

%% Pump power

%% Simulation Parameter 
c = 299792.458;                      	% speed of light nm/ps

% Input Field Paramenters
N2 = 1^2;                             	% Soliton Order
tfwhm = 30;                             % ps
lamda_pulse = 1030;                  	% pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

% Numerical Parameters
nt = 2^12;                              % number of spectral points
time = 100;                            	% ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps

df=1/(nt*dt);                           % frequencies separation (Thz)
f=-(nt/2)*df:df:(nt/2-1)*df;            % frequencies vector (en THz)
lambda = c./(f + c/lamda_pulse);    	% lambdas vector (nm)
w = 2*pi*f;                             % angular frequencies vector (in THz)
dz = 0.00001;                        	% longitudinal step (km)
tol = 1e-6; % photon error number, or local error, depending on the used method.
N_trip = 300;
% N_trip = 1;
% if set to 1, all figures will be plotted after simulation hence the 
% computation time of the last trip would be considerable
plot_all_fig = 0;
reinject_flag = 0;

%% Component Parameter 
% Predefine of fiber parameter
PM980.Amod = pi*3.47^2;                  % Mode Field Area (um^2) (um^2)
PM980.n2 = 2.6;                      	% Kerr coefficient (10^-16*cm^2/W)
PM980.gamma = 2*pi*PM980.n2/lamda_pulse/PM980.Amod*1e4;	% W^-1 * km^-1
PM980.alpha = 0;                       	% atenuation coef. (km^-1)
PM980.betaw = [0 0 24.5864 26.1949e-3];	% beta coefficients (ps^n/nm)

PM980.raman = 0;                      	% 0 = exclude raman effect
PM980.ssp = 0;                        	% 0 = disable self sharpen effect

PM1025 = PM980;
PM1025.Amod = pi*4.97^2;
PM1025.gamma = 2*pi*PM1025.n2/lamda_pulse/PM1025.Amod*1e4;	% W^-1 * km^-1
PM1025.betaw = [0 0 20.8634 34.9621e-3];

F10125 = PM1025;
F10125.Amod = pi*5.69^2;
F10125.gamma = 2*pi*F10125.n2/lamda_pulse/F10125.Amod*1e4;	% W^-1 * km^-1
F10125.betaw = [0 0 20.1814 36.8057e-3];
% F10125.betaw = [0 0 20.1814];

%% YDF10125
YDF10125 = F10125;
YDF10125.Amod = pi*5.69^2;                   % Mode Field Area (um^2)
YDF10125.gamma = 2*pi*YDF10125.n2/lamda_pulse/YDF10125.Amod*1e4;	% W^-1 * km^-1
YDF10125.betaw = [0 0 20.1814 36.8057e-3];
% YDF10125.betaw = [0 0 20.1814];

YDF10125.a = 10;                % Core Diameter(um)                   
YDF10125.b = 125;               % Cladding Diameter(um)
YDF10125.absorb_p_0 = 4.8;      % Cladding Absorption(dB/m) @lambda_p_0
YDF10125.lambda_p_0 = 976;      % Pump lambda of Cladding Absorption(dB/m)
% YDF10125.absorb_s = 0.02;   	% Core Attenuation(dB/m)

YDF10125.tf = 0.8e-3;          	% Upper level life time(s)

YDF10125.lambda_p = 976;                % Pump lambda (nm)
YDF10125.lambda_s = lamda_pulse;      	% Signal lambda (nm)

YDF10125 = AmpParaCal(YDF10125,nt,fo,df,c*1e3);    % Calculation of other parameters

amf1 = YDF10125;
amf1.L = 0.001;
sigE_w = GetGuassianSpectrum_w(nt,fo,df,c);
amf1.gain_w = sigE_w/max(sigE_w(:));
amf1.pump_power = pump_power;           % pump power (W)

[gssdb,PsatdBm ] = EqualAmpSimpPara(amf1,amf1.pump_power)
