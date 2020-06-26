clear
clc
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
% if set to 1, all figures will be plotted after simulation hence the 
% computation time of the last trip would be considerable
plot_all_fig = 1;

%% Component Parameter 
% Predefine of fiber parameter
PM980.Amod = pi*3.3^2;                  % Mode Field Area (um^2) (um^2)
PM980.n2 = 2.6;                      	% Kerr coefficient (10^-16*cm^2/W)
PM980.gamma = 2*pi*PM980.n2/lamda_pulse/PM980.Amod*1e4;	% W^-1 * km^-1
PM980.alpha = 0;                       	% atenuation coef. (km^-1)
PM980.betaw = [0 0 28.9311 25.4837e-3];	% beta coefficients (ps^n/nm)

PM980.raman = 0;                      	% 0 = exclude raman effect
PM980.ssp = 0;                        	% 0 = disable self sharpen effect

YDF10125 = PM980;
YDF10125.a = 5.5;               % Core Diameter(um)                   
YDF10125.b = 125;               % Cladding Diameter(um)
YDF10125.absorb_p_0 = 4.8;      % Cladding Absorption(dB/m) @lambda_p_0
YDF10125.lambda_p_0 = 976;      % Pump lambda of Cladding Absorption(dB/m)
YDF10125.absorb_s = 0.02;       % Core Attenuation(dB/m)

YDF10125.tf = 0.8e-3;          	% Upper level life time(s)

YDF10125.lambda_p = 976;            % Pump lambda (nm)
YDF10125.lambda_s = 1030;           % Signal lambda (nm)

YDF10125 = AmpParaCal(YDF10125,nt,fo,df,c*1e3);    % Calculation of other parameters

% YDF10125.pump_power = 20;           % pump power (W)

P_pump = 50;                 % W
P_signal = 0.001;           % W
L = 4;                      % m

dz = 1e-1;
L = dz:dz:L;

gain = zeros(1,length(L));
P_p_z = zeros(1,length(L));
P_s_z = zeros(1,length(L));

P_p_z(1) = P_pump;
P_s_z(1) = P_signal;
for ii = 1:length(L)
    [gain_t,~,P_p_z_t,P_s_z_t] = AmpGain(YDF10125,P_p_z(ii),P_s_z(ii),dz,'forward');
    gain(ii) = gain_t;
    P_p_z(ii+1) = P_p_z_t;
    P_s_z(ii+1) = P_s_z_t;
   
end
%  N2 = gain*1e-3/YDF10125.sigE_s/YDF10125.dopN;
 N2 = gain*1e-3/YDF10125.sigE_s;
 
% N2 = 0.1:0.1:1;
% N0 = 1- N2;
% for ii = 1:length(N2)
%     gain_w(ii,:) = YDF10125.sigE_w*N2(ii) - YDF10125.sigA_w*N0(ii);
% %     figure(1),plot(c./(f+fo),gain_w(ii,:)/max(gain_w(ii,:)));hold on;axis tight
%     figure(1),plot(c./(f+fo),gain_w(ii,:));hold on;axis tight
% end 