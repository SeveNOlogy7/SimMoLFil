clear
clc

%% General Parameters
c = 299792458;       	% speed of light (m/s)
h = 6.626e-34;          % Planck Constant (J/s)

lamda_pulse = 1030;                  	% pulse central lambda (nm)
fo=c/lamda_pulse;                       % central pulse frequency (THz)

nt = 2^11;                              % number of spectral points
time = 60;                            	% ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps

df=1/(nt*dt);                           % frequencies separation (Thz)

%% Active Fiber Parameters
% Predefine of fiber parameter
PM980.Amod = pi*3.47^2;                  % Mode Field Area (um^2) (um^2)
PM980.n2 = 2.6;                      	% Kerr coefficient (10^-16*cm^2/W)
PM980.gamma = 2*pi*PM980.n2/lamda_pulse/PM980.Amod*1e4;	% W^-1 * km^-1
PM980.alpha = 0;                       	% atenuation coef. (km^-1)
PM980.betaw = [0 0 24.5864 26.1949e-3];	% beta coefficients (ps^n/nm)

PM980.raman = 0;                      	% 0 = exclude raman effect
PM980.ssp = 0;                        	% 0 = disable self sharpen effect

% From Manufacturer DataSheet
YDF10125 = PM980;
YDF10125.a = 10;               % Core Diameter(um)                   
YDF10125.b = 125;               % Cladding Diameter(um)
YDF10125.Amod = pi*5.69^2;                   % Mode Field Area (um^2)
YDF10125.absorb_p_0 = 4.8;      % Cladding Absorption(dB/m) @lambda_p_0
YDF10125.lambda_p_0 = 976;      % Pump lambda of Cladding Absorption(dB/m)
% YDF10125.absorb_s = 0.02;       % Core Attenuation(dB/m)
YDF10125.betaw = [0 0 20.1814 36.8057e-3];

YDF10125.tf = 0.8e-3;          	% Upper level life time(s)

YDF10125.lambda_p = 976;         % Pump lambda (nm)
YDF10125.lambda_s = lamda_pulse;        % Signal lambda (nm)

YDF10125 = AmpParaCal(YDF10125,nt,fo,df,c);    % Calculation of other parameters

P_pump = 10;                   % W
P_signal = 1e-7;             % W
L = 5;                          % m

dz = 1e-2;
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
temp = abs(gain - max(gain)/2);
PsatdBm = 10*log10(P_s_z(temp == min(temp))*1e3)
gssdb = 10*log10(max(gain))

gain2 = gain_saturated2(P_s_z(1:end-1),gssdb,PsatdBm);

gain_s1 = diff(P_s_z)./P_s_z(1:end-1)/dz;
gain_s2 = diff(P_s_z)./P_s_z(2:end)/dz;
gain_s = (gain_s1+gain_s2)/2;

plot(L,gain);hold on;
plot(L,gain2);hold off
