%% General Parameters
c = 299792458;       	% speed of light (m/s)
h = 6.626e-34;          % Planck Constant (J/s)

%% Active Fiber Parameters
% From Manufacturer DataSheet
YDF10125.a = 5.5;               % Core Diameter(um)                   
YDF10125.b = 125;               % Cladding Diameter(um)
YDF10125.absorb_p_0 = 4.8;      % Cladding Absorption(dB/m) @lambda_p_0
YDF10125.lambda_p_0 = 976;      % Pump lambda of Cladding Absorption(dB/m)
YDF10125.absorb_s = 0.02;       % Core Attenuation(dB/m)

YDF10125.tf = 0.8e-3;          	% Upper level life time(s)

YDF10125.lambda_p = 976;         % Pump lambda (nm)
YDF10125.lambda_s = 1030;        % Signal lambda (nm)

% Calculation
YDF10125.Aeff = pi*YDF10125.a^2;                % Effective Area (doped) (um^2)
YDF10125.Amod = pi*5.25^2;                      % Mode Field Area (um^2)
YDF10125.Gamma_s = YDF10125.Amod/YDF10125.Aeff; % Confinement Factor for signal
YDF10125.Gamma_p = (YDF10125.a/YDF10125.b)^2;   % Confinement Factor for pump

% Get Emission and Absorbtion Cross Section of Pump and Signal Light
[sigA,sigE,~] = GetYbSpectrum([YDF10125.lambda_p,YDF10125.lambda_s]);
YDF10125.sigE_p = sigE(1);
YDF10125.sigA_p = sigA(1);
YDF10125.sigE_s = sigE(2);
YDF10125.sigA_s = sigA(2);

% Calculate Absorption Coefficient @lambda_p
% Dopant Concentration (m^-3)
k0 = 4.343;                     % Constant, See 10.1364/OE.24.009237
YDF10125.dopN = YDF10125.absorb_p_0/k0/GetYbSpectrum(976)/YDF10125.Gamma_p;
YDF10125.absorb_p = k0*YDF10125.dopN*YDF10125.Gamma_p*YDF10125.sigA_p;  
% Fit Well with the Manufacturer Data 1.6 dB/m @915 nm

% Calculate Intrinsic and Cross Saturation Power
YDF10125.P_IS_p = YDF10125.Aeff*1e-12/YDF10125.Gamma_p/YDF10125.tf/...
    (YDF10125.sigA_p+YDF10125.sigE_p)*h*c/(YDF10125.lambda_p*1e-9);
YDF10125.P_CS_p = YDF10125.P_IS_p;
YDF10125.P_IS_s = YDF10125.Aeff*1e-12/YDF10125.Gamma_s/YDF10125.tf/...
    (YDF10125.sigA_s+YDF10125.sigE_s)*h*c/(YDF10125.lambda_s*1e-9);
YDF10125.P_CS_s = YDF10125.P_IS_s;

P_pump = 10;               % W
P_signal = 0.001;             % W
Z = 4;                      % m

% AmpEquation = @(P,S)([P-P_pump*exp(-10^(-YDF10125.absorb_p/10)*Z+(P_pump-P)...
%     /YDF10125.P_IS_p+(P_signal-S)/YDF10125.P_CS_p);...
%     S-P_signal*exp(-10^(-YDF10125.absorb_s/10)*Z+(P_signal-S)...
%     /YDF10125.P_IS_s+(P_pump-P)/YDF10125.P_CS_s)]);
% AmpEquation = @(X)(AmpEquation(X(1),X(2)));
AmpEquation = @(P,S)([P-P_pump*exp(10^(-YDF10125.absorb_p/10)*Z+(P_pump-P)...
    /YDF10125.P_IS_p-(P_signal-S)/YDF10125.P_CS_p);...
    S-P_signal*exp(-10^(-YDF10125.absorb_s/10)*Z+(P_signal-S)...
    /YDF10125.P_IS_s-(P_pump-P)/YDF10125.P_CS_s)]);
AmpEquation = @(X)(AmpEquation(X(1),X(2)));

Options = optimoptions('fsolve');
Options.Display = 'none';
if Z<1e-5
    Options.TolFun = 1e-6*Z/1e-5;
    Options.TolX = Options.TolFun;
end
[P,fval,exitflag] = fsolve(AmpEquation,[P_pump,P_signal],Options);
P_pump_z = P(1);
P_signal_z = P(2);
g = log(P_signal_z/P_signal)/Z;     % m^-1

[gain,~,P_p_z,P_s_z] = AmpGain(YDF10125,P_pump,P_signal,Z,'forward')


[gain,~,P_p_0,P_s_z] = AmpTBV(YDF10125,P_p_z,P_signal,Z);
% AmpEquation = @(P_pump,S)([P_p_z-P_pump*exp(10^(-YDF10125.absorb_p/10)*Z+(P_pump-P_p_z)...
%     /YDF10125.P_IS_p-(P_signal-S)/YDF10125.P_CS_p);...
%     S-P_signal*exp(-10^(-YDF10125.absorb_s/10)*Z+(P_signal-S)...
%     /YDF10125.P_IS_s-(P_pump-P_p_z)/YDF10125.P_CS_s)]);
% AmpEquation = @(X)(AmpEquation(X(1),X(2)));
% 
% Options = optimoptions('fsolve');
% Options.Display = 'none';
% if Z<1e-5
%     Options.TolFun = 1e-6*Z/1e-5;
%     Options.TolX = Options.TolFun;
% end
% [P,fval,exitflag] = fsolve(AmpEquation,[P_p_z,P_signal],Options);
% 
% P_pump = P(1)
% P_signal_z = P(2)