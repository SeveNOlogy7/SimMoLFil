function [ mod_cal ] = AmpParaCal( mod,nt,fo,df,c )
%AMPPARACAL calculate active fiber parameters according to Manufacturer Data

% c = 299792458;       	% speed of light (m/s)
h = 6.626e-34;          % Planck Constant (J/s)

% Calculation
mod.Aeff = pi*mod.a^2/4;                  % Effective Area (doped) (um^2)
mod.Gamma_s = mod.Aeff/mod.Amod;        % Confinement Factor for signal
mod.Gamma_p = (mod.a/mod.b)^2;          % Confinement Factor for pump

% Get Emission and Absorbtion Cross Section of Pump and Signal Light
[sigA,sigE,~] = GetYbSpectrum([mod.lambda_p,mod.lambda_s]);
mod.sigE_p = sigE(1);
mod.sigA_p = sigA(1);
mod.sigE_s = sigE(2);
mod.sigA_s = sigA(2);

% Calculate Absorption Coefficient @lambda_p
% Dopant Concentration (m^-3)
mod.k0 = 4.343;                     % Constant, See 10.1364/OE.24.009237
mod.dopN = mod.absorb_p_0/mod.k0/GetYbSpectrum(mod.lambda_p_0)/mod.Gamma_p;    % dopant concentration m-3
mod.absorb_p = mod.k0*mod.dopN*mod.Gamma_p*mod.sigA_p;  
% Fit Well with the Manufacturer Data 1.6 dB/m @915 nm

mod.absorb_s = mod.k0*mod.dopN*mod.Gamma_s*mod.sigA_s;  

% Calculate Intrinsic and Cross Saturation Power
mod.P_IS_p = mod.Aeff*1e-12/mod.Gamma_p/mod.tf/...
    (mod.sigA_p+mod.sigE_p)*h*c/(mod.lambda_p*1e-9);
mod.P_CS_p = mod.P_IS_p;
mod.P_IS_s = mod.Aeff*1e-12/mod.Gamma_s/mod.tf/...
    (mod.sigA_s+mod.sigE_s)*h*c/(mod.lambda_s*1e-9);
mod.P_CS_s = mod.P_IS_s;

[sigA_w,sigE_w,~] = GetYbSpecturm_w(nt,fo,df,c*1e-3);
[sigA_m,sigE_m,~] = GetYbSpectrum(975);
mod.sigA_w = sigA_w/sigA_m;
mod.sigE_w = sigE_w/sigE_m;
mod.sigA_w = sigA_w;
mod.sigE_w = sigE_w;
mod.gain_w = 0;                 % to indicate gain spectrum evelotion
mod.gssdB = 0;                  % to indicate an active fiber

mod.freq_s = c*1e-3/mod.lambda_s;


mod_cal = mod;

end

