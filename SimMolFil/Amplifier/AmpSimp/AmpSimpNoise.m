function [Eout,gain] = AmpSimpNoise(Ein,GssdB,PoutsatdB,NF) 
% simple model of optical amplifier. The model includes the gain
% saturation with noise 
%
% source code by Lam Quoc Huy
% modified by CJH

% Amplifier parameters: 
% small signal gain: GssdB (dB) 
% output saturation power: PoutsatdB (dBm) 
%
% The input is a column vector containing block N samples of the optical signal sampling at the 
% rate 1/Ts
% The output is calculated using 
% Eout = Ein*sqrt(G)
% where: G is the saturated gain 
%       G = Gss*exp(-(G-1)Pin/Psat) (eq1) 
global Ts 
f = 193.1e12;
hplank = 6.6261*1e-34; 
Gss = 10^(GssdB/10);
Poutsat = (10^(PoutsatdB/10))*1e-3; 
Psat = Poutsat*(Gss-2)/Gss/log(2);
% Pinsat = 2* Poutsat/Gss;
N = size(Ein,1);
% Pin = (sum(Ein.*conj(Ein))/N);
Pin = mean( (Ein.*conj(Ein)) );
% numerical calculation of G from the equation G = (Gss - lnG)*Psat/Pin + 1 
tol = 0.05;     % tolerance for G calculation 
step = Gss/2; 
G = Gss; 
err = 10;
while (err > tol)
    G1 = Gss*exp(-(G-1)*Pin/Psat); 
    err = G1 - G;
    if err>0 
        if step <0 
            step = -step/2;
        end
    else
        if step >0
            step = -step/2;
        end
        err = -err;
    end
    G = G + step; 
end
G = G - step;
% Eout = sqrt(G)*Ein; 
gain = G;
Egain = sqrt(G)*Ein; 
dt = Ts;
Bsim = 1/dt;
FigNoise = 10^(NF/10); 
nsp = (FigNoise*G-1)/(2*(G-1)); 
% Pase = hplank.*opfreq.*nsp*(OGain-1)*Bsim 
Pase = hplank.*f.*nsp*(G-1)*Bsim/1000; 
PasedB = 10*log10(Pase);
% afout = fft(Egain) + (randn(size(Egain))+i*randn(size (Egain)) )*sqrt(Pase)./sqrt(2);
% Eout = ifft(afout);
% Eout = Egain + (randn(size(Egain))+i*randn(size (Egain)) )*sqrt(Pase)./sqrt(2)./1; 
Eout = Egain + wgn(N,1,PasedB,'complex');
% Eout = Egain;