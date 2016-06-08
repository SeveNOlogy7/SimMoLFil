function gain = gain_saturated3(Pin,GssdB,PsatdBm) 
% calculate the gain of the amplifier given the input power, the small
% signal gain and saturation power 
%
% Pin: input power
% GssdB: small signal gain(dB)
% PsatdBm : saturation input power(dBm)
%
% source code by CJH

% where: G is the saturated gain 
%       G = Gss*exp(-(G-1)Pin/Psat) (eq1)

Gss = 10^(GssdB/10);
Psat = (10^(PsatdBm/10))/1000;

% numerical calculation of G
G = fzero(@(G)(G-Gss*exp(-(G-1)*Pin/Psat)),Gss/10);

gain = G;