function gain = gain_saturated2(Pin,gssdB,PsatdBm) 
% calculate the gain coefficient of the amplifier given the input power, the small
% signal gain coefficient and saturation power
%
% Input parameters
%   Pin: input average power (W)
%   gssdB: small signal gain coefficient(dB/km)
%   PsdBm : saturation input power(dBm)
% 
% Outputs
%   gain: gain saturated gain coefficient (km^-1)
%
% source code by CJH

% where: G is the saturated gain 
%       G = Gss/(1+Pin/Psat)

gss = 10^(gssdB/10);
Psat = (10^(PsatdBm/10))/1000; 

gain = gss./(1+Pin./Psat);  %(km-1) ############