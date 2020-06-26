function [Eout,gain] = AmpSimpNonoise(Ein,L,gssdB,PsatdBm) 
% simple model of optical amplifier. The model includes the gain
% saturation without noise 
%
% source code by Lam Quoc Huy
% modified by CJH
%
% Inputs: 
% small signal gain: GssdB (dB) 
% input saturation power: PsatdBm (dBm) 
% active fiber length: L 
% The input is a column vector containing block N samples of the optical signal sampling at the 
% rate 1/dt
% The output is calculated using
% Eout = Ein*sqrt(G)
% where: G is the saturated gain 

N = size(Ein,1); 
Pin = (sum(Ein.*conj(Ein))/N);
gain = gain_saturated2(Pin,gssdB,PsatdBm);

Eout = sqrt(exp(gain*L))*Ein; 
