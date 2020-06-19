function [ gain_w ] = AmpGainWCal( mod, gain,dz)
%
%       gain (km-1)

% c = 299792458;       	% speed of light (m/s)
% h = 6.626e-34;          % Planck Constant (J/s)

% N2 = gain*1e-3/mod.sigE_s/mod.dopN;
% if N2<0.1
%     N2 = 0.1;
% end
% N0 = 0.8*(1- N2);
% N0 = 1- N2;
% gain_w = mod.sigE_w*N2 - mod.sigA_w*N0;

n2 = gain*1e-3/(mod.k0*mod.Gamma_s*mod.dopN)/(mod.sigE_s+mod.sigA_s)+...
    mod.sigA_s/(mod.sigE_s+mod.sigA_s);
n0 = 1 - n2;
gain_w = mod.sigE_w*n2 - mod.sigA_w*n0;

end