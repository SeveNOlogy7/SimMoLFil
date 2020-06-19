function [ gain,exitflag,P_p_0,P_s_z ] = AmpTBV( mod,P_p_z,P_s_0,z)
% calculate the gain coefficient of the amplifier induced from Active fiber
% parameter in module descripter, pump and signal optical power
% USING TWO-POINT BOUNDARY VALUE CONDITION
% ONLY BACKWARD PUMPING
%
% Input parameters
%   mod: amplifier parameters
%   P_p_z: pump power at z distance after the beginning of the active fiber
%   P_s_0 : signal power at the beginning of the active fiber
%   z: propogating distance, should be comparatively small
% 
% Outputs
%   gain: gain coefficient (m^-1)
%   exitflag: output of fsolve
%   P_p_0 : pump power at the beginning of the active fiber
%   P_s_z: signal power at z distance after the beginning of the active fiber
%

AmpEquation = @(P_pump,S)([P_p_z-P_pump*exp(10^(-mod.absorb_p/10)*z+(P_pump-P_p_z)...
    /mod.P_IS_p-(P_s_0-S)/mod.P_CS_p);...
    S-P_s_0*exp(-10^(-mod.absorb_s/10)*z+(P_s_0-S)...
    /mod.P_IS_s-(P_pump-P_p_z)/mod.P_CS_s)]);
AmpEquation = @(X)(AmpEquation(X(1),X(2)));

Options = optimoptions('fsolve');
Options.Display = 'none';
if z<1e-4
    Options.TolFun = 1e-6*z/1e-4;
    Options.TolX = Options.TolFun;
end
[P,~,exitflag] = fsolve(AmpEquation,[P_p_z,P_s_0],Options);

P_p_0 = P(1);
P_s_z = P(2);

gain = log(P_s_z/P_s_0)/z;     % m^-1
% 'gain' is averaged over z. Only if z <<1, 'gain' can be regarded
% as gain coefficient at this point.

end

