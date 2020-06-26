function [ gain,exitflag,P_p_z,P_s_z ] = AmpGain( mod,P_p_0,P_s_0,z,dir)
% calculate the gain coefficient of the amplifier induced from Active fiber
% parameter in module descripter, pump and signal optical power
%
% Input parameters
%   mod: amplifier parameters
%   P_p_0 : pump power at the beginning of the active fiber
%   P_s_0 : signal power at the beginning of the active fiber
%   dir: pump direction,'forward'/'backward'
%   z: propogating distance, should be comparatively small (m)
% 
% Outputs
%   gain: gain coefficient (km^-1)
%   exitflag: output of fsolve
%   P_p_z: pump power at z distance after the beginning of the active fiber
%   P_s_z: signal power at z distance after the beginning of the active fiber
%

u_p = 1;
if strcmp(dir,'backward')
    u_p = -1;
end

AmpEquation = @(P,S)([P-P_p_0*exp(u_p*(-10^(-mod.absorb_p/10)*z+u_p*(P_p_0-P)...
    /mod.P_IS_p+(P_s_0-S)/mod.P_CS_p));...
    S-P_s_0*exp(-10^(-mod.absorb_s/10)*z+(P_s_0-S)...
    /mod.P_IS_s+u_p*(P_p_0-P)/mod.P_CS_s)]);
AmpEquation = @(X)(AmpEquation(X(1),X(2)));

Options = optimoptions('fsolve');
Options.Display = 'none';
% if z<1e-4
%     Options.TolFun = 1e-6*z/1e-4;
%     Options.TolX = Options.TolFun;
% end
Options.TolFun = 1e-12;
Options.TolX = Options.TolFun;

% Using bordinary condition as the initial value for the solution,
% for the Distance Z tends to be small,
% though a better way is to forecast the solution
[P,~,exitflag] = fsolve(AmpEquation,[P_p_0,P_s_0],Options);

% Forecasting coarse solution
% maximum gain
% Gmax = exp(z*(10^(-mod.absorb_p/10)/mod.P_CS_s*mod.P_IS_p-10^(-mod.absorb_s/10)));
% P_s_z_coarse = Gmax*P_s_0;
% if u_p == 1
%     P_p_z_coarse = P_p_0-(P_s_z_coarse - P_s_0)*mod.lambda_p/mod.lambda_s;
% elseif u_p == -1
%     P_p_z_coarse = P_p_0+(P_s_z_coarse - P_s_0)*mod.lambda_p/mod.lambda_s;
% end
% 
% [P,~,exitflag] = fsolve(AmpEquation,[P_p_z_coarse,P_s_z_coarse],Options);

P_p_z = P(1);
P_s_z = P(2);

if u_p == 1 && (P_p_z<0 || P_p_z>P_p_0 || P_s_z<0 || (P_s_z+P_p_z)>(P_s_0+P_p_0))
    fprintf('\nAmpGain Calculation failed\n');
elseif u_p == -1 && (P_p_z<0 || P_p_z<P_p_0 || P_s_z<0 || (P_s_z+P_p_0)>(P_s_0+P_p_z))
    fprintf('\nAmpGain Calculation failed\n');
end

if exitflag == -2
    fprintf('Algorithm appears to be converging to a point that is not a root.');
end

gain = log(P_s_z/P_s_0)/z*1e3;     % ;km^-1   
% 'gain' is averaged over z. Only if z <<1, 'gain' can be regarded
% as gain coefficient at this point.

end

