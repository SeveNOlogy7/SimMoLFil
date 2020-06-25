function [u0] = rand_sech(nt, time_window, beta2,gamma)
%RAND_SECH Generate a rand sech profile as the initial condition
% Initial condtion will not affect the stable solution
% Init the light field to a pulse-like shape would reduce convergence time,
% compaired with plain random numbers

if(nargin == 2)
    beta2 = 25;
    gamma = 3;
end

dt = time_window/nt;
t = -time_window/2:dt:(time_window/2-dt);

N2 = 1^2;                             	% Soliton Order
tfwhm = time_window/4;                  	% ps
P_peak = 2*N2*abs(beta2)/gamma/tfwhm^2;	% Peak power of the initial pulse (W)
u0 =sqrt(P_peak)*sech(t/tfwhm);         %initial field shape in W^0.5
rng(0);
u0 = abs(wgn(nt,1,25))'.*u0;

end

