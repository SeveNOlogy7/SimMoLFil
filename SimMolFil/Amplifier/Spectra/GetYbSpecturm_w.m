function [ sigA_w,sigE_w,alpha_w ] = GetYbSpecturm_w( nt,fo,df,c )
%GETYBSPECTURM_W_ return the emission and absoption spectra of Yb

% Input parameters
%  	nt:
%   fo: central pulse frequency (THz)
%   df: frequencies separation (THz)
%   c:
% Output parameters
%   

f = (-(nt/2)*df:df:(nt/2-1)*df) + fo;           % frequencies vector (THz)
lam = c./f;                                     % nm
[sigA_w,sigE_w,alpha_w] = GetYbSpectrum(lam);
% plot(lam,[sigA_w;sigE_w;])

end

