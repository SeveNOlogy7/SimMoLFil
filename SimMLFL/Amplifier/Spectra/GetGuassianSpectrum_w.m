function [ sigE_w ] = GetGuassianSpectrum_w( nt,fo,df,c )

f = (-(nt/2)*df:df:(nt/2-1)*df) + fo;           % frequencies vector (THz)
lam = c./f;                                     % nm

sigE_w = 1e-27* (340*exp(-(((lam-1030)/100).^2)));


end

