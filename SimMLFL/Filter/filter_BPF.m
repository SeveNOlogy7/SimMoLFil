function uo = filter_BPF( ui,mod,fo,df )
%FILTER_BPF filter the pulse by a band-pass filter

Ui = fft(ui);
N = size(Ui,2);

c = 299792.458;                      	% speed of light nm/ps

f = fftshift((-(N/2)*df:df:(N/2-1)*df) + fo); 	% frequencies vector (THz)
lambda = c./(f + fo);    	% lambdas vector (nm)

delta_phi = (2*pi./lambda*1e9)*mod.L*mod.bf;
Tf = cos(delta_phi/2).^2;

uo = ifft(Ui.*Tf);       % apply BP filter

end

