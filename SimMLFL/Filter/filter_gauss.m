function uo = filter_gauss(ui,f3dB,fc,n,fo,df) 
% filter_gauss(ui,f3dB,n)
% filter the input signal with the n order gaussian filter 
%
% Input parameters
%   ui: input field amplitude (row vector)
%   f3dB: gaussian filter 3dB bandwith
%   fc: gaussian filter central frequency (THz)
%   n: gaussian filter order
%   fo: central pulse frequency (THz)
%   df: frequencies separation (THz)
% 
% Outputs
%   uo: output field amplitude
%

Ui = fft(ui);
N = size(Ui,2);

f = fftshift((-(N/2)*df:df:(N/2-1)*df) + fo); 	% frequencies vector (THz)

% n order gaussian filter
Tf = exp(-log(sqrt(2))*(2/f3dB*(f-fc)).^(2*n)); 

uo = ifft(Ui.*Tf);       % apply n order gaussian filter
