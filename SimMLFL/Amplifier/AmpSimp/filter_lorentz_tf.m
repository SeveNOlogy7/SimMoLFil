function tf = filter_lorentz_tf(ui,fbw,fc,fo,df) 
% filter_lorentz_tf(fbw,fc,fo,df) 
% calculate the transfer function of lorentz filter
%
% Input parameters
%  	ui: input field amplitude (row vector)
%   fbw: lorentz filter bandwith (FWHM)
%   fc: lorentz filter central frequency (THz)
%   fo: central pulse frequency (THz)
%   df: frequencies separation (THz)
% Output parameters
%   tf: transfer function of lorentz filter
%
% 
% Outputs
% source code by CJH

N = size(ui,2);

% The Ui's corresponding frequency
f = (-(N/2)*df:df:(N/2-1)*df) + fo;      % frequencies vector (THz)
w = 2*pi*f;

% lorentz filter
tf = (fbw)/2/pi./((f-fc).^2+(fbw/2)^2);
tf = tf/max(tf(:));