function [width,I_l,I_r] = fwzm(x,eps)
% return the Full Width at Zero Maximum of the pulse x 

nsize = length(x); 
[peak,ind_peak] = max(x); 
zero_peak = peak*eps; % 1e-4

x_l = x(1:ind_peak);
x_r = x(1+ind_peak:end);
I_l = find(fliplr(x_l)<=zero_peak,1,'first');
I_r = find(x_r<=zero_peak,1,'first');

if isempty(I_l)
    I_l = find(fliplr(x_r)<=zero_peak,1,'first');
    width = I_l+I_r+ind_peak;
    I_l = nsize - I_l;
elseif isempty(I_r)
    I_r = find(x_l<=zero_peak,1,'first');
    width = I_l+I_r+nsize-ind_peak;
else
    width = I_l+I_r;
    I_l = ind_peak - I_l;
    I_r = ind_peak + I_r;
end