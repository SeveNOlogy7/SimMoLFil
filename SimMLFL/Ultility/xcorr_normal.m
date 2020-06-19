function [ out ] = xcorr_normal( x,y,lag )

nt = size(x,1);
out = xcorr(x,y,lag);
out = out/((sqrt(sum(x.^2)*sum(y.^2)))/nt^2);

end

