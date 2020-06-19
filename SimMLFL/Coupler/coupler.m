function [ u1o,u2o ] = coupler( u1i,u2i,rho )
%COUPLER simulate an four-end optical coupler
% 
% rho: couple ratio, u1o = sqrt(rho)*u1i, u1i and u1o are at the same side

if rho>1
    rho = 1;
elseif rho <0
    rho = 0;
end

u1o = sqrt(rho)*u1i + 1i*sqrt(1-rho)*u2i;
u2o = 1i*sqrt(1-rho)*u1i + sqrt(rho)*u2i;

end

