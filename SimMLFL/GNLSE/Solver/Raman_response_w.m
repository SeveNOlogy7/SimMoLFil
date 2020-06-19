function [ hrw,fr] = Raman_response_w( t,mod)
%RAMAN_RESPONSE
global Raman;
if isfield(mod,'raman') && mod.raman==0
    hrw = 0;
    fr = 0;
elseif isempty(Raman) || length(Raman.hrw)~=length(t)
    % Raman parameters
    t1 = 12.2e-3;                   % raman parameter t1 [ps]
    t2 = 32e-3;                     % raman parameter t2 [ps]
    tb = 96e-3;                     % ps
    fc = 0.04;
    fb = 0.21;
    fa = 1 - fc - fb;
    fr = 0.245;                     % fraccion de respuesta retardada Raman
    
    tres = t-t(1);                  % time starting in 0
    
    ha =((t1^2+t2^2)/(t1*t2^2)).*exp(-tres/t2).*sin(tres/t1);
    hb = ((2*tb - tres)./tb^2).*exp(-tres/tb);
    hr = (fa + fc)*ha + fb*hb;      % hr(t) Raman responce function (ps^-1)
    
    hrw = fft(hr);
    
    Raman.hrw = hrw;
    Raman.fr = fr;
else 
    hrw = Raman.hrw;
    fr = Raman.fr;
end

end

