function [ gssdb,PsatdBm ] = EqualAmpSimpPara(mod,P_pump)
% Retrive small signal gain coefficient and saturation power by numerically 
% solving the amplifier rate equation

% gssdb :db/km

% gss = pump_p*mod.absorb_p*mod.P_IS_p/mod.P_CS_s/(mod.P_IS_p+pump_p)-mod.absorb_s; % m^-1
% gss = gss*1e3
% gssdb = 10*log10(gss);

P_signal = 0.00001;             % W, small signal
L = 5;                          % m

dz = 1e-2;
L = dz:dz:L;

gain = zeros(1,length(L));
P_p_z = zeros(1,length(L));
P_s_z = zeros(1,length(L));

P_p_z(1) = P_pump;
P_s_z(1) = P_signal;
for ii = 1:length(L)
    [gain_t,~,P_p_z_t,P_s_z_t] = AmpGain(mod,P_p_z(ii),P_s_z(ii),dz,'forward');
    gain(ii) = gain_t;
    P_p_z(ii+1) = P_p_z_t;
    P_s_z(ii+1) = P_s_z_t;  
end
temp = abs(gain - max(gain)/2);
PsatdBm = 10*log10(P_s_z(temp == min(temp))*1e3);
gssdb = 10*log10(max(gain));


end

