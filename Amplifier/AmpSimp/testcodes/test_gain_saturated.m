% Amplifier parameters: 
GssdB = 10;          % (dB) 
PsatdB = 40;      % (dBm) 
NF = 7;              % (dB)

Pin = linspace(0,20,100);

gain2 = [];
for ii = 1:100
    gain2 = [gain2,gain_saturated2(Pin(ii),GssdB,PsatdB)]; 
end
figure(2);plot(Pin,gain2);

gain3 = [];
for ii = 1:100
    gain3 = [gain3,gain_saturated3(Pin(ii),GssdB,PsatdB)]; 
end
figure(3);plot(Pin,gain3);
