clear
% Active fiber parameters: 
gssdB = 10;             % (dB/m) 
PsatdBm = 40;           % (dBm) 
LossdB = 0;             % (dB)
loss = 10^(LossdB/10);
Length =0.5;           	% (m)
d_z = 0.001;       	% (m)
N_z = Length/d_z;

Pmax = 10;
N_P = 100;

Pin = linspace(0,Pmax,N_P);           
Pout = Pin;
tic
for jj = 1:N_P;
    for ii = 1:N_z
        gain2 = gain_saturated2(Pout(jj),gssdB,PsatdBm);
        Pout(jj) = exp((gain2-loss)*d_z)*Pout(jj);
    end
end
toc

figure(1);plot(Pin,Pout);