Pin =1;

GssdB = 20;          % (dB) 
PoutsatdB = 40;      % (dBm) 

Gss = 10^(GssdB/10);
Poutsat = (10^(PoutsatdB/10))/1000; 
 
G0 = Gss/100;
G = fzero(@(G)(G-Gss*exp(-G+1)),G0);

Psat = Poutsat/G;

G = fzero(@(G)(G-Gss*exp(-(G-1)*Pin/Psat)),Gss);

