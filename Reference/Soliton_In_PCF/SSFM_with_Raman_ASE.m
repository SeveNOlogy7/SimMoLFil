function u1 = SSFM_with_Raman_ASE(u0,dt,dz,nz,alpha,betap,gamma,fo)

% This function solves the nonlinear Schrodinger equation for
% pulse propagation in an optical fiber using the split-step
% Fourier method and a fourth-order RK method for the non-linear part
%
% INPUT
%
% u0 - starting field amplitude (vector)
% dt - time step
% dz - propagation stepsize
% nz - number of steps to take, ie, ztotal = dz*nz
% alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
% betap - dispersion polynomial coefs, [beta_0 ... beta_m], or beta(w)
% gamma - nonlinearity coefficient
% fo - central frequency of the simulation (THz)
%
% OUTPUT
%
% u1 - field at the output
%
% NOTES  The dimensions of the input and output quantities can
% be anything, as long as they are self consistent.  E.g., if
% |u|^2 has dimensions of Watts and dz has dimensions of
% meters, then gamma should be specified in W^-1*m^-1.
% Similarly, if dt is given in picoseconds, and dz is given in
% meters, then beta(n) should have dimensions of ps^(n-1)/m.

nt = length(u0);                            % number of sample points
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);  % angular frequencies
t = -(nt/2)*dt:dt:(nt/2-1)*dt;              %vector temporal (en ps)

% Raman parameters and hr(t)
t1 = 12.2e-3;                 % raman parameter t1 [ps]
t2 = 32e-3;                   % raman parameter t2 [ps]
tb = 96e-3;                   % ps
fc = 0.04;
fb = 0.21;
fa = 1 - fc - fb;
fr = 0.245;                    % fraccion de respuesta retardada Raman
tres = t-t(1);                % time starting in 0

ha =((t1^2+t2^2)/(t1*t2^2)).*exp(-tres/t2).*sin(tres/t1); %funcion respuesta Raman (ps^-1)
hb = ((2*tb - tres)./tb^2).*exp(-tres/tb); 
hr = (fa + fc)*ha + fb*hb;
hrw = fft(hr);

%ASE nopise term
hrwi = imag(hrw);
UH = zeros(1,nt);
UH(nt/2:1:nt)=1;
plank = (6.62606896*10^-34)*(10^12)/(2*pi); % J*ps
kT = 300*1.3806504*10^-23; % J
boltz = (exp(plank*abs(w)/(kT) - 1)).^-1;
ASEdp = (2*fr*plank*2*pi*fo/gamma)*abs(hrwi).*(boltz' + UH);

% constructing linear operator
linearoperator = -alpha/2;
if (length(betap) == nt)     % If the user manually specifies beta(w)
    linearoperator = linearoperator - 1i*betap;
    linearoperator = fftshift(linearoperator);
    else
    for ii = 0:length(betap)-1;
        linearoperator = linearoperator - 1i*betap(ii+1)*(w).^ii/factorial(ii);
    end
    linearoperator = conj(linearoperator');
end


% Performig the SSFM ***********************************************
fprintf(1, '\nSimulation running...      ');
ufft = fft(u0);
fiberlength = nz*dz;
propagedlength = 0;

while propagedlength < fiberlength,
      
  if (dz + propagedlength) > fiberlength,
      dz = fiberlength - propagedlength;
  end    
  
  % First linear part, until dz/2
  halfstep = exp(linearoperator*dz/2);
  uhalf = ifft(halfstep.*ufft);

  % Non-Linear part
  ASEterm = sqrt(ifft(ASEdp)*10^12).*randn(1,nt);
  
  integral_convoluida = ifft( hrw.*fft( abs(uhalf).^2 ) );
  k1 = -1i*gamma*(uhalf.*(1-fr).*abs(uhalf).^2 + dt*uhalf.*fr.*integral_convoluida + 1i*uhalf.*ASEterm);
  k1 = k1 - (gamma/(2*pi*fo))*(1/dt)*gradient(uhalf.*(1-fr).*abs(uhalf).^2 + dt*uhalf.*fr.*integral_convoluida + 1i*uhalf.*ASEterm ) ;
  
  uhalf2 = uhalf + dz*k1/2;
  integral_convoluida = ifft( hrw.*fft( abs(uhalf2).^2 ) );
  k2 = -1i*gamma*(uhalf2.*(1-fr).*abs(uhalf2).^2 + dt*uhalf2.*fr.*integral_convoluida + 1i*uhalf.*ASEterm);
  k2 = k2 - (gamma/(2*pi*fo))*(1/dt)*gradient(uhalf2.*(1-fr).*abs(uhalf2).^2 + dt*uhalf2.*fr.*integral_convoluida + 1i*uhalf.*ASEterm ) ;

  uhalf3 = uhalf + dz*k2/2;
  integral_convoluida = ifft( hrw.*fft( abs(uhalf3).^2 ) );
  k3 = -1i*gamma*(uhalf3.*(1-fr).*abs(uhalf3).^2 + dt*uhalf3.*fr.*integral_convoluida + 1i*uhalf.*ASEterm );
  k3 = k3 - (gamma/(2*pi*fo))*(1/dt)*gradient(uhalf3.*(1-fr).*abs(uhalf3).^2 + dt*uhalf3.*fr.*integral_convoluida + 1i*uhalf.*ASEterm) ;

  uhalf4 = uhalf + dz*k3;
  integral_convoluida = ifft( hrw.*fft( abs(uhalf4).^2 ) );
  k4 = -1i*gamma*(uhalf4.*(1-fr).*abs(uhalf4).^2 + dt*uhalf4.*fr.*integral_convoluida + 1i*uhalf.*ASEterm );
  k4 = k4 - (gamma/(2*pi*fo))*(1/dt)*gradient(uhalf4.*(1-fr).*abs(uhalf4).^2 + dt*uhalf4.*fr.*integral_convoluida + 1i*uhalf.*ASEterm) ;
  
  u1 = uhalf + dz*(k1 + 2*k2 + 2*k3 + k4)/6;% + gamma* dz*nt*uhalf.*ifft(sqrt(ASEdp).*randn(1,nt));
  
  % Second linear part (from dz/2 to dz)
  ufft = halfstep.*fft(u1);
  
  propagedlength = propagedlength + dz;
  fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /fiberlength );
  
end

% output parameter
u1 = ifft(ufft);


