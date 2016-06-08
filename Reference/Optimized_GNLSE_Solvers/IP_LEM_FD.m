function [u1,nf] = IP_LEM_FD(u0,dt,L,dz,alpha,betap,gamma,fo,tol)

% This function solves the Generalized Nonlinear Schrodinger Equation with 
% thwe complete Raman response for pulse propagation in an optical fiber
% using the Interaction Picture Method in combination with the Local Error 
% method for the step-size determination and the frequency domain
% integration of the nonlinear operator
%
% INPUT
%
% u0 - starting field amplitude (vector)
% dt - time step
% L - propagation distance
% dz - initial step size
% alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
% betap - dispersion polynomial coefs, [beta_0 ... beta_m], or beta(w)
% gamma - nonlinearity coefficient
% fo - central frequency of the simulation (THz)
% tol - relative photon error
%
% OUTPUT
%
% u1 - field at the output
% nf - number of FFTs performed


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
fr = 0.245;                   % fraccion de respuesta retardada Raman
tres = t-t(1);                % time starting in 0

ha =((t1^2+t2^2)/(t1*t2^2)).*exp(-tres/t2).*sin(tres/t1);
hb = ((2*tb - tres)./tb^2).*exp(-tres/tb);
hr = (fa + fc)*ha + fb*hb;   %Raman responce function (ps^-1)

hrw = fft(hr);

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
propagedlength = 0;
u1 = u0;
nf = 1;

% Performig the SSFM according to the LEM spatial-step size
while (propagedlength < L),
    if (dz + propagedlength) > L,
        dz = L - propagedlength;
    end
    
    halfstep = exp(linearoperator*dz/2);
    quarterstep=exp(linearoperator*dz/4);
    uhalf = halfstep.*ufft; 
    uquarter = quarterstep.*ufft;
    
    
    %%%%% NON LINEAR OPERATOR COARSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uip = uhalf;
    
    k1 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( u1.*((1-fr)*abs(u1).^2) + fr*dt*u1.*ifft(hrw.*fft( abs(u1).^2 )));
    k1 = halfstep.*k1;
    
    uhalf2 = ifft(uip + k1/2);
    k2 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf2.*((1-fr)*abs(uhalf2).^2) + fr*dt*uhalf2.*ifft(hrw.*fft( abs(uhalf2).^2 )));
    
    uhalf3 = ifft(uip + k2/2);
    k3 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf3.*((1-fr)*abs(uhalf3).^2) + fr*dt*uhalf3.*ifft(hrw.*fft( abs(uhalf3).^2 )));
    
    uhalf4 = ifft(halfstep.*(uip + k3));
    k4 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf4.*((1-fr)*abs(uhalf4).^2) + fr*dt*uhalf4.*ifft(hrw.*fft( abs(uhalf4).^2 )));
    
    uc = halfstep.*(uip + k1./6 + k2./3 + k3./3) + k4./6;
    
    nf = nf + 15;
    %%%%%%%%% END COARSE%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%FIRST NON LINEAR FINE%%%%%%%%%%%%%%%%%%%%%%%%
    uip = uquarter;
    
    k1 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( u1.*((1-fr)*abs(u1).^2) + fr*dt*u1.*ifft(hrw.*fft( abs(u1).^2 )));
    k1 = quarterstep.*k1;
    
    uhalf2 = ifft(uip + k1/2);
    k2 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf2.*((1-fr)*abs(uhalf2).^2) + fr*dt*uhalf2.*ifft(hrw.*fft( abs(uhalf2).^2 )));
    
    uhalf3 = ifft(uip + k2/2);
    k3 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf3.*((1-fr)*abs(uhalf3).^2) + fr*dt*uhalf3.*ifft(hrw.*fft( abs(uhalf3).^2 )));
    
    uhalf4 = ifft(quarterstep.*(uip + k3));
    k4 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf4.*((1-fr)*abs(uhalf4).^2) + fr*dt*uhalf4.*ifft(hrw.*fft( abs(uhalf4).^2 )));
    
    uf = quarterstep.*(uip + k1./6 + k2./3 + k3./3) + k4./6;

    nf = nf + 15;
    
    %%%%%%%SECOND NON LINEAR FINE----------------
    uip = quarterstep.*uf;
    uf = ifft(uf);
    
    k1 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uf.*((1-fr)*abs(uf).^2) + fr*dt*uf.*ifft(hrw.*fft( abs(uf).^2 )));
    k1 = quarterstep.*k1;
    
    uhalf2 = ifft(uip + k1/2);
    k2 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf2.*((1-fr)*abs(uhalf2).^2) + fr*dt*uhalf2.*ifft(hrw.*fft( abs(uhalf2).^2 )));
    
    uhalf3 = ifft(uip + k2/2);
    k3 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf3.*((1-fr)*abs(uhalf3).^2) + fr*dt*uhalf3.*ifft(hrw.*fft( abs(uhalf3).^2 )));
    
    uhalf4 = ifft(quarterstep.*(uip + k3));
    k4 = -0.5*dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf4.*((1-fr)*abs(uhalf4).^2) + fr*dt*uhalf4.*ifft(hrw.*fft( abs(uhalf4).^2 )));
    
    uf = quarterstep.*(uip + k1./6 + k2./3 + k3./3) + k4./6;
    
    nf = nf + 16;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uf = ifft(uf); uc = ifft(uc); nf = nf +2;
    
    delta = sqrt(sum((abs(uf-uc)).^2))/sqrt(sum((abs(uf)).^2));
    
    if (delta < (tol/2))
        u1 = (16/15)*uf-(1/15)*uc;
        ufft = fft(u1); nf = nf +1;
        propagedlength = propagedlength + dz;
        dz = dz*(2^(1/5));
        
    elseif ( tol <= delta <= (2*tol))
        u1 = (16/15)*uf-(1/15)*uc;
        ufft = fft(u1); nf = nf +1;
        propagedlength = propagedlength + dz;
        dz = dz/(2^(1/5));
        
    elseif delta > (2*tol)
        dz = dz/2;
    else
        u1 = (16/15)*uf-(1/15)*uc;
        ufft = fft(u1); nf = nf +1;
        propagedlength = propagedlength + dz;
    end
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /L);
    
end

% giving output parameters
u1 = ifft(ufft);

