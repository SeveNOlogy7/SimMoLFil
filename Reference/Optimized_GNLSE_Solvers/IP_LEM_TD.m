function [u1,nf] = IP_LEM_TD(u0,dt,L,dz,alpha,betap,gamma,fo,tol)

% This function solves the Generalized Nonlinear Schrodinger Equation with 
% thwe complete Raman response for pulse propagation in an optical fiber
% using the Interaction Picture Method in combination with the Local Error 
% method for the step-size determination and the time domain
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
    uhalf = ifft(halfstep.*ufft); nf = nf +1;
    uquarter=ifft(quarterstep.*ufft); nf = nf +1;
    
    
    %%%%% NON LINEAR OPERATOR COARSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uip = uhalf;
    
    convolution1 = ifft( hrw.*fft( abs(u1).^2 ) );
    k1 = -dz*1i*gamma*(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution1 );
    k1 = k1 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(u1*(1-fr).*abs(u1).^2 + dt*u1.*fr.*convolution1 );
    k1 = ifft(halfstep.*fft(k1));
    
    uhalf2 = uip + k1/2;
    convolution = ifft( hrw.*fft( abs(uhalf2).^2 ) );
    k2 = -dz*1i*gamma*(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution );
    k2 = k2 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution ) ;
    
    uhalf3 = uip + k2/2;
    convolution = ifft( hrw.*fft( abs(uhalf3).^2 ) );
    k3 = -dz*1i*gamma*(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    k3 = k3 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    
    uhalf4 = ifft(halfstep.*fft(uip + k3));
    convolution = ifft( hrw.*fft( abs(uhalf4).^2 ) );
    k4 = -dz*1i*gamma*(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    k4 = k4 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    
    uc = ifft(halfstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;
    
    nf = nf + 14;
    %%%%%%%%% END COARSE%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%FIRST NON LINEAR FINE%%%%%%%%%%%%%%%%%%%%%%%%
    uip = uquarter;
    
    k1 = -0.5*dz*1i*gamma*(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution1 );
    k1 = k1 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(u1*(1-fr).*abs(u1).^2 + dt*u1.*fr.*convolution1 );
    k1 = ifft(quarterstep.*fft(k1));
    
    uhalf2 = uip + k1/2;
    convolution = ifft( hrw.*fft( abs(uhalf2).^2 ) );
    k2 = -0.5*dz*1i*gamma*(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution );
    k2 = k2 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution ) ;
    
    uhalf3 = uip + k2/2;
    convolution = ifft( hrw.*fft( abs(uhalf3).^2 ) );
    k3 = -0.5*dz*1i*gamma*(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    k3 = k3 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    
    uhalf4 = ifft(quarterstep.*fft(uip + k3));
    convolution = ifft( hrw.*fft( abs(uhalf4).^2 ) );
    k4 = -0.5*dz*1i*gamma*(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    k4 = k4 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    
    uf = ifft(quarterstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;
    ufftf = fft(uf);

    nf = nf + 13;
    
    %%%%%%%SECOND NON LINEAR FINE----------------
    uip = ifft(quarterstep.*ufftf);
    
    convolution = ifft( hrw.*fft( abs(uf).^2 ) );
    k1 = -0.5*dz*1i*gamma*(uf*(1-fr).*abs(uf).^2 + dt*uf*fr.*convolution );
    k1 = k1 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uf*(1-fr).*abs(uf).^2 + dt*uf.*fr.*convolution );
    k1 = ifft(quarterstep.*fft(k1));
    
    uhalf2 = uip + k1/2;
    convolution = ifft( hrw.*fft( abs(uhalf2).^2 ) );
    k2 = -0.5*dz*1i*gamma*(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution );
    k2 = k2 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution ) ;
    
    uhalf3 = uip + k2/2;
    convolution = ifft( hrw.*fft( abs(uhalf3).^2 ) );
    k3 = -0.5*dz*1i*gamma*(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    k3 = k3 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    
    uhalf4 = ifft(quarterstep.*fft(uip + k3));
    convolution = ifft( hrw.*fft( abs(uhalf4).^2 ) );
    k4 = -0.5*dz*1i*gamma*(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    k4 = k4 - (0.5*dz*gamma/(2*pi*fo))*(1/dt)*gradient(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    
    uf = ifft(quarterstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;
    
    nf = nf + 15;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

