function [u1,nf] =  IP_CQEM_TD(u0,dt,L,dz,alpha,betap,gamma,fo,tol)

% This function solves the Generalized Nonlinear Schrodinger Equation with 
% thwe complete Raman response for pulse propagation in an optical fiber
% using the Interaction Picture Method in combination with the Conserved 
% Quantity Error method for the step-size determination and the time 
% domain integration of the nonlinear operator
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

while propagedlength < L,
    
    if (dz + propagedlength) > L,
        dz = L - propagedlength;
    end
    
    PhotonN = sum( (abs(ufft).^2)./(w' + 2*pi*fo) );
    
    halfstep = exp(linearoperator*dz/2);
    uip = ifft(halfstep.*ufft);
    convolution = ifft( hrw.*fft( abs(u1).^2 ) );
    k1 = -dz*1i*gamma*(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution );
    k1 = k1 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(u1*(1-fr).*abs(u1).^2 + dt*u1.*fr.*convolution );
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
    
    uaux = ifft(halfstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;    
    
    ufft = fft(uaux);
    nf = nf + 15;
    
    propagedlength = propagedlength + dz;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /L );
    
    % set dz for the next step
    error = abs(   sum( (abs(ufft).^2)./(w' + 2*pi*fo)) - PhotonN   ) / PhotonN;
    if error > 2 * tol,
        propagedlength = propagedlength - dz;
        dz = dz/2;
        ufft = fft(u1); nf = nf +1;
    else
        if error > tol,
            u1 = uaux;
            dz = dz/(2^0.2);
        else
            if error < 0.1*tol,
                u1 = uaux;
                dz = dz*(2^0.2);
            else
                u1 = uaux;
            end
        end
    end
    
    
end


