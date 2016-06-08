function [u1,nf,Plotdata] =  IP_CQEM_FD(u0,dt,L,dz,alpha,betap,gamma,fo,tol,dplot)
% This function solves the Generalized Nonlinear Schrodinger Equation with 
% the complete Raman response for pulse propagation in an optical fiber
% using the Interaction Picture Method in combination with the Conserved 
% Quantity Error method for the step-size determination and the frequency 
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
%       - Note(CJH) the coefs are Tylor expansion at the central frequency fo
% gamma - nonlinearity coefficient
% fo -  central frequency of the simulation (THz)
% tol - relative photon error
%
% OUTPUT
%
% u1 - field at the output
% nf - number of FFTs performed

nt = length(u0);                            % number of sample points
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);  % angular frequencies (which is FFTshifted!!! CJH)
t = -(nt/2)*dt:dt:(nt/2-1)*dt;              % vector temporal (en ps)

hrw = Raman_response(t);

% constructing linear operator
linearoperator = -fftshift(alpha'/2);     % alpha 可被扩展为 alpha(w)--(CJH)
if (length(betap) == nt)        % If the user manually specifies beta(w)
    linearoperator = linearoperator - 1i*betap;
    linearoperator = fftshift(linearoperator);
else
    for ii = 0:length(betap)-1;
        linearoperator = linearoperator - 1i*betap(ii+1)*(w).^ii/factorial(ii);
    end
    linearoperator = conj(linearoperator');
%     linearoperator = fftshift(linearoperator);  % !!!!!!
end

% Performig the IPM ***********************************************
fprintf(1, '\nSimulation running...      ');
ufft = fft(u0);
propagedlength = 0;
u1 = u0;
nf = 1;

if dplot ==1
    z_all = [];
    ufft_z = [];
    u_z = [];
end

while propagedlength < L,
    
    if (dz + propagedlength) > L,
        dz = L - propagedlength;
    end
    
    PhotonN = sum( (abs(ufft).^2)./(w' + 2*pi*fo) );
    
    halfstep = exp(linearoperator*dz/2);
    uip = halfstep.*ufft;

    k1 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( u1.*((1-fr)*abs(u1).^2) + fr*dt*u1.*ifft(hrw.*fft( abs(u1).^2 )));
    k1 = halfstep.*k1;
    
    uhalf2 = ifft(uip + k1/2);
    k2 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf2.*((1-fr)*abs(uhalf2).^2) + fr*dt*uhalf2.*ifft(hrw.*fft( abs(uhalf2).^2 )));
    
    uhalf3 = ifft(uip + k2/2);
    k3 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf3.*((1-fr)*abs(uhalf3).^2) + fr*dt*uhalf3.*ifft(hrw.*fft( abs(uhalf3).^2 )));
    
    uhalf4 = ifft(halfstep.*(uip + k3));
    k4 = -dz*1i*gamma*(1 + w'/(2*pi*fo)).*fft( uhalf4.*((1-fr)*abs(uhalf4).^2) + fr*dt*uhalf4.*ifft(hrw.*fft( abs(uhalf4).^2 )));
    
    uaux = halfstep.*(uip + k1./6 + k2./3 + k3./3) + k4./6;    
    
    propagedlength = propagedlength + dz;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /L );
    
    % set dz for the next step
    error = abs(   sum( (abs(uaux).^2)./(w' + 2*pi*fo)) - PhotonN   ) / PhotonN;
    
    if error > 2 * tol,
        % discard this calculation and roll back
        propagedlength = propagedlength - dz;
        dz = dz/2;
    else
        ufft = uaux;
        u1 = ifft(ufft);
        if error > tol,
            dz = dz/(2^0.2);
        else
            if error < 0.1*tol,
                dz = dz*(2^0.2);
            end
        end
        if dplot ==1
            z_all = [z_all;propagedlength];
            ufft_z = [ufft_z;fftshift(ufft)];
            u_z = [u_z;u1];
        end
    end
    nf = nf + 16;
end
Plotdata.z = z_all;
Plotdata.ufft = abs(ufft_z);
Plotdata.u = abs(u_z);