function [u1,nf,Plotdata] =  SSFM_with_Raman(u0,dt,dz,nz,alpha,betap,gamma,fo,dplot)

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

hrw = Raman_response(t);

% constructing linear operator
linearoperator = -fftshift(alpha'/2);
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
% ufft = fft(u0);
fiberlength = nz*dz;
propagedlength = 0;
u1 = u0;
nf = 1;

if dplot ==1
    z_all = [];
    ufft_z = [];
    u_z = [];
end

while propagedlength < fiberlength,
    
    if (dz + propagedlength) > fiberlength,
        dz = fiberlength - propagedlength;
    end
    
    halfstep = exp(linearoperator*dz/2);
    uip = ifft(halfstep.*fft(u1));
    
    convolution = ifft( hrw.*fft( abs(u1).^2 ) );
    k1 = -dz*1i*gamma*(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution );
    k1 = k1 - (dz*gamma/(2*pi*fo))*(1/dt)*gradient(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution );
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
    
    u1 = ifft(halfstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;
    nf = nf + 16;
    
    ufft = fft(u1);
    
    propagedlength = propagedlength + dz;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /fiberlength );
    
    if dplot ==1
        z_all = [z_all;propagedlength];
        ufft_z = [ufft_z;fftshift(ufft)];
        u_z = [u_z;u1];
    end
    
end

% output parameter
% u1 = u1;

Plotdata.z = z_all;
Plotdata.ufft = abs(ufft_z);
Plotdata.u = abs(u_z);