function [u1,nf,Plotdata] =  SSFM_with_Raman(u0,dt,dz,mod,fo,dplot)
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
% meters, then mod.gamma should be specified in W^-1*m^-1.
% Similarly, if dt is given in picoseconds, and dz is given in
% meters, then beta(n) should have dimensions of ps^(n-1)/m.
nz = mod.L/dz;
nt = length(u0);                            % number of sample points
w = fftshift(2*pi*(-nt/2:nt/2-1)/(dt*nt));  % angular frequencies
t = (-nt/2:1:nt/2-1)*dt;                    % vector temporal (en ps)

[hrw,fr] = Raman_response_w(t,mod);

% Preparation before the SSFM
% ufft = fft(u0);
propagedlength = 0;
u1 = u0;
nf = 0;

if isfield(mod,'gssdB')
    gain_w = filter_lorentz_tf(u1,mod.fbw,mod.fc,fo,1/dt/nt); 
    alpha_0 = mod.alpha;
end

if dplot ==1
    z_all = [];
    ufft_z = [];
    u_z = [];
end

% Performig the SSFM ***********************************************
fprintf(1, '\nSimulation running...      ');

while propagedlength < mod.L,
    
    if (dz + propagedlength) > mod.L,
        dz = mod.L - propagedlength;
    end
    
    % (re)constructing linear operator
    if isfield(mod,'gssdB')
        % modify alpha to include gain
        Pin0 = (sum(u1.*conj(u1))/nt);
        gain = gain_saturated2(Pin0,mod.gssdB,mod.PsatdBm).*gain_w;
        mod.alpha = alpha_0-gain;
    end
    LOP = Linearoperator_w(mod.alpha,mod.betaw,w);
    
    halfstep = exp(LOP*dz/2);
    uip = ifft(halfstep.*fft(u1));
    
    convolution = ifft( hrw.*fft( abs(u1).^2 ) );
    k1 = -dz*1i*mod.gamma*(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution );
    k1 = k1 - (dz*mod.gamma/(2*pi*fo))*(1/dt)*gradient(u1*(1-fr).*abs(u1).^2 + dt*u1*fr.*convolution );
    k1 = ifft(halfstep.*fft(k1));
    
    uhalf2 = uip + k1/2;
    convolution = ifft( hrw.*fft( abs(uhalf2).^2 ) );
    k2 = -dz*1i*mod.gamma*(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution );
    k2 = k2 - (dz*mod.gamma/(2*pi*fo))*(1/dt)*gradient(uhalf2*(1-fr).*abs(uhalf2).^2 + dt*uhalf2*fr.*convolution ) ;
    
    uhalf3 = uip + k2/2;
    convolution = ifft( hrw.*fft( abs(uhalf3).^2 ) );
    k3 = -dz*1i*mod.gamma*(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    k3 = k3 - (dz*mod.gamma/(2*pi*fo))*(1/dt)*gradient(uhalf3*(1-fr).*abs(uhalf3).^2 + dt*uhalf3*fr.*convolution );
    
    uhalf4 = ifft(halfstep.*fft(uip + k3));
    convolution = ifft( hrw.*fft( abs(uhalf4).^2 ) );
    k4 = -dz*1i*mod.gamma*(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    k4 = k4 - (dz*mod.gamma/(2*pi*fo))*(1/dt)*gradient(uhalf4*(1-fr).*abs(uhalf4).^2 + dt*uhalf4*fr.*convolution );
    
    u1 = ifft(halfstep.*fft(uip + k1./6 + k2./3 + k3./3)) + k4./6;
    nf = nf + 16;
    
    ufft = fft(u1);
    
    propagedlength = propagedlength + dz;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /mod.L );
    
    if dplot ==1
        z_all = [z_all;propagedlength];
        ufft_z = [ufft_z;fftshift(ufft)];
        u_z = [u_z;u1];
    end
    
end

% pack the output field
Plotdata.z = z_all;
Plotdata.ufft = abs(ufft_z);
Plotdata.u = abs(u_z);