function [u1,Plotdata] = UPM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol,dplot)

% This function solves the nonlinear Schrodinger equation for
% pulse propagation in an optical fiber using the split-step
% Fourier method described in:
%
%  A. A. Rieznik,  T. Tolisano,  F. A. Callegari,  D. F. Grosz, and H. L. Fragnito,
% "Uncertainty relation for the optimization of optical-fiber transmission
% systems simulations," Opt. Express 13, 3822-3834 (2005).
% http://www.opticsexpress.org/abstract.cfm?URI=OPEX-13-10-3822
%
% USAGE
%
% u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma);
% u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma,tol);
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
% tol - convergence tolerance (default = 1e-5)
%
% OUTPUT
%
% u1 - field at the output
% number_of_FFTs - number of Fast Fourier Transforms performed during the
% propagation
%
% NOTES  The dimensions of the input and output quantities can
% be anything, as long as they are self consistent.  E.g., if
% |u|^2 has dimensions of Watts and dz has dimensions of
% meters, then gamma should be specified in W^-1*m^-1.
% Similarly, if dt is given in picoseconds, and dz is given in
% meters, then beta(n) should have dimensions of ps^(n-1)/m.

% if (nargin<8)
%   tol = 1e-3;
% end
nt = length(u0);
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);

% constructing linear operator
linearoperator = -fftshift(alpha'/2);
if (length(betap) == nt)     % If the user manually specifies beta(w)
    linearoperator = linearoperator - 1j*betap;
    linearoperator = fftshift(linearoperator);
else
    for ii = 0:length(betap)-1;
        linearoperator = linearoperator - 1j*betap(ii+1)*(w).^ii/factorial(ii);
    end
    linearoperator = conj(linearoperator');
%     linearoperator = fftshift(linearoperator);  %!!!!!
end

u1 = u0;
ufft = fft(u0);

fiberlength = nz*dz;
propagedlength =0;

if dplot ==1
    z_all = [];
    ufft_z = [];
    u_z = [];
end

% Performig the SSFM according to the UPM spatial-step size
fprintf(1, '\nSimulation running...      ');
while propagedlength < fiberlength,
    % calculating dz at each interaction
    Et = sum(abs(u1).^2);
    meanN = sum( gamma*(abs(u1).^4) ) / Et;
    aux = ( gamma * (abs(u1).^2) - meanN ).^2;
    deltaN = sqrt(sum( aux.*abs(u1).^2 ) / Et);
    
    Ef = sum(abs(ufft).^2);
    meanD = sum( 1j*linearoperator.*(abs(ufft).^2) ) / Ef;
    aux = ( 1j*linearoperator - meanD ).^2;
    deltaD = sqrt(sum( aux.*abs(ufft).^2 ) / Ef);
    
    
    dz =  (tol^(1/3)) * sqrt(1 / (deltaD*deltaN));
    % end of dz calculation at each inteaction
    
    if (dz + propagedlength) > fiberlength,
        dz = fiberlength - propagedlength;
    end
    
    % First linear parte (until dz/2)
    halfstep = exp(linearoperator*dz/2);
    uhalf = ifft(halfstep.*ufft);
    % non linear part
    u1 = uhalf .* exp(-j*gamma*(abs(uhalf).^2 )*dz);
    % Second linear part
    ufft = halfstep.*fft(u1);
    
    propagedlength = propagedlength + dz;
    fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
    
    if dplot ==1
        z_all = [z_all;propagedlength];
        ufft_z = [ufft_z;fftshift(ufft)];
        u_z = [u_z;u1];
    end
end

% giving output parameters
u1 = ifft(ufft);

Plotdata.z = z_all;
Plotdata.ufft = abs(ufft_z);
Plotdata.u = abs(u_z);

