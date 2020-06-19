function [ NLOP ] = NonLinearoperator_w(u_t,gamma,w,fo,fr,hrw,dt,mod)
%   This function calculate the right half side of frequency domain GNLSE

    % Implement the right half side of equation (3) in the "Optimum 
    % Integration Procedures for Supercontinuum Simulation", which is a 
    % frequency domain formulation of equation (1)
    if ~isfield(mod,'ssp') || mod.ssp == 1
        NLOP = -1i*gamma*(1 + w/(2*pi*fo)).*fft(((1-fr)*u_t.*u_t.*conj(u_t))...
            + fr*dt*u_t.*ifft(hrw.*fft(u_t.*conj(u_t))));
    elseif fr ~= 0
        NLOP = -1i*gamma*fft(((1-fr)*u_t.*u_t.*conj(u_t))...
            + fr*dt*u_t.*ifft(hrw.*fft(u_t.*conj(u_t))));
    else
%         NLOP = -1i*gamma*fft((1-fr)*u_t.*abs(u_t).^2);
        NLOP = -1i*gamma*fft((1-fr)*u_t.*u_t.*conj(u_t));   % faster
    end
    % The factor dt in the second term is due to the discretization of the 
    % convolution between hr and u_t calculated in frequency domain using 
    % 2-FFT method
    % Please refer to section "Nonlinear fibre optics overview"of the book 
    % "SUPERCONTINUUM GENERATION IN OPTICAL FIBERS" 
    
end

