function [] = plot_model(config,input)
%PLOT_MODEL A function passed to Show operator to plot its input
%           config, the current configuration, TYPE. Configuration 
%           input, anything

%---------plot temporal input and output pulses
figure(1);
subplot(2,1,1)
plot (config.t,abs(input).^2,'r.-');axis tight;
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Pulse Shape');
%---------plot output spectrum------------------
spec = fftshift(fft(input));
specnorm = spec.*conj(spec);
specnorm = specnorm/max(specnorm);
subplot(2,1,2)
plot(simulation.Constants.c./(config.f + config.f0),specnorm,'r.-');axis tight;
grid on;
xlabel ('Wavelength (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');
pause(0.1)
end

