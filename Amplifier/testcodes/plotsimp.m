% PLOT RESULTS ************************************************************
fprintf('\n\nPloting Results');

%---------plot temporal input and output pulses
figure(1);
plot (t,abs(u0).^2,'*-',t,abs(u).^2,'o-');axis tight;
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Initial (blue) and Final (green) Pulse Shapes');


%---------plot output spectrum------------------

spec = fftshift(abs(fft(u)).^2);
specnorm = spec ./lambda.^2;
specnorm = specnorm/max(specnorm);

figure(3)
plot(c./(f + fo),specnorm,'r.-');axis tight;
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');

fprintf(1,'\n\n----------------------------------------------\n');

PeakPower = max(abs(u).^2);
fprintf('Output Peak Power (W) = %5.2f\n',PeakPower);

b = dt*sum(abs(u).^2);
fprintf('Output Pulse Energy in pJ = %5.2f\n', b );