main_SMF;
figure(1);
plot(t,abs(u).^2);axis tight;
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Pulse Shapes');

spec = fftshift(abs(fft(u)).^2);
specnorm = spec ./lambda.^2;
specnorm = specnorm/max(specnorm);

figure(3)
plot(c./(f + fo),specnorm);axis tight;
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');

load_simulator_data;
figure(1),hold on;
plot(t/1000,inten);

figure(3),hold on;
plot(c./f,specf);
% axis([1450,1650,0,1]);