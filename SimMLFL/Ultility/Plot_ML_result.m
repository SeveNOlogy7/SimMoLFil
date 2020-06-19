figure(4),pcolor(t,1:ii,abs(u_z(1:ii,:)).^2);
colorbar;axis tight;
shading flat;
ylabel ('run trip');
xlabel ('t (ps)');
zlabel ('|u(z,t)|^2 (W)');
title ('Output Optical Field Evolution');
view(0,90);

figure(5),pcolor(c./(f + fo),1:ii,spec_z(1:ii,:));
colorbar;axis tight;
shading flat;
ylabel ('run trip');
xlabel ('lambda (nm)');
zlabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum Evolution');
view(0,90);

repmat_length = round(size(Plotdata_smf5.u,1)/2);
ut_fft = repmat(abs(fftshift(fft(ut))),repmat_length,1);
u_f_fft = repmat(abs(fftshift(fft(u_f))),repmat_length,1);
uout_fft = repmat(abs(fftshift(fft(uout))),repmat_length,1);
u_nolm_fft = [Plotdata_smf1.ufft;Plotdata_smf2.ufft];
Plotdata.ufft = [u_f_fft;Plotdata_smf4.ufft;Plotdata_amf1.ufft;Plotdata_smf5.ufft;u_nolm_fft;ut_fft;uout_fft];

spec = abs(Plotdata.ufft').^2;
specnorm = spec/max(spec(:));
figure(7),pcolor(lambda,1:size(Plotdata.ufft,1),specnorm');
colorbar;axis tight;
shading flat;
xlabel ('lambda (nm)');
zlabel ('Normalized Spectrum (a.u.)');
title ('Last Trip Spectrum Evolution');
view(0,90);
ufft_length = [size(u_f_fft);size(Plotdata_smf4.ufft);size(Plotdata_amf1.ufft);...
    size(Plotdata_smf5.ufft);size(u_nolm_fft);size(ut_fft);size(uout_fft)];
text_all = {'Filter';'SMF1';'AF';'SMF2';'NOLM-SMF';'NOLM-out';'OC'};
current_height = 0;
for ii = 1:length(text_all)
    current_height = current_height+ufft_length(ii)/2;
    text(min(lambda),current_height,0,text_all(ii),'color','yellow','FontSize',10);
    current_height = current_height+ufft_length(ii)/2;
    line([min(lambda),max(lambda)],[current_height,current_height]+1,[0,0],...
        'LineWidth',0.5,'color','white');
end

ut_t = repmat(ut,repmat_length,1);
u_f_t = repmat(u_f,repmat_length,1);
uout_t = repmat(uout,repmat_length,1);
u_nolm_t = [Plotdata_smf1.u;Plotdata_smf2.u];
Plotdata.u = [u_f_t;Plotdata_smf4.u;Plotdata_amf1.u;Plotdata_smf5.u;u_nolm_t;ut_t;uout_t];
Plotdata.u = (abs(Plotdata.u)).^2;
figure(8),pcolor(t,1:size(Plotdata.u,1),Plotdata.u);
colorbar;axis tight;
shading flat;
xlabel ('t (ps)');
zlabel ('|u(z,t)|^2 (W)');
title ('Last Trip Pulse Evolution');
view(0,90);
current_height = 0;
for iii = 1:length(text_all)
    current_height = current_height+ufft_length(iii)/2;
    text(min(t),current_height,0,text_all(iii),'color','yellow','FontSize',10);
    current_height = current_height+ufft_length(iii)/2;
    line([min(t),max(t)],[current_height,current_height]+1,[0,0],...
        'LineWidth',0.5,'color','white');
end


NOLM_ufft = [u_nolm_fft(end,:);Plotdata_smf2_f.ufft(end,:);ut_fft(end,:)];
NOLM_ufft = abs(NOLM_ufft').^2;
NOLM_ufft = NOLM_ufft ./(lambda'*ones(1,size(NOLM_ufft,2))).^2;
NOLM_ufft = NOLM_ufft./max(NOLM_ufft(:));

figure(9);plot(lambda,NOLM_ufft);grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('NLM Specturm Interference');

NOLM_u = [u_nolm_t(end,:);Plotdata_smf2_f.u(end,:);ut(end,:)];
NOLM_u = abs(NOLM_u').^2;
phase_out = [diff(phase(u_nolm_t(end,:)));diff(phase(Plotdata_smf2_f.u(end,:)));diff(phase(ut(end,:)))];

figure(10);plot(t,NOLM_u);grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('NLM Pulse Interference');

AC = xcorr(Eout,Eout);
AC = AC/max(AC);
figure(10),plot(-time+dt:dt:time-dt,AC);
xlabel ('Time Delay (ps)');
ylabel ('AC Trace(a.u.)');
