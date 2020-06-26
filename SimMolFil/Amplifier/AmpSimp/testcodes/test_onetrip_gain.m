clear
% Active fiber parameters: 
gssdB = 30;             % (dB/m) 
PsatdBm = 60;           % (dBm) 
LossdB = 0;             % (dB)
loss = 10^(LossdB/10);
Length = 0.75;           	% (m)

d_z = 0.001;          	% (m)
N_z = Length/d_z;
Pin = 0.001;            % (W)

P_l2 = Pin;
gain2 = [];
tic
for ii = 1:N_z
    gain2 = [gain2,gain_saturated2(P_l2(ii),gssdB,PsatdBm)];
    P_l2 = [P_l2,exp((gain2(ii)-loss)*d_z)*P_l2(ii)];
end
gain2 = [gain2,gain_saturated2(P_l2(ii+1),gssdB,PsatdBm)];
toc

L = linspace(0,Length,N_z+1);
figure(1);
subplot(2,1,1),plot(L,P_l2);
subplot(2,1,2),plot(L,gain2);

% P_l1 = Pin;
% gain1 = [];
% GssdB = 10*log10(exp(10^(gssdB/10)*d_z));
% tic
% for ii = 1:N_z
%     gain1 = [gain1,gain_saturated3(P_l1(ii),GssdB,PsatdBm)];
%     P_l1 = [P_l1,gain1(ii)*exp(-loss*d_z)*P_l1(ii)];
% end
% gain1 = [gain1,gain_saturated3(P_l1(ii+1),GssdB,PsatdBm)];
% toc
% gain_1 = log(gain1)/d_z;
% 
% figure(2);
% subplot(2,1,1),plot(L,P_l1);
% % subplot(2,1,2),plot(L,gain1);
% subplot(2,1,2),plot(L,gain_1);