Inten_i = abs(u0).^2;
Inten_t = abs(u).^2;
Pin_min = 2.5e0;
[Pin_max,I_Pin_max] = max(Inten_i);
I_Pin_min = find(Inten_i>Pin_min,1)

trans = Inten_t(I_Pin_min:I_Pin_max)./Inten_i(I_Pin_min:I_Pin_max);
Pin_axis = Inten_i(I_Pin_min:I_Pin_max);

figure(10),plot(Pin_axis,trans),hold on;