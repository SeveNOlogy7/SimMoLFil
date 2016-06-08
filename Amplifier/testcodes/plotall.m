figure(79),surf(t,Plotdata.z,abs(Plotdata.u).^2);
colorbar;axis tight;shading interp;
ylabel ('x (km)');
xlabel ('t (ps)');
zlabel ('|u(z,t)|^2 (W)');
% axis([-0.5,0.5,0,max(Plotdata.z),0,max(max(abs(Plotdata.u).^2))])
view(0,90);
 
spec = abs(Plotdata.ufft').^2;
specnorm = spec ./(lambda'*ones(1,size(spec,2))).^2;
% specnorm = specnorm./(ones(size(spec,1),1)*max(specnorm));
specnorm = specnorm/max(specnorm(:));
% figure(80),waterfall(c./(f + fo),Plotdata.z,specnorm');
figure(80),surf(c./(f + fo),Plotdata.z,specnorm');
colorbar;axis tight;
shading interp;
ylabel ('x (km)');
xlabel ('lambda (nm)');
zlabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');
axis([1500,1600,0,max(Plotdata.z),0,1])
view(0,90);