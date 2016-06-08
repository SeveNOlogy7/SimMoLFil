function [Z, AT, AW, W] = gnlse(T, A, w0, gamma, betas,...
            loss, fr, RT, flength, nsaves)
% Propagate an optical field using the generalised NLSE
% This code integrates Eqs. (3.13), (3.16) and (3.17).
% For usage see the exampe of test_Dudley.m (below)
% Written by J.C. Travers, M.H. Frosz and J.M. Dudley (2009)
% Please cite this chapter in any publication using this code.
% Updates to this code are available at www.scgbook.info
n = length(T); dT = T(2)-T(1); % grid parameters
V = 2*pi*(-n/2:n/2-1)'/(n*dT); % frequency grid
alpha = log(10.^(loss/10)); % attenuation coefficient
B=0;
for i = 1:length(betas)     % Taylor expansion of betas
	B=B+ betas(i)/factorial(i+1).*V.^(i+1);
end
L = 1i*B - alpha/2;         % linear operator
if abs(w0) > eps            % if w0>0 then include shock
    gamma = gamma/w0;
    W=V+w0;                 %for shock W is true freq
else
    W = 1;                  % set W to 1 when no shock
end
RW = n*ifft(fftshift(RT.'));        % frequency domain Raman
L = fftshift(L); W = fftshift(W);   % shift to fft space
% === define function to return the RHS of Eq. (3.13)
function R = rhs(z, AW)
    AT = fft(AW.*exp(L*z));     % time domain field
    IT = abs(AT).^2;            % time domain intensity
    if (length(RT) == 1) || (abs(fr) < eps)     % no Raman case
        M = ifft(AT.*IT);       % response function
    else
        RS = dT*fr*fft(ifft(IT).*RW); % Raman convolution
        M = ifft(AT.*((1-fr).*IT + RS));% response function
    end
    R = 1i*gamma*W.*M.*exp(-L*z); % full RHS of Eq. (3.13)
end
% === define function to print ODE integrator status
function status = report(z, y, flag) %
    status = 0;
    if isempty(flag)
        fprintf('%05.1f %% complete\n', z/flength*100);
    end
end
% === setup and run the ODE integrator
Z = linspace(0, flength, nsaves); % select output z points
% === set error control options
options = odeset('RelTol', 1e-5, 'AbsTol',1e-12,...
    'NormControl','on',...
    'OutputFcn', @report);
[Z, AW] = ode45(@rhs, Z, ifft(A), options); % run integrator
% === process output of integrator
AT = zeros(size(AW(1,:)));
for i = 1:length(AW(:,1))
    AW(i,:) = AW(i,:).*exp(L.'*Z(i));       % change variables
    AT(i,:) = fft(AW(i,:));         % time domain output
    AW(i,:) = fftshift(AW(i,:))./dT;    % scale
end
W=V+w0; %the absolute frequency grid
end