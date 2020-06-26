function [u1,nf,Plotdata] =  IP_CQEM_FD(u0,dt,dz,mod,fo,tol,dplot,quiet)
% This function solves the Generalized Nonlinear Schrodinger Equation with 
% the complete Raman response for pulse propagation in an optical fiber
% using the Interaction Picture Method in combination with the Conserved 
% Quantity Error method for the step-size determination and the frequency 
% domain integration of the nonlinear operator
%
% INPUT
%
% u0 - starting field amplitude (vector)
% dt - time step
% dz - initial step size
% mod - propogation module parameters field
%       mod.L - propagation distance
%       mod.alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
%       mod.betap - dispersion polynomial coefs, [beta_0 ... beta_m], or beta(w)
%               - the coefs are Tylor expansion at the central frequency fo
%       mod.gamma - nonlinearity coefficient
% fo -  central frequency of the simulation (THz)
% tol - relative photon error
% dplot - if set to 1, Plot data will be packed
%
% OUTPUT
%
% u1 - field at the output
% nf - number of FFTs performed

nt = length(u0);                            % number of sample points
w = fftshift(2*pi*(-nt/2:nt/2-1)/(dt*nt));  % angular frequencies (fftshifted)
t = (-nt/2:1:nt/2-1)*dt;                    % vector temporal (in ps)
df = 1/(nt*dt);
f=-(nt/2)*df:df:(nt/2-1)*df;                % frequencies vector (in THz)

% calculate the raman response fuction of frequency domain
[hrw,fr] = Raman_response_w(t,mod);

% preparation before the IPM
ufft = fft(u0);
propagedlength = 0;
u1 = u0;
nf = 1;

if isfield(mod,'gssdB') || isprop(mod,'gssdB')
    if isfield(mod,'gain_w')|| isprop(mod,'gain_w')
        gain_w = mod.gain_w;
    else
        gain_w = filter_lorentz_tf(u1,mod.fbw,mod.fc,fo,1/dt/nt); 
    end
    alpha_0 = mod.alpha;
end
LOP_0 = Linearoperator_w(mod.alpha,mod.betaw,w);

if dplot ==1
    dz_all = zeros(100,1);
    z_all = dz_all;
    ufft_z = zeros(100,nt);
    u_z = ufft_z;
end

% Performig the IPM ***********************************************
if ~quiet
    fprintf(1, '\nSimulation running...      ');
end

count = 0;
while propagedlength < mod.L
    
    if (dz + propagedlength) > mod.L
        dz = mod.L - propagedlength;
    end
    
    % (re)constructing linear operator
    if isfield(mod,'gssdB')|| isprop(mod,'gssdB')
        % modify alpha to include saturated gain and gain spectral
%         Pin0 = (sum(u1.*conj(u1))*nt);
        Pin0 = (sum(u1.*conj(u1))*dt)*mod.RepeatFre*1e-12;
%         Pin0 = (sum(u1.*conj(u1)))*mod.RepeatFre*1e-12;
        if isfield(mod,'pump_power')|| isprop(mod,'pump_power')
            if ~exist('pump_power_t','var')
                pump_power_t = mod.pump_power;
            end
            [gain,~,pump_power_t,~ ] = AmpGain( mod,pump_power_t,Pin0,dz*1e-3,'forward');
            if gain_w == 0
                gain_w = AmpGainWCal(mod, gain,dz);        % need normorlisation
                gain_s = interp1(f+fo,gain_w,mod.freq_s,'linear');
                gain_w(gain_w>0) = gain.*gain_w(gain_w>0)/abs(gain_s);  % gain section using amplifier model, with normorlisation
                gain_w(gain_w<=0) = -(1-10.^((mod.k0*mod.dopN*mod.Gamma_s.*gain_w(gain_w<=0))/10))*1e3; % loss section
                gain = gain_w;
            else
                gain = gain.*gain_w;
            end
        else
            gain = gain_saturated2(Pin0,mod.gssdB,mod.PsatdBm).*gain_w;
        end
        mod.alpha = alpha_0-gain;
        LOP = LOP_0 + fftshift(gain/2);             % fftshift
    elseif ~exist('LOP','var') || isempty(LOP)
        % In passive fiber, only calculate LOP once
        LOP = LOP_0;
    end
    
    % Calculate the real phonton number before and after dz propogation 
    % neglect the frequency dependant of S(w) of the equation (14) in 
    % "Efficient Adaptive Step Size Method for the Simulation of 
    % Supercontinuum Generation in Optical Fibers"
    % (Instead of calculating the partial differential of PN in (14), 
    % here we derectly calculate PN(z))
%     PhotonN = sum((abs(ufft).^2)./(w + 2*pi*fo));
    PhotonN_z = abs(sum(exp(-dz*fftshift(mod.alpha)).* ufft.*conj(ufft)./(w + 2*pi*fo)));

    % Applying the Runge-Kutta method in the interaction picture
    % Implement equations (5) in the "Optimum Integration Procedures for Supercontinuum Simulation"
    halfstep = exp(dz/2*LOP);
    uip = halfstep.*ufft;       % calculate A in interaction picture
    k1 = halfstep*dz.*NonLinearoperator_w(u1,mod.gamma,w,fo,fr,hrw,dt,mod);
    
    uhalf2 = ifft(uip + k1/2);
    k2 = dz*NonLinearoperator_w(uhalf2,mod.gamma,w,fo,fr,hrw,dt,mod);
    
    uhalf3 = ifft(uip + k2/2);
    k3 = dz*NonLinearoperator_w(uhalf3,mod.gamma,w,fo,fr,hrw,dt,mod);
    
    uhalf4 = ifft(halfstep.*(uip + k3));
    k4 = dz*NonLinearoperator_w(uhalf4,mod.gamma,w,fo,fr,hrw,dt,mod);
    
    uaux = halfstep.*(uip + k1./6 + k2./3 + k3./3) + k4./6;    
    
    propagedlength = propagedlength + dz;
    
    if ~quiet
        fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /mod.L );
    end
    
    % set dz for the next step
    error = abs(sum(uaux.*conj(uaux)./(w+2*pi*fo))-PhotonN_z)/PhotonN_z;
    if error > 2 * tol
        % error exceeds the double of tolerance, discard this calculation 
        % and roll back
        propagedlength = propagedlength - dz;
        dz = dz/2;  % reduce the step size by half 
    else
        % accept this step and optimise the next step size
        ufft = uaux;
        u1 = ifft(ufft);
        if error > tol
            % error exceeds tolerance, reduce step size
            dz = dz/(2^0.2);
        else
            % error is too small, increase step size to accelerate
            % change the factor from 0.1 to 0.5 to accelerate. Now both the 
            % speed and precision are mainly controled by tolerance
            if error < 0.5*tol
                dz = dz*(2^0.2);
            end
        end
        count = count +1;
        
        if dplot ==1       
            % save the plot of this step
            if count>size(dz_all,1)
                temp = zeros(100,1);            % buffered auto grow method to save time
                dz_all = [dz_all;temp];         %#ok<AGROW>
                z_all = [z_all;temp];           %#ok<AGROW>
            end
            if count>size(ufft_z,1)
                temp = zeros(200,nt);
                ufft_z = [ufft_z;temp];         %#ok<AGROW>
                u_z = [u_z;temp];               %#ok<AGROW>
            end
            dz_all(count) = dz;
            z_all(count) = propagedlength;
            ufft_z(count,:) = ufft;
            u_z(count,:) = u1;
        end
    end
    
    nf = nf + 16;
end

% pack the output struct
if dplot ==1
    Plotdata.dz = dz_all(1:count);
    Plotdata.z = z_all(1:count);
    Plotdata.ufft = abs(fftshift(ufft_z(1:count,:),2));
    Plotdata.u = u_z(1:count,:);
else
    Plotdata = [];
end
