% This Matlab script generates Figure 7 in the paper:
%
% Emil Bjornson, Luca Sanguinetti, “Power Scaling Laws and Near-Field
% Behaviors of Massive MIMO and Intelligent Reflecting Surfaces,” IEEE Open
% Journal of the Communications Society, to appear.
%
% Download article: https://arxiv.org/pdf/2002.04960
%
% This is version 1.0 (Last edited: 2020-08-29)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.


close all;
clear;

%Wavelength
lambda = 0.1;

%Propagation distances and angles for source and destination
d = 25;
delta = 2.5;
eta = pi/6;
omega = -pi/6;

%Number of IRS elements
NvaluesNonInteger = logspace(0,6,100);
Nsqrt = ceil(sqrt(NvaluesNonInteger));
Nvalues = Nsqrt.^2;

%Area of isotropic antenna
A = (lambda/4)^2;

%Side length of each element
a = sqrt(A);

%Computing the channel gain like expressions in Eq. (30)
varsigma_d_eta = cos(eta)*A/(4*pi*d^2);
varsigma_delta_omega = cos(omega)*A/(4*pi*delta^2);


%Compute locations of the source and destination
p_t = [d*sin(eta); 0; d*cos(eta)];
p_r = [delta*sin(omega); 0; delta*cos(omega)];


%Prepare to save simulation results
channelGain_IRS = zeros(length(Nvalues),1);
channelGain_IRS_farfield = zeros(length(Nvalues),1);
channelGain_IRS_UB = zeros(length(Nvalues),1);
channelGain_mMIMO = zeros(length(Nvalues),1);


%% Go through the different number elements
for j = 1:length(Nvalues)
    
    %Extract the number of elements/antennas
    N = Nvalues(j);
    
    %Prepare to store channel gains for individual elements/antennas
    betaHn = zeros(N,1);
    betaGn = zeros(N,1);
    
    %Go through each element/antenna and compute channel gains
    for n = 1:N
        
        %Compute location using Eqs. (22)-(23)
        x = -a*(sqrt(N)-1)/2 + a*mod(n-1,sqrt(N));
        y = a*(sqrt(N)-1)/2 - a*floor((n-1)/sqrt(N));
        
        %Compute channel gain for the n:th element
        betaHn(n) = channelgainGeneral(p_t,[x; y; 0],a);
        betaGn(n) = channelgainGeneral(p_r,[x; y; 0],a);
        
    end
    
    %Compute the exact total channel gain with the IRS using Eq. (42),
    %by removing the P/sigma^2 term
    channelGain_IRS(j) = sum(sqrt(betaHn.*betaGn)).^2;
    
    %Compute the far-field approximation of the total channel gain with the
    %IRS using Eq. (48), by removing the P/sigma^2 term
    channelGain_IRS_farfield(j) = N^2 * varsigma_d_eta * varsigma_delta_omega;
    
    %Compute the upper bound on the total channel gain with the IRS using
    %Eq. (43), by removing the P/sigma^2 term
    channelGain_IRS_UB(j) = channelGainArray(d,eta,N,A)*channelGainArray(delta,omega,N,A);
    
    %Compute the total channel gain with the mMIMO receiver using Eq. (26),
    %by removing the P/sigma^2 term
    channelGain_mMIMO(j) = channelGainArray(d,eta,N,A);

end


%% Plot the simulation results
figure;
hold on; box on;
plot(Nvalues,channelGain_mMIMO,'r-','LineWidth',2);
plot(Nvalues,channelGain_IRS_farfield,'k-.','LineWidth',2);
plot(Nvalues,channelGain_IRS_UB,'b:','LineWidth',2);
plot(Nvalues,channelGain_IRS,'b--','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of antennas/elements ($N$)','Interpreter','Latex');
ylabel('Total channel gain','Interpreter','Latex');
legend({'mMIMO','IRS (far-field approx)','IRS (upper bound)','IRS (exact)'},'Interpreter','Latex','Location','SouthEast');
set(gca,'fontsize',18);
