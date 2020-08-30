% This Matlab script generates Figure 8 in the paper:
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

%Set transmit power divided by noise power
Ptxsigma2 = 0.01 / 1e-8;
Prelaysigma2 = 0.01 / 1e-8;

%Number of IRS elements
NvaluesNonInteger = logspace(0,6,100);
Nsqrt = ceil(sqrt(NvaluesNonInteger));
Nvalues = Nsqrt.^2;

%Area of isotropic antenna
A = (lambda/4)^2;

%Side length of each element
a = sqrt(A);

%Compute locations of the source and destination
p_t = [d*sin(eta); 0; d*cos(eta)];
p_r = [delta*sin(omega); 0; delta*cos(omega)];


%Prepare to save simulation results
SNR_IRS = zeros(length(Nvalues),1);
SNR_relaying = zeros(length(Nvalues),1);
SNR_mMIMO = zeros(length(Nvalues),1);


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
    
    %Compute the exact total channel gain with the IRS using Eq. (42)
    SNR_IRS(j) = Ptxsigma2 * sum(sqrt(betaHn.*betaGn)).^2;
    
    %Compute the total channel gain with the mMIMO receiver using Eq. (26)
    SNR_mMIMO(j) = Ptxsigma2*channelGainArray(d,eta,N,A);
    
    %Compute the total channel gain with the MIMO relay using Eq. (37)
    SNR_relaying(j) = min([SNR_mMIMO(j) Prelaysigma2*channelGainArray(delta,omega,N,A)]);
    
end


%Compute the information rates in the different setups
rate_IRS = log2(1+SNR_IRS);
rate_mMIMO = log2(1+SNR_mMIMO);
rate_relaying = (1/2)*log2(1+SNR_relaying);



%% Plot the simulation results
figure;
hold on; box on;
plot(rate_relaying,Nvalues,'b:','LineWidth',2);
plot(rate_IRS,Nvalues,'r-','LineWidth',2);
plot(rate_mMIMO,Nvalues,'k-.','LineWidth',2);
set(gca,'YScale','log');
xlabel('Information rate [bit/s/Hz]','Interpreter','Latex');
ylabel('Number of antennas/elements ($N$)','Interpreter','Latex');
legend({'Relaying','IRS','mMIMO'},'Interpreter','Latex','Location','SouthEast');
set(gca,'fontsize',18);
xlim([0 10]);
