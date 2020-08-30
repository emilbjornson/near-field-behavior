% This Matlab script generates Figure 6 in the paper:
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

%Propagation distance
d = 25;

%Number of antennas
N = logspace(0,10,100);

%Area of isotropic antenna
A = (lambda/4)^2;

%Computing free-space channel gain as defined in Eq. (1)
beta_d = A/(4*pi*d^2);


%Computing exact total channel gain as in Corollary 1
PrxPtx_exact = (N*beta_d)./(3*(N*beta_d*pi+1).*sqrt(2*N*beta_d*pi+1)) + 2/(3*pi)*atan(N*beta_d*pi./sqrt(2*N*beta_d*pi+1));


%Compute SNR scaling to get 0 dB SNR with N = 1
SNR0scaling = 1/min(PrxPtx_exact);


%Define the scaling exponents
rho = [0 1/2 1];


%Compute the SNRs with the different power scaling laws
SNR = zeros(length(PrxPtx_exact),length(rho));

for j = 1:length(rho)
    
    SNR(:,j) = (SNR0scaling./N.^rho(j)).*PrxPtx_exact;
    
end


%% Plot the simulation results
figure;
hold on; box on;
plot(N,SNR(:,1),'r','LineWidth',2);
plot(N,SNR(:,2),'b--','LineWidth',2);
plot(N,SNR(:,3),'k-.','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of receive antennas ($N$)','Interpreter','Latex');
ylabel('SNR','Interpreter','Latex');
legend({'$\rho=0$','$\rho=1/2$','$\rho=1$'},'Interpreter','Latex','Location','SouthWest');
set(gca,'fontsize',18);
