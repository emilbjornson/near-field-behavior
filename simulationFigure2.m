% This Matlab script generates Figure 2 in the paper:
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

%Computing the approximate total channel gain in the far-field
PrxPtx_approx = N*beta_d;

%Compute the maximum value that satisfies the rule-of-thumb
N_thumb = d^2/A/9;
PrxPtx_thumb = (N_thumb*beta_d)./(3*(N_thumb*beta_d*pi+1).*sqrt(2*N_thumb*beta_d*pi+1)) + 2/(3*pi)*atan(N_thumb*beta_d*pi./sqrt(2*N_thumb*beta_d*pi+1));


%% Plot the simulation results
figure;
hold on; box on;
plot(N,PrxPtx_approx,'b--','LineWidth',2);
plot(N,PrxPtx_exact,'r','LineWidth',2);
plot(N_thumb,PrxPtx_thumb,'ko','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of receive antennas ($N$)','Interpreter','Latex');
ylabel('Total channel gain','Interpreter','Latex');
legend({'Approximation','Exact'},'Interpreter','Latex','Location','NorthWest');
set(gca,'fontsize',18);
ylim([1e-6 1]);
