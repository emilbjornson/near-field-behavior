% This Matlab script generates Figure 3 in the paper:
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
N = logspace(4,12,100);

%Area of isotropic antenna
A = (lambda/4)^2;

%Computing free-space channel gain as defined in Eq. (1)
beta_d = A/(4*pi*d^2);

%Computing exact total channel gain as in Corollary 1
PrxPtx_exact = (N*beta_d)./(3*(N*beta_d*pi+1).*sqrt(2*N*beta_d*pi+1)) + 2/(3*pi)*atan(N*beta_d*pi./sqrt(2*N*beta_d*pi+1));


%Computing total channel gains by solving integrals that contain all or
%only some of the factors from (69).

fun1 = @(x,y) 1./(4*pi*(x.^2 + y.^2 + d^2));
fun2 = @(x,y) d./sqrt( x.^2 + y.^2 + d^2) .* 1./(4*pi*(x.^2 + y.^2 + d^2));
fun3 = @(x,y) d./sqrt( x.^2 + y.^2 + d^2) .* (x.^2 + d^2)./(x.^2 + y.^2 + d^2) .* 1./(4*pi*(x.^2 + y.^2 + d^2));

PrxPtx_exact1 = zeros(1,length(N));
PrxPtx_exact2 = zeros(1,length(N));
PrxPtx_exact3 = zeros(1,length(N));

for n = 1:length(N)
    
    PrxPtx_exact1(n) = integral2(fun1,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    PrxPtx_exact2(n) = integral2(fun2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    PrxPtx_exact3(n) = integral2(fun3,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    
end



%% Plot the simulation results
figure;
hold on; box on;
plot(N,PrxPtx_exact1,'b--','LineWidth',2);
plot(N,PrxPtx_exact2,'k-.','LineWidth',2);
plot(N,PrxPtx_exact,'r','LineWidth',2);
plot(N,PrxPtx_exact3,'r','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of receive antennas ($N$)','Interpreter','Latex');
ylabel('Total channel gain','Interpreter','Latex');
legend({'First property','First two properties','All three properties'},'Interpreter','Latex','Location','SouthEast');
set(gca,'fontsize',18);
ylim([1e-2 1e1]);
