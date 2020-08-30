% This Matlab script generates Figure 10 in the paper:
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
eta = 0;
omega = 0;

%Distance to the destination point
deltaFocus1 = 5;
deltaFocus2 = 25;
deltaRange = 1:0.1:100;

%Number of IRS elements
N = 1e4;

%Area of isotropic antenna
A = (lambda/4)^2;

%Side length of each element
a = sqrt(A);

%Compute locations of the source and destination
p_t = [d*sin(eta); 0; d*cos(eta)];
p_r1 = [deltaFocus1*sin(omega); 0; deltaFocus1*cos(omega)];
p_r2 = [deltaFocus2*sin(omega); 0; deltaFocus2*cos(omega)];

%Compute locations of the potential focus points for the destination
p_r_range = [deltaRange*sin(omega); zeros(size(deltaRange)); deltaRange*cos(omega)];


%Prepare to save simulation results
channelGain_IRS_focusing1 = zeros(length(deltaRange),1);
channelGain_IRS_focusing2 = zeros(length(deltaRange),1);
channelGain_IRS_optimized = zeros(length(deltaRange),1);
channelGain_IRS_mirror = zeros(length(deltaRange),1);



%Prepare to store channel gains for individual elements/antennas
betaHn = zeros(N,1);
betaGn = zeros(N,1);
betaGn_Range = zeros(N,1);
phaseHn = zeros(N,1);
phaseGn1 = zeros(N,1);
phaseGn2 = zeros(N,1);
phaseGn_Range = zeros(N,length(deltaRange));


%Go through each element/antenna and compute channel gains and
for n = 1:N
    
    %Compute location using Eqs. (22)-(23)
    x = -a*(sqrt(N)-1)/2 + a*mod(n-1,sqrt(N));
    y = a*(sqrt(N)-1)/2 - a*floor((n-1)/sqrt(N));
    
    %Compute channel gain for the n:th element
    betaHn(n) = channelgainGeneral(p_t,[x; y; 0],a);
    
    
    %Compute phase-shift for the n:th element to the destination
    phaseHn(n) = mod(norm(p_t-[x; y; 0])/lambda,1)*2*pi;
    phaseGn1(n) = mod(norm(p_r1-[x; y; 0])/lambda,1)*2*pi;
    phaseGn2(n) = mod(norm(p_r2-[x; y; 0])/lambda,1)*2*pi;
    
    %Compute phase-shift for the n:th element to different locations of the
    %destination
    for j = 1:length(deltaRange)
        
        betaGn_Range(n,j) = channelgainGeneral(p_r_range(:,j),[x; y; 0],a);
        phaseGn_Range(n,j) = mod(norm(p_r_range(:,j)-[x; y; 0])/lambda,1)*2*pi;
        
    end
    
end




for j = 1:length(deltaRange)
    
    %Compute the exact total channel gain with the IRS using Eq. (42) when
    %focusing on a destination at the wrong distance
    channelGain_IRS_focusing1(j) = abs(sum(sqrt(betaHn.*betaGn_Range(:,j)).*exp(-1i*(phaseGn_Range(:,j)-phaseGn1)))).^2;
    channelGain_IRS_focusing2(j) = abs(sum(sqrt(betaHn.*betaGn_Range(:,j)).*exp(-1i*(phaseGn_Range(:,j)-phaseGn2)))).^2;
    
    %Compute the exact total channel gain with the IRS using Eq. (42)
    channelGain_IRS_optimized(j) = sum(sqrt(betaHn.*betaGn_Range(:,j))).^2;
    
    %Compute the exact total channel gain with the IRS using Eq. (21) when
    %mimicking a mirror using theta_n=0 and mu_n=1
    channelGain_IRS_mirror(j) = abs(sum(sqrt(betaHn.*betaGn_Range(:,j)).*exp(-1i*(phaseHn+phaseGn_Range(:,j))))).^2;
    
end


%% Plot the simulation results
figure;
hold on; box on;
plot(deltaRange,channelGain_IRS_optimized,'b--','LineWidth',2);
plot(deltaRange,channelGain_IRS_focusing1,'k-','LineWidth',2);
plot(deltaRange,channelGain_IRS_focusing2,'r-.','LineWidth',2);
plot(deltaRange,channelGain_IRS_mirror,'b:','LineWidth',2);
set(gca,'YScale','log');
xlabel('Distance $\delta$ to the destination [m]','Interpreter','Latex');
ylabel('Total channel gain','Interpreter','Latex');
legend({'IRS (optimal)','IRS (focusing at 5 m)','IRS (focusing at 25 m)','IRS (mirror)'},'Interpreter','Latex','Location','NorthEast');
set(gca,'fontsize',18);
ylim([1e-9 1e-4]);
