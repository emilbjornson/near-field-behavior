% This Matlab script generates Figure 11 in the paper:
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


close all
clear;

%Carrier frequency
f_c = 3e9;

%Wavelength
lambda = 3e8/f_c;

%Coordinates of point source
pt1 = [0,0,10];

%Values of x_n for the nth antenna element
xn_values = [0 5 10];


%Function to be integrated with right-hand side of (65), achievable by
%continuous matched filtering
continuous_MF_fun = @(x,y) (1/4/pi).*(pt1(3).*((pt1(1)-x).^2+pt1(3).^2))./(((pt1(1)-x).^2+(pt1(2)-y).^2+pt1(3).^2).^(5/2));

%Function to be integrated with left-hand side of (65), achievable by
%discrete matched filtering
discrete_MF_fun = @(x,y) sqrt((1/4/pi).*(pt1(3).*((pt1(1)-x).^2+pt1(3).^2))./(((pt1(1)-x).^2+(pt1(2)-y).^2+pt1(3).^2).^(5/2))).*exp(-1j.*2.*pi.*sqrt(((pt1(1)-x).^2+(pt1(2)-y).^2+pt1(3).^2))./lambda);


%Area of isotropic antenna
A = lambda^2/(4*pi);

%Width/length of the antenna element
width = lambda*logspace(-2, 1, 5000);

%Prepare to save results
continuous_MF_gain = zeros(length(width),length(xn_values));
discrete_MF_gain = zeros(length(width),length(xn_values));


for m = 1:length(xn_values)
    
    
    %Coordinates of nth antenna element
    xn = xn_values(m);
    yn = 0;
    
    for i = 1:length(width)
        
        %Extract antenna size
        a = width(i);
        
        %Channel gain with discrete MF
        discrete_MF_gain(i,m) = abs(integral2(discrete_MF_fun, xn-a/2,xn+a/2,yn-a/2,yn+a/2)).^2/(a*a);
        
        %Channel gain with continuos MF
        continuous_MF_gain(i,m) = integral2(continuous_MF_fun, xn-a/2,xn+a/2,yn-a/2,yn+a/2);
        
    end
    
end

%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(width/lambda, 10*log10(discrete_MF_gain(:,1)./continuous_MF_gain(:,1)),'k-', 'Linewidth', 2);
plot(width/lambda, 10*log10(discrete_MF_gain(:,2)./continuous_MF_gain(:,2)),'r-.', 'Linewidth', 2);
plot(width/lambda, 10*log10(discrete_MF_gain(:,3)./continuous_MF_gain(:,3)),'b:', 'Linewidth', 2);
set(gca,'XScale','log');
ylim([-50 0]);
legend({'$x_n=0$ m', '$x_n=5$ m','$x_n=10$ m'},'Interpreter','Latex','Location','NorthWest');
xlabel('Normalized antenna element size $a/\lambda$','Interpreter','Latex');
ylabel('Normalized channel gain [dB]','Interpreter','Latex');
set(gca,'fontsize',18);
