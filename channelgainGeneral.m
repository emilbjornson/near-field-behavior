function zetaVal = channelgainGeneral(p_t,p_n,a)
% Computes the channel gain zeta_{p_t,p_n,a} in Eq. (4) when
% transmitting from an omni-directional antenna located at p_t to a
% receive antenna located at p_n in the XY-plane.
%
% This function was developed as a part of the paper:
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
%
% INPUT:
% p_t = Location (x_t,y_t,d) of the transmitter
% p_n = Location (x_n,y_n,0) of the receiver in the XY-plane
% a   = Area of the square-shaped receive antenna
%
% OUTPUT:
% zetaVal = Free-space channel gain


%Compute the X and Y sets defined in Lemma 1
X = [a/2+p_n(1)-p_t(1); a/2-p_n(1)+p_t(1)]';
Y = [a/2+p_n(2)-p_t(2); a/2-p_n(2)+p_t(2)]';

%Compute square distance in Z-dimension
d2 = (p_t(3) - p_n(3))^2;

%Compute the channel gain using Eq. (4)
zetaVal = 0;

for x = X
    
    x2 = x^2;
    
    for y = Y
        
        y2 = y^2;
        
        zetaVal = zetaVal + x*y/d2 / ( 3*(y2/d2+1)*sqrt(x2/d2+y2/d2+1)) + (2/3)*atan( x*y/d2 /  sqrt(x2/d2+y2/d2+1) );
        
    end
    
end

zetaVal = zetaVal/(4*pi);
