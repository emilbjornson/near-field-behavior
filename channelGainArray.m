function etaVal = channelGainArray(d,eta,N,A)
% Computes the total channel gain eta_{d,eta,N} in Eq. (27) when
% transmitting from an omni-directional antenna to a planar square array,
% under the conditions defined in Assumptions 1 and 2. It can also be used
% for transmission from the array to an omni-directional receiver.
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
% d   = Distance from the transmitter to the center of the array
% eta = Angle in the XZ plane from transmitter to the array.
% N   = Number of square-elements in the array
% A   = Area of each square-element in the array
%
% OUTPUT:
% etaVal = Total free-space channel gain for the entire array


%Compute the variable B as defined in Proposition 1
B = N*A/(4*d^2*cos(eta)^2);


%Compute eta based on the summation expression in Eq. (27)
etaVal = 0;

for i = [1 2]
    
    etaVal = etaVal + (B+(-1)^i*sqrt(B)*tan(eta))/(6*pi*(B+1)*sqrt(2*B+tan(eta)^2+1+2*(-1)^i*sqrt(B)*tan(eta) )) + 1/(3*pi) * atan( (B+(-1)^i*sqrt(B)*tan(eta))/sqrt(2*B+tan(eta)^2+1+2*(-1)^i*sqrt(B)*tan(eta) ) );
    
end
