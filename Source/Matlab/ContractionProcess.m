function [TwistedLab,TwistedXY] = ContractionProcess(Lab,XY,LapMatrix,lambdaLab,lambdaXY)
% syntax:
%       [TwistedLab,TwistedXY] = ContractionProcess(Lab,XY,LapMatrix,lambdaLab,lambdaXY)
% description:
%       Pixel-based Contraction process
%
% Input arguments:
    % Lab : L*a*b* Image with I_Height*I_Width*3 matrix
    % XY : spatial coordinate with I_Height*I_Width*2 matrix
    % LapMatrix : laplacian matrix
    % lambdaLab : regularization paramters of Lab
    % lambdaXY : regularization paramters of XY
% Output arguments:
    % TwistedLab : Twisted color value with I_Height*I_Width*3 matrix
    % TwistedXY : Twiested spatial coordinate with I_Height*I_Width*2 matrix
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

[I_Height,I_Width,~] = size(Lab); 
% location and R G B value
x = XY(:,:,1);
y = XY(:,:,2);
L = Lab(:,:,1);
a = Lab(:,:,2);
b = Lab(:,:,3);

% I_Height by I_Width matrix ==> (I_Height*I_Width) by 1 matrix
x_bold = x(:);
y_bold = y(:);
L_bold = L(:);
a_bold = a(:);
b_bold = b(:);

% Contraction Process of XY
A = LapMatrix+ diag(sparse(ones(size(x(:)))*lambdaXY));
b = lambdaXY*[x_bold y_bold];

XY_CP = A\b;

% Contraction Process of Lab
A = LapMatrix+ diag(sparse(ones(size(x(:)))*lambdaLab));
b = lambdaLab*[L_bold(:) a_bold(:) b_bold(:)];

Lab_CP = A\b;

% (I_Height*I_Width) by 1 matrix ==> I_Height by I_Width matrix 
TwistedLab(:,:,1) = reshape(Lab_CP(:,1), [I_Height, I_Width]);
TwistedLab(:,:,2) = reshape(Lab_CP(:,2), [I_Height, I_Width]);
TwistedLab(:,:,3) = reshape(Lab_CP(:,3), [I_Height, I_Width]);

TwistedXY(:,:,1) = reshape(XY_CP(:,1), [I_Height, I_Width]);
TwistedXY(:,:,2) = reshape(XY_CP(:,2), [I_Height, I_Width]);
