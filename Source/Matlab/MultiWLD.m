function [F1,F2,F3,Excitation,Excitation2,Excitation3] = MultiWLD(Para)
% syntax:
%       [F1,F2,F3,Excitation,Excitation2,Excitation3] = MultiWLD(Para)
% description:
%       Compute multiscale texture feature
%
% Input arguments:
%       Para : Including T and D which indicates excitation and direction

% Output arguments:
%       F1 :  Texture feature of scale 1
%       F2 :  Texture feature of scale 2
%       F3 :  Texture feature of scale 3
%       Excitation : Excitation of scale 1
%       Excitation2 : Excitation of scale 2
%       Excitation2 : Excitation of scale 3
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016
% Original WLD code is available in http://www.cse.oulu.fi/CMV/Research/NewTextureDescriptors

image = Para.OriImg*255; % 0~1 => 0~255
grayimage = double(rgb2gray(uint8(image))); % Compute gray image

% parameters
T = Para.Texture.T; % excitation
D = Para.Texture.D; % direction
% Windows = Para.Texture.Window; % overlap window size 3 for (2*3+1)*(2*3+1) = 7*7

% Multi-WLD feature is computed within a local window 
% local window 3*3
R = 1;
f00 = [1, 1, 1; 1, -8, 1; 1, 1, 1];

[d_differential_excitation d_gradient_orientation] = WLDFeature(grayimage,f00,R,size(Para.OriImg,1),size(Para.OriImg,2));
d_differential_excitation = d_differential_excitation(1:end-R,1:end-R);
d_gradient_orientation = d_gradient_orientation(1:end-R,1:end-R);

% local window 5*5
R = 2;
f00 = [1, 1, 1 1 1;1 0 0 0 1;1 0 -16 0 1;1 0 0 0 1;1 1 1 1 1];
[d_differential_excitation2 d_gradient_orientation2] = WLDFeature(grayimage,f00,R,size(Para.OriImg,1),size(Para.OriImg,2));
d_differential_excitation2 = d_differential_excitation2(1:end-R,1:end-R);
d_gradient_orientation2 = d_gradient_orientation2(1:end-R,1:end-R);

% local window 7*7 
R = 3;
f00=[1 1 1 1 1 1 1;1 0 0 0 0 0 1;1 0 0 0 0 0 1;1 0 0 -24 0 0 1;1 0 0 0 0 0 1;1 0 0 0 0 0 1;1 1 1 1 1 1 1];

[d_differential_excitation3 d_gradient_orientation3] = WLDFeature(grayimage,f00,R,size(Para.OriImg,1),size(Para.OriImg,2));
d_differential_excitation3 = d_differential_excitation3(1:end-R,1:end-R);
d_gradient_orientation3 = d_gradient_orientation3(1:end-R,1:end-R);

Excitation = padarray(d_differential_excitation(2:end,2:end),[1 1],'replicate','both');
Orientation = padarray(d_gradient_orientation(2:end,2:end),[1 1],'replicate','both');

Excitation2 = padarray(d_differential_excitation2(3:end,3:end),[2 2],'replicate','both');
Orientation2 = padarray(d_gradient_orientation2(3:end,3:end),[2 2],'replicate','both');

Excitation3 = padarray(d_differential_excitation3(4:end,4:end),[3 3],'replicate','both');
Orientation3 = padarray(d_gradient_orientation3(4:end,4:end),[3 3],'replicate','both');

% Transform T and D to one dimension value T*D
[F1] = WLDToFeature(Excitation,Orientation,T,D); 
[F2] = WLDToFeature(Excitation2,Orientation2,T,D); 
[F3] = WLDToFeature(Excitation3,Orientation3,T,D); 
