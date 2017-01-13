function ColorMap = RegionLabelToColor20131029(Label)
% syntax:
%       ColorMap = RegionLabelToColor20131029(Label)
% description:
%       Labeling image(segments) transform to color image
%
% Input arguments:
    % Label : Labeling Image(segments) (MxN)
    
% Output arguments:
    % ColorMap : Color Image  (MxNx3)
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

% random assign the color to each segments
R = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));
G = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));
B = sparse(ones(max(Label(:)),1),1:max(Label(:)),rand([max(Label(:)) 1]));

% Mapping to original image space with random color segments
ColorMap = full(R(Label));
ColorMap(:,:,2) = full(G(Label));
ColorMap(:,:,3) = full(B(Label));