function [Output] = LabelToMeanImageRGB(data,Label)
% syntax:
%       [Output] = LabelToMeanImageRGB(data,Label)
% description:
%       Compute mean image with Labeling image and original image
%
% Input arguments:
    % data :  original image (MxNx3)
    % Label : Labeling Image(segments) (MxN)
    
% Output arguments:
    % Output : Mean color Image  (MxNx3)
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016
R = data(:,:,1);
G = data(:,:,2);
B = data(:,:,3);

% Number of pixel of each segments
Num = LabelSum(Label(:),ones(size(Label(:),1),1),size(Label,1),size(Label,2),max(Label(:)))';
% Summation of R,G,B of each segments
SumR = LabelSum(Label(:),R(:),size(Label,1),size(Label,2),max(Label(:)))';
SumG = LabelSum(Label(:),G(:),size(Label,1),size(Label,2),max(Label(:)))';
SumB = LabelSum(Label(:),B(:),size(Label,1),size(Label,2),max(Label(:)))';

% Mean R,G,B of each segments
MeanR = SumR./Num;
MeanG = SumG./Num;
MeanB = SumB./Num;

% Mapping to original image space with Mean R,G,B
Output = MeanToLabel(Label,MeanR,size(Label,1),size(Label,2),max(Label(:)));
Output(:,:,2) = MeanToLabel(Label,MeanG,size(Label,1),size(Label,2),max(Label(:)));
Output(:,:,3) = MeanToLabel(Label,MeanB,size(Label,1),size(Label,2),max(Label(:)));

