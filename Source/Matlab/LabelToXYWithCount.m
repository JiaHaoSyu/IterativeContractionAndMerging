function [MeanValxy,Num] = LabelToXYWithCount(data,LabelT,NumCount)
% syntax:
%       [MeanValxy Num] = LabelToXYWithCount(data,Label,NumCount)
% description:
%       Compute mean x,y of merging segments
%
% Input arguments:
    % data :  x and y with Mx2 matrix
    % LabelT :  merging information of regions, same number in the LabelT will merge together 
    % NumCount : number of pixels in each segments before merging
    
% Output arguments:
    % MeanValxy : Mean x,y after merging
    % Num : number of pixels in each segments after merging
    
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016
x = data(:,1);
y = data(:,2);

% Number of pixel of each segments
Num = sparse(LabelT(:),ones(size(LabelT(:))),NumCount);

% Summation of x,y of each segments
Sumx = sparse(LabelT(:),ones(size(LabelT(:))),x(:).*NumCount);
Sumy = sparse(LabelT(:),ones(size(LabelT(:))),y(:).*NumCount);
% Mean x,y of each segments
Meanx = Sumx./Num;
Meany = Sumy./Num;

MeanValxy = full(Meanx);
MeanValxy(:,2) = full(Meany);

