function [MeanValLab,Num] = LabToMeanWithCount(data,Label,NumCount)
% syntax:
%       [MeanValLab,Num] = LabToMeanWithCount(data,Label,NumCount)
% description:
%       Compute mean L,a,b of merging segments
%
% Input arguments:
    % data :  L,a,b with Mx3 matrix
    % LabelT :  merging information of regions, same number in the LabelT will merge together 
    % NumCount : number of pixels in each segments before merging
    
% Output arguments:
    % MeanValLab : Mean L,a,b after merging regions
    % Num : number of pixels in each segments after merging regions
    
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

L = data(:,1);
a = data(:,2);
b = data(:,3);

% Number of pixel of each segments
Num = full(sparse(Label(:),ones(size(Label(:))),NumCount));

% Summation of L,a,b of each segments
SumL = sparse(Label(:),ones(size(Label(:))),L(:).*NumCount);
Suma = sparse(Label(:),ones(size(Label(:))),a(:).*NumCount);
Sumb = sparse(Label(:),ones(size(Label(:))),b(:).*NumCount);

% Mean L,a,b of each segments
MeanL = SumL./Num;
Meana = Suma./Num;
Meanb = Sumb./Num;

MeanValLab = full(MeanL);
MeanValLab(:,2) = full(Meana);
MeanValLab(:,3) = full(Meanb);
