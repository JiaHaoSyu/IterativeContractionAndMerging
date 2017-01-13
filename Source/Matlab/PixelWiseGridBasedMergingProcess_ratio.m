function [Label] = PixelWiseGridBasedMergingProcess_ratio(OperationImg,XY,RatioSpatial,BinColor,ColorI)
% Syntax:
%   [Label] = PixelWiseGridBasedMergingProcess_ratio(OperationImg,XY,RatioSpatial,BinColor,ColorI)
% Description:
%   Perform grid-based merging process with twisted color and spatial coordinate
%
% Input arguments:
%   OperationImg : original L,a,b image with MxNx3 matrix
%   XY : Twisted spatial coordinate with MxNx2 matrix
%   RatioSpatial : discretization invervals of x,y 
%   BinColor : discretization invervals of L,a,b
%   ColorI : Twisted L,a,b with MxNx3 matrix

% Output arguments:
%   Label : Output the segments after grid-based merging process

%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

[I_Height,I_Width,Dim] = size(OperationImg);

% Color information for grid-based method
OperationImgMin = zeros(Dim,1);
for k = 1:Dim
    OperationImgMin(k) = min(min(min(OperationImg(:,:,k))),min(min(ColorI(:,:,k))));
end

OperationImgMax = zeros(Dim,1);
for k = 1:Dim
    OperationImgMax(k) = max(max(max(OperationImg(:,:,k))),max(max(ColorI(:,:,k))));
end

ColorBinRange = (OperationImgMax-OperationImgMin)/ BinColor;

%% Grid-based method of (x,y,L,a,b)
BinSpatialRangeX = round(I_Width/RatioSpatial);
BinSpatialRangeY = round(I_Height/RatioSpatial);

I_Height = size( OperationImg, 1 );
I_Width = size( OperationImg, 2 );

ii = XY(:,:,1);
jj = XY(:,:,2);
downsampledWidth = ceil((I_Width)/BinSpatialRangeX +0.5) + 1;
downsampledHeight = ceil((I_Height)/BinSpatialRangeY +0.5) + 1;
downsampledDepth = ceil( (OperationImgMax-OperationImgMin)./ColorBinRange +0.5) + 1;

di = ( ii / BinSpatialRangeY ) + 0.5;
dj = ( jj / BinSpatialRangeX ) + 0.5;

dz1 = zeros(size(OperationImg));

for i = 1:Dim
    dz1(:,:,i) = (( ColorI(:,:,i)-OperationImgMin(i) )/ColorBinRange(i)) + 0.5;
end

gridSize = [downsampledHeight;downsampledWidth;downsampledDepth]';

gridNodeNumber = prod(gridSize);

ImageToGrid_Index = sub2ind(gridSize,ceil(di),ceil(dj),ceil(dz1(:,:,1)),ceil(dz1(:,:,2)),ceil(dz1(:,:,3)));
 
gridWeightsC = sparse(ImageToGrid_Index(:), 1, ones(numel(ImageToGrid_Index),1) , gridNodeNumber,1);

%% Non Empty OperationImg Construction
[NonEmptyGridIndex, ~, ~] = find(gridWeightsC);

NonEmptyNum = numel(NonEmptyGridIndex);

% Compute the mapping vector from cell index to Re-labeling index
grid2NonemptyIndex = SparseMatrix(NonEmptyGridIndex,ones(size(NonEmptyGridIndex,1),1),1:NonEmptyNum,size(NonEmptyGridIndex,1),gridNodeNumber,1);
%                             1110  2576  14304  32489
% grid2NonemptyIndex = [0 0 ... 4 ... 3 ... 1 ... 2 ...] mapping vector                                                         
% Re-labeling to the cell index ex. [14304 32489 2576 1110] =================> [1 2 3 4]
%                                      ImageToGrid_Index                         Label      
Label = SparseVector1(ImageToGrid_Index,grid2NonemptyIndex,size(ImageToGrid_Index,1),size(ImageToGrid_Index,2));





