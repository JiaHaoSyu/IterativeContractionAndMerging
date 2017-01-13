function [mid_inds,adj_inds,NeighborInf] = NeighborInformation(InputImg,Window)
% syntax:
%       [row_inds,col_inds,NeighborInf] = NeighborInformation(InputImg,Windows)
% description:
%       Sliding window to input image for Neighbor(adjacent) information
%       of pair pixels
%
% Input arguments:
%       InputImg : original input image, with MxNx3 matrix
%       Window : sliding window size, 1 for 3x3 window, 2 for 5x5 window         

% Output arguments:
%       mid_inds :  Index of middle pixel, 1x(total number of pair pixels) matrix
%       adj_inds :  Index of adjacent pixel, 1x(total number of pair pixels) matrix
%       NeighborInf :  [mid_inds' adj_inds'],(total number of pair pixels)x2 matrix
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016

[I_Height,I_Width,~] = size(InputImg);
OneFilter = ones(2*Window+1,2*Window+1); % WinFil are 3*3 all 1 value
AllOne = ones(I_Height,I_Width); % Image mask with all one value
ComputeTimes = conv2(AllOne,OneFilter,'same'); % Compute the computing times in each pixels, ex. (1,1) will compute the connection value to (1,1),(1,2),(2,1)and(2,2) Compute times = 4
TotalLength = int32(sum(ComputeTimes(:))); % Total Compute times of all pixels
[I_y,I_x] = meshgrid(1:I_Width,1:I_Height);
Index = sub2ind([I_Height I_Width],I_x,I_y); % Index value of each pixel
% Index = [1    I_Height+1  ...   (I_Width-1)*I_Height+1]
%         [2    I_Height+2  ...              .          ]
%           .        .           .           .          ]
%           .        .           .           .
%           .        .           .           .
%         [I_Height 2*I_Height ...   I_Width*I_Heigh    ]
% Compute the Middle index and Neighbor index for quickly computing
% D(i,j) of all connection with two vectors
%     [row_inds col_inds]
% Ex. [   1        1    ]
%     [   1        2    ]    pixel index 1 and neighbor pixel index 2
%     [   1   I_Height+1]
%     [   1   I_Height+2]
%     [   ...      ...  ]
[mid_inds,adj_inds] = GetNeighborInformation(ComputeTimes,Index,TotalLength,I_Height,I_Width,Window);

NeighborInf = [mid_inds' adj_inds'];
% Remove the same index (1,1),(2,2),....(N,N), because it didn't
% have connection with their own
MaskN = (NeighborInf(:,1)-NeighborInf(:,2))~=0;
NeighborInf = NeighborInf(MaskN,:);
mid_inds = NeighborInf(:,1)';
adj_inds = NeighborInf(:,2)';