%  description:
%       The code is created based on the following paper
%       [1] Jia-Hao Syu, Sheng-Jyh Wang, and Li-Chun Wang, "Hierarchical Image Segmentation based 
%           on Iterative Contraction and Merging Process," IEEE Transaction on Image Process, 2017.
%
%       Code Author : Jia-Hao Syu <jiahaushiu4@gmail.com>
%       Data : 11/19/2016

clear all;
close all;
clc

% Add path of matlab and C++ mex function
addpath('Source/C++/')
addpath('Source/Matlab/')

%% Parameters
Para = [];

% adaptive adjusting affinity matrix for pixel-based and region-based
% Large value gives strong conhension force
Para.PercentPixelWise = 70; % pixel-based contraction-and-merging 70% affinity value larger than 0.01
Para.PercentRegionWise = 10; % region-based contraction-and-merging 10% affinity value larger than 0.01

% The use of features : size, textuer, border and spatial
% intertwining, 1 : turn on , 0 : turn off
Para.Property.RegionSizeOn = 1;  % equation (18)
Para.Property.TextureOn = 1; % equation (17) 
Para.Property.BorderAreaOn = 1; % equation (13)
Para.Property.SpatialIntertwiningOn = 1; % equation (19)

% Just noticeable difference(JND) in equation(10)
Para.phi = 0.4; %  phi

% Initial lambda value
Para.LambdaXY = 0.001;
Para.LambdaLab = 0.01;

% parameters of grid based merging process
Para.Grid.BinSpatial = 25; % quantizing bins of (x,y) 25
Para.Grid.BinColor = 15; % quantizing bins of (L,a,b) 15

% Learned paramters based on BSDS500
% size Parameter 
Para.t = 1.7; % size parameter t , ex. N^(1/t) in equation (18)
% weighting of similarity matrix in equation (20)
Para.alpha = 1; % color
Para.beta = 3; % weighting of texture
Para.gamma = 3; % weighting of border 

%% Three types of output
% 'Segs' : Superpixel-like Segmentation, user gives the number of segments, and our model provides the close segments 
% 'Hier' : Hierarchcial image segmentation
% 'Video' : Whole process result to Video, for 321*481 image, requiring 3GB Memory for saving all process result
%===== You can change the type =====
Type = 'Video'; % 'Segs','Hier','Video'

% Type 1 : Parameters for 'Segs'
NumSegment = 200; % Number of Segments (choosing by users)

% Type 2 : Parameters for 'Hier'
HierSegment = 60; % Number of Hierarchical image Segmentation, (1~HierSegment) segments 

% Display for Type 1 and 2 ('Segs','Hier')
Show = '2'; % '0','1', or '2', 0 for segs, 1 for segs and color image, 2 for all , default 1

% Type 3 : Parametrees for 'Video'
FrameRate = 40; 

for i = 1  
    
    DataPath = [pwd '/'];
    
    % Input image
    DataFilename = [DataPath num2str(i) '.jpg'];
    File = [];
    File.InputPath = DataPath;
    File.InputFilename = [num2str(i) '.jpg'];
    File.OutputPath = 'Output';
    File.OutputFilename = [num2str(i) '.jpg'];
    mkdir(File.OutputPath)
    
    switch Type
        case 'Segs'
            % ===========  image segmentation =================
            TypeInput(File,Type,NumSegment,Show,Para);
        case 'Hier'
            % ===========  hierarchical image segmentation======
            TypeInput(File,Type,HierSegment,Show,Para);
        case 'Video'
            % =========== Video of process result of hierarchical image segmentation==
            TypeInput(File,Type,Para,FrameRate);
    end
end
