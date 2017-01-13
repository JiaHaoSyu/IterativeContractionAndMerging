function TypeInput(varargin)

if nargin==1
    error('Wrong number of input arguments')
elseif nargin==2
    error('Wrong number of input arguments')
elseif nargin==3
    error('Wrong number of input arguments')
elseif nargin==4
    File = varargin{1};
    if strcmp('Segs',varargin{2})
        error('Wrong number of input arguments')
    elseif strcmp('Hier',varargin{2})
        error('Wrong number of input arguments')
    elseif strcmp('Video',varargin{2})
        Type = 'Video';
        Show = '2';
        Para = varargin{3};
        FrameRate = varargin{4};
    end
elseif nargin==5
    File = varargin{1};
    if strcmp('Segs',varargin{2})
        Type = varargin{2};
        Segments = varargin{3};
        Show = varargin{4};
        Para = varargin{5};
    elseif strcmp('Hier',varargin{2})
        Type = varargin{2};
        HirSegs = varargin{3};
        Show = varargin{4};
        Para = varargin{5};
    elseif strcmp('Video',varargin{2})
        error('Wrong number of input arguments');
    end
else
    error('Wrong number of input arguments')
end

% Load RGBImage from input path
InputImg = im2double(imread([File.InputPath File.InputFilename]));

Para.RemnantNum = floor(size(InputImg,1)*size(InputImg,2)/5000);

% Weber's Local Descriptor(WLD) Texture Parameter
Para.Texture.T = 6;   % excitation
Para.Texture.D = 8;   % direction
Para.Texture.Window = 3; % window size of overlapping 3 for 7*7 window

% Para.SaveHir = 60;
% Para.Videoflag = 0; % Show Segmentation Result to Video
% Para.ShowAllorOff = -1; % 1 for All , 0 for SaveSeg, -1 for SaveHir

if strcmp(Type,'Segs')
%     [Para] = ICM_Segs(InputImg,Para,Segments);
    Para.Videoflag = 0; % no Video
    Para.SaveHir = 0; % no hierarchical
    Para.Segments = Segments; % specific number of segments
    [Para] = ICM(InputImg,Para);
    %     CostTime = toc(RGID)
    CostTime = Para.Phase1Time + Para.Phase2Time
    segs = Para.PickSeg;
    Display_Segs(segs,File,Show)
    save([File.OutputPath '/' File.OutputFilename(1:end-4) '.mat'],'segs') % for matlab
elseif strcmp(Type,'Hier')
    Para.Videoflag = 0; % no Video
    Para.SaveHir = HirSegs; % hierarchical image segmentation
    Para.Segments = 0; % no specific number of segments
    [Para] = ICM(InputImg,Para);
    CostTime = Para.Phase1Time + Para.Phase2Time
    
    ICPM = [];
    ICPM.Phase1Time = Para.Phase1Time;
    ICPM.Phase2Time = Para.Phase2Time;
    ICPM.Iter = Para.Phase2Iter;
    ICPM.Height = size(Para.OriImg,1);
    ICPM.Width = size(Para.OriImg,2);
    ICPM.P = Para.P;
    ICPM.Ps = Para.Ps;
      
    segs = {};
    len = size(Para.LabelCCombineTemp,1);
    if len >= Para.iter
        Count = 1;
        for j = Para.iter-1:-1:1
            LabelTC = Para.LabelCCombineTemp(j,:);
%             Label = LabelTC(Para.LabelIni);
            Label = reshape(LabelTC,[size(Para.OriImg,1) size(Para.OriImg,2)]);
            segs{1,Count} = Label;
            Count = Count + 1;
        end
    else
        Count = 1;
        for j = len:-1:1
            LabelTC = Para.LabelCCombineTemp(j,:);
            Label = reshape(LabelTC,[size(Para.OriImg,1) size(Para.OriImg,2)]);
%             Label = LabelTC(Para.LabelIni);
            segs{1,Count} = Label;
            Count = Count + 1;
        end
    end
    
    Display_Hier(segs,File,Show)
    
    save([File.OutputPath '/' File.OutputFilename(1:end-4) '.mat'],'segs') % for matlab
    mkdir([File.OutputPath '/ICPM/'])
    save([File.OutputPath '/ICPM/' File.OutputFilename(1:end-4) '_ICPM.mat'],'ICPM') % for matlab
elseif strcmp(Type,'Video')
    Para.Videoflag = 1; % video on
    Para.SaveHir = 0; % no hierarchical
    Para.Segments = 0; % no specific number of segments
    [Para] = ICM(InputImg,Para);
    CostTime = Para.Phase1Time + Para.Phase2Time
    
    Display_Video(Para,File,FrameRate)
       
    ICPM = [];
    ICPM.Phase1Time = Para.Phase1Time;
    ICPM.Phase2Time = Para.Phase2Time;
    ICPM.Iter = Para.Phase2Iter;
    ICPM.Height = size(Para.OriImg,1);
    ICPM.Width = size(Para.OriImg,2);
    ICPM.P = Para.P;
    ICPM.Ps = Para.Ps;
    mkdir([File.OutputPath '/ICPM/'])
    save([File.OutputPath '/ICPM/' File.OutputFilename(1:end-4) '_ICPM.mat'],'ICPM') % for matlab
    
    disp('done!!')
end

return

