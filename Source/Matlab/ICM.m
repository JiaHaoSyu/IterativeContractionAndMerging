function [Para] = ICM(InputImg,Para)
% syntax:
%       [Para] = ICM(InputImg,Para)
% description:
%       Perform iterative contraction-and-merging process with input image
%
% Input arguments:
%       InputImg : original input image (color image with MxNx3 matrix
%       gray image with MxN matrix)
%       Para : Initial parameters of our proposed method

% Output arguments:
%       Para : Hierarchical image segmentation and computational time
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016
%% Color Space Conversion : RGB to Lab
Para.OriImg = InputImg;
colorTransform = makecform('srgb2lab');
LabImg = applycform(InputImg, colorTransform);

[I_Height,I_Width,~] = size(InputImg);
LabImgR = reshape(LabImg,[I_Height*I_Width 3]); % MxNx3 maxtrix to (MxN)x3 matrix

%% Phase 1
Phase1 = tic;

% % Neighbor Information with a local sliding windows
[mid_inds,adj_inds,NeighborInf] = NeighborInformation(InputImg,1);

iter = 1; % first iterative contracton-and-merging process
Count = 1;

Ci = LabImgR(NeighborInf(:,1),:); % Middle pixel L*a*b* Color
Cj = LabImgR(NeighborInf(:,2),:); % Neighboring pixel L*a*b* Color
% compute color distance with middle pixels to their neighboring pixels
Dij = (sum(abs(Ci-Cj).^2,2).^(1/2));

DijT =  Dij';
% ====== Adaptive adjusting affinity matrix (reference the paper,how to derive tho)
ThoN = (DijT+1e-5)/(-log(0.01)); % eqation (8)
[ThoVal,~] = sort(ThoN,'descend');
% Compute the location that makes Para.PercentPixelWise percent of affinity value larger than 0.01
Leng = round((1-Para.PercentPixelWise/100)*size(DijT,2));
if Leng==0
    Leng = 1;
end
Tho = ThoVal(Leng); % Decided the Tho value
AffValue12 = exp(-(DijT+1e-5)/Tho); % eq(7)

% eq (10), if Neighboring regions with small distance < phi*JND Threshold, give large conhension force 10000
Mask = double(Dij<Para.phi*2.3);
AffValue12(Mask>0) = 10000;

% Transform neighbor information and affinity value to
% Affinity matrix and Compute Laplacian matrix
[LapMatrix] = ComputeLaplacianMatrix(size(LabImgR,1),mid_inds,adj_inds,AffValue12);
[y,x] = meshgrid(1:I_Width,1:I_Height);
XY = x;
XY(:,:,2) = y;

% Contraction process of xy and Lab , eq(6)
[TwistedLab,TwistedXY] = ContractionProcess(LabImg,XY,LapMatrix,Para.LambdaLab,Para.LambdaXY);

% Prevent the location XY after contraction process outside the image
TwistedXY(:,:,1) = min(max(floor(TwistedXY(:,:,1)+0.5),1),size(InputImg,1));
TwistedXY(:,:,2) = min(max(floor(TwistedXY(:,:,2)+0.5),1),size(InputImg,2));

% Pixel-Wise merging process
[Label] = PixelWiseGridBasedMergingProcess_ratio(LabImg,TwistedXY,Para.Grid.BinSpatial,Para.Grid.BinColor,TwistedLab);

LabelC = Label(:)';

%% update the neighboring information and feature of region information
LabelTC = LabelC;
NumSeg = max(LabelC(:))
Para.Ps = max(LabelC(:));
% Skewed lab(MeanValLab),skewed xy(XY) information from M-by-N-by-D to M*N-by-D
% Skewed : after contraction process
TwistedLab = reshape(TwistedLab,[I_Height*I_Width 3]);
TwistedXY = reshape(TwistedXY,[I_Height*I_Width 2]);

%  update the neighboring information with mapping vector LabelC
%  ex.  mapping vector = [50 60 60 50] for following case
%  Original row_indsM col_indsM  =>  new  row_indsM col_indsM
%              1          2    mapping       50        60
%              1          3                  50        60
%              1          4                  50        50
row_indsM = LabelTC(mid_inds);
col_indsM = LabelTC(adj_inds);
% Remove the repeated neighbor information ex.(50,60)
Index = sub2ind([max(row_indsM) max(col_indsM)],row_indsM,col_indsM);
[iia2,~,~] = unique(Index);
[row_indsMM,col_indsMM] = ind2sub([max(row_indsM) max(col_indsM)],iia2);
%  Remove the same index (1,1),(2,2),....(N,N), because connection didn't
% connect with their own, ex.(50,50)
NeighborInf = [row_indsMM' col_indsMM'];
MaskN = (NeighborInf(:,1)-NeighborInf(:,2))~=0;
NeighborInf = NeighborInf(MaskN,:);
row_indsMM = NeighborInf(:,1)';
col_indsMM = NeighborInf(:,2)';


% SizeRegion : number pixel in region, In pixel-based all value equals to 1
SizeRegion = ones(I_Height*I_Width,1); % Number pixel of region, in pixel-based contraction-and-merging all are 1

% update the (x,y) information with mapping vector, mapping vector
% knows whose pixels will merge togethor,
% TwistedXY : pixel of TwistedXY (x_s,y_s) merge to region with mean value of TwistedXY
[MeanValxy] = LabelToXYWithCount(TwistedXY,LabelC,SizeRegion);
MeanValxyM = MeanValxy;
SizeRegionT = SizeRegion; % for circle representation
% Update skewed the (L*,a*,b*) and SizeRegion(Output : not all 1, Input : all 1)
[MeanValLab,SizeRegion] = LabToMeanWithCount(TwistedLab,LabelC,SizeRegion);

% Circle Representation and only display in Type of video
% img for Showing the circle as initial background cyan color
img = zeros(I_Height,I_Width,3);
BackColor = [66 230 238]./255; % cyan color
img(:,:,1) = BackColor(1);
img(:,:,2) = BackColor(2);
img(:,:,3) = BackColor(3);

InputImgI = reshape(InputImg,[I_Height*I_Width 3]); % RGB or Gray Image for print to cell representation
[MeanValRGB] = LabToMeanWithCount(InputImgI,LabelC,SizeRegionT);
InputImgI = MeanValRGB;
if (Para.Videoflag==1)
    % Circle Parameter with (x,y,size) and we input the updated((x_s_r,y_s_r),(SizeRegion).^(1/4))
    circle = int32([MeanValxy(:,2) MeanValxy(:,1) (SizeRegion).^(1/4)]);
    % Matlab fuction for setting the Circle parameters
    shapeInserter = vision.ShapeInserter('Shape','Circles','Fill',true,'FillColor','Custom','CustomFillColor',MeanValRGB,'Opacity', 1);
    % Print the Circle to img(background with 0.5 value image)
    CellContract = step(shapeInserter,img,circle);
    LabelCombine(:,:,iter) = Label;
    XYCellCombine(:,:,:,iter) = CellContract;
else
    LabelCCombineTemp = zeros(Para.SaveHir,I_Height*I_Width);
    if NumSeg<=Para.SaveHir
       LabelCCombineTemp(Count,:) = LabelTC;
       Count = Count + 1;
    end
end

%% Remnant Removal
iter = iter + 1
while(sum(double(SizeRegion<=Para.RemnantNum))>0)
    
    % Find the region smaller than the Para.RemnantNum
    SmallCount = SizeRegion<=Para.RemnantNum;
    
    % Compute the Distance between Middle reginos and neighboring regions
    C1 = MeanValLab(NeighborInf(:,1),:);
    C2 = MeanValLab(NeighborInf(:,2),:);
    Dist12 = sum((C1-C2).^2,2).^(1/2);
    
    % Initial matrix for finding neighboring region with smallest distance
    Aff12Val = 9999*ones(size(SizeRegion,1),1);
    Aff12Ind = 9999*ones(size(SizeRegion,1),1);
    
    % Find neighboring region with smallest distance
    for j = 1:size(NeighborInf,1)
        if Aff12Val(NeighborInf(j,1))>Dist12(j,1)
            Aff12Ind(NeighborInf(j,1)) = NeighborInf(j,2);
            Aff12Val(NeighborInf(j,1)) = Dist12(j,1);
        end
    end
    
    % Assign small region(SmallCount) with strong affinity value 10000 to thier neighboring region with smallest distance(Aff12Ind)
    RandomAssign = sparse(1:size(SizeRegion,1),Aff12Ind,SmallCount.*10000,max(LabelC(:)),max(LabelC(:)));
    RandomAssign = RandomAssign + RandomAssign';
    % Matlab function of graphconncomp can find the connection of
    % graph ex. point1-point2, point2-point3, based on this
    % function, it will know point1-point2-piont3 connect together
    % and assign same labeling index(equal to mapping vector:LabelC)
    [~,LabelC] = conncomp(RandomAssign);
    
    % LabelInitial : original labeling index
    % LabelC : mapping vector
    % LabelTC : mapping vector for original labeling index, in
    % our hierarchical image segmentation process, i only save the LabelTC, and can direct
    % compute current labeling index(segments) with Label = LabelTC(LabelInitial)
    LabelTC = SparseVector1(LabelTC,LabelC,size(LabelTC,1),size(LabelTC,2)); %  LabelTC = LabelC(LabelTC);
   
    NumSeg = max(LabelC(:))
    
    % update the neighboring information with mapping vector LabelC
    % ex.  mapping vector = [50 60 60 50] for following case
    % Original row_indsM col_indsM   =>  new  row_indsM col_indsM
    %             1          2    mapping        50        60
    %             1          3                   50        60
    %             1          4                   50        50
    row_indsMT = LabelC(NeighborInf(:,1));
    col_indsMT = LabelC(NeighborInf(:,2));
    % Remove the repeated neighbor information
    Index = sub2ind([max(row_indsMT) max(col_indsMT)],row_indsMT,col_indsMT);
    [iia2,~,~] = unique(Index);
    [row_indsMM col_indsMM] = ind2sub([max(row_indsMT) max(col_indsMT)],iia2);
    %  Remove the same index (1,1),(2,2),....(N,N), because it didn't
    %  have connection with their own
    NeighborInf = [row_indsMM' col_indsMM'];
    MaskN = (NeighborInf(:,1)-NeighborInf(:,2))~=0;
    NeighborInf = NeighborInf(MaskN,:);
    row_indsMM = NeighborInf(:,1)';
    col_indsMM = NeighborInf(:,2)';
    
    % update the (x,y) information with mapping vector, mapping vector
    % knows whose pixels will merge togethor,
    % MeanValxy : skew region of (x_s,y_s) to region of (x_s_r,y_s_r)
    % SizeRegion : number pixel in region, (here not all 1)
    [MeanValxy] = LabelToXYWithCount(MeanValxy,LabelC,SizeRegion);
    [MeanValxyM] = LabelToXYWithCount(MeanValxyM,LabelC,SizeRegion);
    SizeRegionT = SizeRegion;
    % update the skewed (L*,a*,b*) and SizeRegion
    [MeanValLab,SizeRegion] = LabToMeanWithCount(MeanValLab,LabelC,SizeRegion);
    % update Mean RGB color of region for circle
    [MeanValRGB] = LabToMeanWithCount(InputImgI,LabelC,SizeRegionT);
    InputImgI = MeanValRGB;
    if (Para.Videoflag==1)
        % Circle Parameter with (x,y,size) and we input the updated((x_r,y_r),(SizeRegion).^(1/4))
        circle = int32([MeanValxy(:,2) MeanValxy(:,1) (SizeRegion).^(1/4)]);
        % Matlab fuction for setting the Circle parameters
        shapeInserter = vision.ShapeInserter('Shape','Circles','Fill',true,'FillColor','Custom','CustomFillColor',MeanValRGB,'Opacity', 1);
        % Print the Circle to img(background with 0.5 value image)
        CellContract = step(shapeInserter,img,circle);
        Label = reshape(LabelTC,[I_Height I_Width]);
        % save data for showing video
        % Mapping the skewed L*a*b* regions to original coordinate
        % Mapping the skewed x,y regions to original coordinate
        L = MeanValLab(:,1);L = L(Label);
        a = MeanValLab(:,2);a = a(Label);
        b = MeanValLab(:,3);b = b(Label);
        Lab(:,:,1) = L;Lab(:,:,2) = a;Lab(:,:,3) = b;
        X = MeanValxy(:,1);X = X(Label);
        Y = MeanValxy(:,2);Y = Y(Label);
        XYM(:,:,1) = X;XYM(:,:,2) = Y;
        LabelCombine(:,:,iter) = Label;
        XYCellCombine(:,:,:,iter) = CellContract;
    else
        if NumSeg<=Para.SaveHir
           LabelCCombineTemp(Count,:) = LabelTC;
           Count = Count + 1;
        end
    end
    iter = iter + 1
end

Label = reshape(LabelTC,[I_Height I_Width]);

Para.Phase1Time = toc(Phase1);
Para.P = max(LabelC(:));
%% Phase 2
Phase2 = tic;
%%  Border of Color eq.(13)
Windows =  1;
WinFil = ones(2*Windows+1,2*Windows+1);
SumLabel = conv2(Label,WinFil,'same');
AllOne = ones(size(LabImg,1),size(LabImg,2));
ComputeTimes = conv2(AllOne,WinFil,'same');
MulLabel = ComputeTimes.*Label;
Mask = double(abs(MulLabel-SumLabel)>0);
TotalLength = sum(sum(ComputeTimes.*Mask));
Length = max(Label(:));
% Only record boundary area of the middle region to neigboring region
[row_indsN,col_indsN,Dist1,BoundaryCount] = GetBorderDist(LabImg,Label,Mask,TotalLength,size(LabImg,1),size(LabImg,2),Windows,Length);
BoundaryCount = BoundaryCount - diag(diag(BoundaryCount));
BoundaryCountIni = BoundaryCount;
BoundaryCount = max(BoundaryCount,BoundaryCount');
% Computing total summation of boundary distance of middle region to neigboring region
SumBoundaryCount = SparseMatrix(row_indsN,col_indsN,Dist1+1e-5,length(row_indsN),max(Label(:)),max(Label(:)));
SumBoundaryCountIni = SumBoundaryCount;
SumBoundaryCount = max(SumBoundaryCount,SumBoundaryCount');
% Mean boundary distance = total summation of boundary distance/total times
MeanBoundary = zeros(size(BoundaryCount));
MeanBoundary(BoundaryCount>0) = SumBoundaryCount(BoundaryCount>0)./BoundaryCount(BoundaryCount>0);
% Transfer from boundary information matrix to vector by using neighboring information for similarity
IndN = sub2ind([size(MeanBoundary,1) size(MeanBoundary,2)],NeighborInf(:,1),NeighborInf(:,2));
Mean12 = full(MeanBoundary(IndN));

%% Spatial Intertwining eq.(19)
if  Para.Property.SpatialIntertwiningOn==1
    Para.Spatial.Window = 2; % spatial interwining 2 for 5*5 window
    % Each pixel cascaded the Label index with a sliding window 5*5
    [LabelCascade] = FreqLabelOfSlidingWindow(Label,Para.Spatial.Window,I_Height,I_Width);
    LabelMode = mode(LabelCascade,2); % Pick out most common index in each pixel
    LabelMode = reshape(LabelMode,I_Height,I_Width);
    % MRI after spatial intertwining of Labeling image
    MRI = Label;
    MRI(1+Para.Spatial.Window:end-Para.Spatial.Window,1+Para.Spatial.Window:end-Para.Spatial.Window) = LabelMode(1+Para.Spatial.Window:end-Para.Spatial.Window,1+Para.Spatial.Window:end-Para.Spatial.Window);
    % eq (20)
    deltaMRI = sparse(MRI(:),Label(:),ones(size(MRI,1)*size(MRI,2),1),max(Label(:)),max(Label(:)));
    % remove diagonal value
    deltaMRI = deltaMRI - diag(diag(deltaMRI));
    deltaMRIIni = full(deltaMRI);
    SI = min(deltaMRI,deltaMRI');
    % Transform matrix to [row_ind col_ind SI12]
    IndN = sub2ind([size(SI,1) size(SI,2)],NeighborInf(:,1),NeighborInf(:,2));
    SI12 = full(SI(IndN));
end

%% Texture eq.(17)
if Para.Property.TextureOn==1
    [F1,F2,F3,~,~,~] = MultiWLD(Para);
    % Accumulation range of texture information
    Windows = Para.Texture.Window; % overlap window size 3 for (2*3+1)*(2*3+1) = 7*7 window
    [I_j,I_i] = meshgrid(1:size(F1,2),1:size(F1,1));
    Index = sub2ind([size(F1,1) size(F1,2)],I_i,I_j);
    IndexT = Index(Windows+1:end-Windows,Windows+1:end-Windows);
    IndexT = padarray(IndexT,[Windows Windows],'replicate','both');
    K = 2*Para.Texture.Window + 1;
    T = Para.Texture.T;
    D = Para.Texture.D;
    % Cascade the texture information with sliding window and Repeat Label information for fast accumulation range of texture information in later operation
    [CascadeF1,CascadeF2,CasecadeF3] = TextureCascadeWithSlidingWindow(F1,F2,F3,Index,IndexT,Windows,size(Para.OriImg,1),size(Para.OriImg,2));
    RepeatLabel = repmat(Label(:),[K*K 1]);
    % Accumulation range of texture information
    W1 = AccumulateRangeofTexture(CascadeF1+1,RepeatLabel,K,size(Para.OriImg,1),size(Para.OriImg,2),max(Label(:)),T*D);
    W1Ini = W1;
    W2 = AccumulateRangeofTexture(CascadeF2+1,RepeatLabel,K,size(Para.OriImg,1),size(Para.OriImg,2),max(Label(:)),T*D);
    W2Ini = W2;
    W3 = AccumulateRangeofTexture(CasecadeF3+1,RepeatLabel,K,size(Para.OriImg,1),size(Para.OriImg,2),max(Label(:)),T*D);
    W3Ini = W3;
    % Normalized [W1 W2 W3]
    W1 = W1./repmat(sum(W1,2),[1 size(W1,2)]);
    W2 = W2./repmat(sum(W2,2),[1 size(W2,2)]);
    W3 = W3./repmat(sum(W3,2),[1 size(W3,2)]);
    W = [W1 W2 W3];
end

ThNc = 0.04*size(InputImg,1)*size(InputImg,2); % 4% of total pixels of N

Percent = Para.PercentRegionWise;
CountJ = 1; % for JND
NumSegOld = 10000;
while(NumSegOld>2)
    if NumSeg<=200
        Percent = 1;  % 1% affininty value larger than 0.01 in region-based contraction process and sequential merging
    end
    if(max(LabelC(:))<1)
        break;
    end
    NormalizeLabelNum = ones(1,max(LabelC(:))); % Compute the repeated connection times, due to the neighboring information remove the repeated connection, the normalized value will be all 1
    
    C1 = MeanValLab(NeighborInf(:,1),:); % Middle region Color (L*a*b*)
    C2 = MeanValLab(NeighborInf(:,2),:); % Neighboring region Color (L*a*b*)
    N1 = SizeRegion(NeighborInf(:,1)); % size of Middle region
    N2 = SizeRegion(NeighborInf(:,2)); % size of neighboring region
    N1 = N1.^(1/Para.t);
    N2 = N2.^(1/Para.t);
    % Size constraint, too larger size will dominant the result
    N1 = min(N1,ThNc.^(1/Para.t));
    N2 = min(N2,ThNc.^(1/Para.t));
    
    AB1 = MeanValLab(NeighborInf(:,1),2:3);  % Middle region (a*b*)
    AB2 = MeanValLab(NeighborInf(:,2),2:3);  % Neighboring region (a*b*)
    DistColor = (sum(abs(C1-C2).^2,2).^(1/2)); % color difference
    MeanDistColor = mean(DistColor);
    if Para.Property.TextureOn==1
        WR1 = W(NeighborInf(:,1),:);  % Middle region's WLD feature
        WR2 = W(NeighborInf(:,2),:);  % Neighboring region's WLD  feature
        DistW = sum(abs(WR1-WR2).^2,2).^(1/2); % WLD difference
        DistAB = sum(abs(AB1-AB2).^2,2).^(1/2); % (a*b*) difference
        %normalize DistW with same mean value of DistAB
        MeanDistAB = mean(DistAB)+1e-8; % mean of (a*b*) difference
        MeanDistW = mean(DistW)+1e-8;% mean of WLD difference
        DistW12N = MeanDistAB/MeanDistW*DistW; % normalize WLD feature with same scale(same difference mean) of (a*b*)
        DistT = DistAB.*DistW12N; % eq(17) % Combine (a*b*) and normalized WLD feature as Texture feature
        % normalize Texture with same scale(same difference mean) of color
        MeanDistT = mean(DistT)+1e-8;
        DistTN = MeanDistColor/MeanDistT*DistT;
    end
    % Border information eq.(13)
    if Para.Property.BorderAreaOn==1
        %  The boundary distance(mean boundary distance) between middle region and their neighboring region
        DistBorder = Mean12;
        % normalize Border with same scale(same difference mean) of color
        MeanDistBorderHat = mean(DistBorder)+1e-8;
        DistBorderHatN = MeanDistColor/MeanDistBorderHat*DistBorder;
    end
    % Similarity with Size,Texture,Border information
    if Para.Property.RegionSizeOn==1
        DistN = ((N1.*N2)./(N1+N2)); % eq(18)
        if Para.Property.TextureOn==1
            if Para.Property.BorderAreaOn==1
                Property = (Para.alpha*DistColor+Para.beta*DistTN+Para.gamma*DistBorderHatN);
            else
                Property = (Para.alpha*DistColor+Para.beta*DistTN);
            end
        else
            if Para.Property.BorderAreaOn==1
                Property = (Para.alpha*DistColor+Para.gamma*DistBorderHatN);
            else
                Property = (Para.alpha*DistColor);
            end
        end
        if Para.Property.SpatialIntertwiningOn == 1
            Temp = Property./((SI12+80).^(1/2));
            MeanTemp = mean(Temp)+1e-8;
            TempN = MeanDistColor/MeanTemp*Temp;
            Property = min(Property,TempN);
        end
        DistR12 = DistN.*Property; % eq (21)
    else
        if Para.Property.TextureOn==1
            if Para.Property.BorderAreaOn==1
                Property = (Para.alpha*DistColor+Para.beta*DistTN+Para.gamma*DistBorderHatN);
            else
                Property = (Para.alpha*DistColor+Para.beta*DistTN);
            end
        else
            if Para.Property.BorderAreaOn==1
                Property = (Para.alpha*DistColor+Para.gamma*DistBorderHatN);
            else
                Property = (Para.alpha*DistColor);
            end
        end
        if Para.Property.SpatialIntertwiningOn == 1
            Temp = Property./((SI12+80).^(1/2));
            MeanTemp = mean(Temp)+1e-8;
            TempN = MeanDistColor/MeanTemp*Temp;
            Property = min(Property,TempN);
        end
        DistR12 = Property;
    end
    DistR12T =  DistR12';
    
    %====== Adaptive adjusting affinity matrix(reference the paper,how to derive tho)
    ThoN = (DistR12T+1e-5)/(-log(0.01));
    [ThoVal,~] = sort(ThoN,'descend');
    % Compute the location that makes Para.PercentPixelWise percent of affinity value larger than 0.01
    Leng = round((1-Percent/100)*size(DistR12T,2));
    if Leng==0 % prevent zero index
        Leng = 1;
    end
    Tho = ThoVal(Leng); % Decided the Tho value
    AffValue12 = exp(-(DistR12T+1e-5)/Tho);  % eq (11)
    AffValue12(AffValue12>1) = 1;
    
    % eq (23), if Neighboring regions with small distance < k*JND Threshold, give large conhension force 10000
    JNDThreshold = min((Para.phi+0.1*CountJ)*2.3,2.3);
    Mask = double((MeanBoundary<JNDThreshold & BoundaryCount>0));
    Mask = full(Mask(IndN));
    AffValue12(Mask>0) = 10000;
    
    if NumSegOld==2
        AffValue12(AffValue12>0) = 1;
    end
    
    if Para.Segments>0
       Para.PickSeg = reshape(LabelTC,[I_Height I_Width]);
        if (max(Para.PickSeg(:))-Para.Segments)<=0;
            break;
        end
     end
    
    % Transform neighbor information and affinity value to
    % Affinity matrix and Compute Laplacian matrix
    [LapMatrixR] = ComputeLaplacianMatrix(size(MeanValLab(:,1),1),row_indsMM,col_indsMM,AffValue12);
    %================ Region-based Contraction process =================
    Likelihood_Precision = NormalizeLabelNum';
    MapI = MeanValxy(:,1);
    MapJ = MeanValxy(:,2);
    MapL = MeanValLab(:,1);
    Mapa = MeanValLab(:,2);
    Mapb = MeanValLab(:,3);
    
    A = LapMatrixR+ diag(sparse(Likelihood_Precision*Para.LambdaXY));
    b = Para.LambdaXY*[MapI(:) MapJ(:)];
    
    TwistedXY = A\b;
    
    A = LapMatrixR+ diag(sparse(Likelihood_Precision*Para.LambdaLab));
    b = Para.LambdaLab*[MapL(:) Mapa(:) Mapb(:)];
    
    TwistedLab = A\b;
    
    if NumSeg<=200
        % ============ Sequential merging ================
        XY1 = TwistedXY(NeighborInf(:,1),:); % skewed coordinate of middle region
        XY2 = TwistedXY(NeighborInf(:,2),:); % skewed coordinate of neighboring region
        Distance = DistR12T;  % Similarity matrix with region information(size,texture,color,.....)
        
        % normalized the skewed coordinate to the same scale of similarity matrix
        DistXY = sum(abs((XY1-XY2)).^2,2).^(1/2);
        DistXYN = mean(Distance)/mean(DistXY)*DistXY;
        Penalty = Distance' + DistXYN; % eq(25)
        
        % find the minimum distance, and merging
        [~,a] = sort(Penalty);
        AffValue12 = zeros(size(NeighborInf(:,1),1),1);
        AffValue12(a(1),1) = 1;
        AffValue12(a(2),1) = 1;
        % LabelC : mapping vector from previous Label to current Label
        [~,~,~,LabelC] = GraphAndConnectedComp(MeanValLab(:,1),NeighborInf(:,1),NeighborInf(:,2),AffValue12);
        % LabelTC : mapping vector from original labeling index to current Label, Label = LabelTC(LabelInitial)
        LabelTC = SparseVector1(LabelTC,LabelC,size(LabelTC,1),size(LabelTC,2)); %  LabelTC = LabelC(LabelTC);
    else
        % =================== Division method ==================
        [LabelT] = RegionWiseGridBasedMergingProcess_ratio(LabImg,TwistedXY,Para.Grid.BinSpatial,Para.Grid.BinColor,TwistedLab);
        % [LabelT] = RegionWiseGridBasedMergingProcess(OperationImg,XY,Para.Grid.BinSpatial,Para.Grid.BinColor,ColorI,NeighborInf);
        % [LabelT] = RegionWiseGridBasedMergingProcess_mod(OperationImg,XY,Para.Grid.BinSpatial,Para.Grid.BinColor,ColorI,NeighborInf,MeanValxyM);
        
        % LabelC : mapping vector from previous Label to current Label
        LabelC = LabelT';
        % LabelTC : mapping vector from original labeling index to current Label, Label = LabelTC(LabelInitial)
        LabelTC = SparseVector1(LabelTC,LabelC,size(LabelTC,1),size(LabelTC,2)); %  LabelTC = LabelC(LabelTC);
        
        %         Label = LabelTC(LabelInitial);
        %         figure(1),imshow(RegionLabelToColor20131029(Label));
    end
    
    % NumSegOld : Number of regions before merging
    % NumSeg : Number of region after merging
    NumSegOld = NumSeg;
    NumSeg = max(LabelC(:))
    
    %%  Update the neighboring information with mapping vector LabelC
    % ex.  mapping vector = [50 60 60 50] for following case
    % Original row_indsM col_indsM   =>  new  row_indsM col_indsM
    %             1          2    mapping        50        60
    %             1          3                   50        60
    %             1          4                   50        50
    row_indsMT = LabelC(NeighborInf(:,1));
    col_indsMT = LabelC(NeighborInf(:,2));
    % Remove the repeated neighbor information
    Index = sub2ind([max(row_indsMT) max(col_indsMT)],row_indsMT,col_indsMT);
    [iia2,~,~] = unique(Index);
    [row_indsMM,col_indsMM] = ind2sub([max(row_indsMT) max(col_indsMT)],iia2);
    
    [MeanValxy] = LabelToXYWithCount(TwistedXY,LabelC,SizeRegion);
    [MeanValxyM] = LabelToXYWithCount(MeanValxyM,LabelC,SizeRegion);
    SizeRegionT = SizeRegion;
    
    if NumSeg~=NumSegOld
        [MeanValLab,SizeRegion] = LabToMeanWithCount(MeanValLab,LabelC,SizeRegion);
    else
        [MeanValLab,SizeRegion] = LabToMeanWithCount(TwistedLab,LabelC,SizeRegion);
    end
    NeighborInf = [row_indsMM' col_indsMM'];
    MaskN = (NeighborInf(:,1)-NeighborInf(:,2))~=0;
    NeighborInf = NeighborInf(MaskN,:);
    row_indsMM = NeighborInf(:,1)';
    col_indsMM = NeighborInf(:,2)';
    % color of Border
    % update with mapping vector LabelC : direct summation the Border color information (do not recompute the Border information)
    BoundaryCount = CountNeiMerge(BoundaryCountIni,LabelC,size(BoundaryCountIni,1),size(BoundaryCountIni,2),max(LabelC(:)));
    BoundaryCountIni = BoundaryCount;
    BoundaryCount = max(BoundaryCount,BoundaryCount');
    SumBoundaryCount = CountNeiMerge(SumBoundaryCountIni,LabelC,size(SumBoundaryCountIni,1),size(SumBoundaryCountIni,2),max(LabelC(:)));
    SumBoundaryCountIni = SumBoundaryCount;
    SumBoundaryCount = max(SumBoundaryCount,SumBoundaryCount');
    MeanBoundary = zeros(size(BoundaryCount));
    MeanBoundary(BoundaryCount>0) = SumBoundaryCount(BoundaryCount>0)./BoundaryCount(BoundaryCount>0);
    IndN = sub2ind([size(MeanBoundary,1) size(MeanBoundary,2)],NeighborInf(:,1),NeighborInf(:,2));
    Mean12 = full(MeanBoundary(IndN));
    
    % spatial intertwining
    if  Para.Property.SpatialIntertwiningOn==1
        % update with mapping vector LabelC : direct summation the Spatialintertwining information (do not recompute the Spatialintertwining information)
        deltaMRI = CountNeiMerge(deltaMRIIni,LabelC,size(deltaMRIIni,1),size(deltaMRIIni,2),max(LabelC(:)));
        deltaMRIIni = full(deltaMRI);
        SI = min(deltaMRI,deltaMRI');
        % Transform matrix to [row_ind col_ind SI12]
        IndN = sub2ind([size(SI,1) size(SI,2)],NeighborInf(:,1),NeighborInf(:,2));
        SI12 = full(SI(IndN));
    end
    
    if Para.Property.TextureOn==1
        % update with mapping vector LabelC : direct summation the texture information (do not recompute the texture information)
        W1 = CountNeiMergeMatrix(W1Ini,LabelC,size(W1Ini,1),size(W1Ini,2),max(LabelC(:)),size(W1Ini,2));
        W2 = CountNeiMergeMatrix(W2Ini,LabelC,size(W2Ini,1),size(W2Ini,2),max(LabelC(:)),size(W1Ini,2));
        W3 = CountNeiMergeMatrix(W3Ini,LabelC,size(W3Ini,1),size(W3Ini,2),max(LabelC(:)),size(W1Ini,2));
        % Ini without normalized for next update
        W1Ini = W1;
        W2Ini = W2;
        W3Ini = W3;
        % normalized texture feature vector
        W1 = W1./repmat(sum(W1,2),[1 size(W1,2)]);
        W2 = W2./repmat(sum(W2,2),[1 size(W2,2)]);
        W3 = W3./repmat(sum(W3,2),[1 size(W3,2)]);
        W = [W1 W2 W3];
    end
    
    % Save for video
    if (Para.Videoflag==1)
        [MeanValRGB] = LabToMeanWithCount(InputImgI,LabelC,SizeRegionT);
        InputImgI = MeanValRGB;
        circle = int32([MeanValxy(:,2) MeanValxy(:,1) (SizeRegion).^(1/4)]);
        shapeInserter = vision.ShapeInserter('Shape','Circles','Fill',true,'FillColor','Custom','CustomFillColor',MeanValRGB,'Opacity', 1);
        CellContract = step(shapeInserter,img,circle);
        %             Label = LabelTC(LabelInitial);
        Label = reshape(LabelTC,[I_Height I_Width]);
        L = MeanValLab(:,1);L = L(Label);
        a = MeanValLab(:,2);a = a(Label);
        b = MeanValLab(:,3);b = b(Label);
        Lab(:,:,1) = L;Lab(:,:,2) = a;Lab(:,:,3) = b;
        X = MeanValxy(:,1);X = X(Label);
        Y = MeanValxy(:,2);Y = Y(Label);
        XYM(:,:,1) = X;XYM(:,:,2) = Y;
        LabelCombine(:,:,iter) = Label;
        XYCellCombine(:,:,:,iter) = CellContract;
    else
        if NumSeg<=Para.SaveHir
            LabelCCombineTemp(Count,:) = LabelTC;
            Count = Count + 1;
        end
    end
    iter = iter + 1
    CountJ = CountJ + 1;
end

if (Para.Videoflag==1)
    Para.LabelCombine = LabelCombine;
    Para.XYCellCombine = XYCellCombine;
else
    Para.LabelCCombineTemp = LabelCCombineTemp;
end

Para.iter = iter;
Para.Phase2Time = toc(Phase2);
Para.Phase2Iter = CountJ;



