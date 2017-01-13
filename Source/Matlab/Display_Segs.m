function Display_Segs(segs,File,Show)
% Display the type of Segments
if Show=='1'
    disp('Write Image Results of Image Segmetation')
    ImageR = tic;
    InputImg = im2double(imread([File.InputPath File.InputFilename]));
    Label = segs;
    Img = RegionLabelToColor20131029(Label);
    [dx dy] = gradient(Label);
    Mask = (abs(dx)+abs(dy))>0;
    ZC = zeros(size(Mask));
%     ZC = ones(size(Mask));
    Mask2 = cat(3,Mask,ZC);
    Mask2 = cat(3,Mask2,ZC);
    [MeanImg] = LabelToMeanImageRGB(InputImg,Label);
    MeanImgMask = repmat((1-Mask),[1 1 3]).*InputImg + repmat((Mask),[1 1 3]).*Mask2;
    MeanImg = repmat((1-Mask),[1 1 3]).*MeanImg + repmat((Mask),[1 1 3]).*Mask2;
    imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Segs_Mean.png'])
    imwrite(Img,[File.OutputPath '/' File.OutputFilename(1:end-4) '_SegsC.png'])
    imwrite(MeanImgMask,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Segs_Mask.png'])
    BuildImageTime = toc(ImageR)
    disp('done!!')
end
if Show=='2'
    disp('Write Labeling and Image Results of Image Segmetation')
    ImageR = tic;
    InputImg = im2double(imread([File.InputPath File.InputFilename]));
    Label = segs;
    Img = RegionLabelToColor20131029(Label);
    [dx dy] = gradient(Label);
    Mask = (abs(dx)+abs(dy))>0;
    ZC = zeros(size(Mask))
%     ZC = ones(size(Mask));
    Mask2 = cat(3,Mask,ZC);
    Mask2 = cat(3,Mask2,ZC);
    [MeanImg] = LabelToMeanImageRGB(InputImg,Label);
    imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) 'Segs_MeanNoBoundaryLine.png'])
    MeanImgMask = repmat((1-Mask),[1 1 3]).*InputImg + repmat((Mask),[1 1 3]).*Mask2;
    MeanImg = repmat((1-Mask),[1 1 3]).*MeanImg + repmat((Mask),[1 1 3]).*Mask2;
    imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Segs_Mean.png'])
    imwrite(Img,[File.OutputPath '/' File.OutputFilename(1:end-4) '_SegsC.png'])
    imwrite(MeanImgMask,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Segs_Mask.png'])
%     imwrite(uint8(Label),[File.OutputPath '/' File.OutputFilename(1:end-4) '_Segs.png'])
    BuildImageTime = toc(ImageR)
    disp('done!!')
end