function Display_Hier(segs,File,Show)
% Display the type of Hierarchical image segmentation
if Show=='1'
    disp('Write Image Results of Hierarchical Image Segmetation')
    ImageR = tic;
    InputImg = im2double(imread([File.InputPath File.InputFilename]));
    for j = 1:size(segs,2)
        Label = segs{j};
        Img = RegionLabelToColor20131029(Label);
        [dx dy] = gradient(Label);
        Mask = (abs(dx)+abs(dy))>0;
        ZC = zeros(size(Mask));
        Mask2 = cat(3,Mask,ZC);
        Mask2 = cat(3,Mask2,ZC);
        [MeanImg] = LabelToMeanImageRGB(InputImg,Label);
        MeanImgMask = repmat((1-Mask),[1 1 3]).*InputImg + repmat((Mask),[1 1 3]).*Mask2;
        MeanImg = repmat((1-Mask),[1 1 3]).*MeanImg + repmat((Mask),[1 1 3]).*Mask2;
        imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Mean.png'])
        imwrite(Img,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Seg.png'])
        imwrite(MeanImgMask,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Mask.png'])
    end
    BuildImageTime = toc(ImageR)
    disp('done!!')
end
if Show=='2'
    disp('Write Labeling and Image Results of Hierarchical Image Segmetation')
    ImageR = tic;
    InputImg = im2double(imread([File.InputPath File.InputFilename]));
    for j = 1:size(segs,2)
        Label = segs{j};
        Img = RegionLabelToColor20131029(Label);
        [dx dy] = gradient(Label);
        Mask = (abs(dx)+abs(dy))>0;
        ZC = zeros(size(Mask));
        Mask2 = cat(3,Mask,ZC);
        Mask2 = cat(3,Mask2,ZC);
        [MeanImg] = LabelToMeanImageRGB(InputImg,Label);
        imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'MeanNoBoundaryLine.png'])
        MeanImgMask = repmat((1-Mask),[1 1 3]).*InputImg + repmat((Mask),[1 1 3]).*Mask2;
        MeanImg = repmat((1-Mask),[1 1 3]).*MeanImg + repmat((Mask),[1 1 3]).*Mask2;
        imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Mean.png'])
        imwrite(Img,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Seg.png'])
        imwrite(MeanImgMask,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(j) 'Mask.png'])
    end
    BuildImageTime = toc(ImageR)
    disp('done!!')
end