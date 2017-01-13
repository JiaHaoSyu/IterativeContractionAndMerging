function Display_Video(Para,File,FrameRate)
% Display the type of video
disp('Write Video Results of Hierarchical Image Segmetation')
ImageR = tic;
InputImg = im2double(imread([File.InputPath File.InputFilename]));
% ==================== Write the result to video
[LabelImage] = MultiLabelToImage(Para.LabelCombine);
writerObj = VideoWriter([File.OutputPath '/' File.OutputFilename(1:end-4) '.avi']);
writerObj.Quality = 100;
writerObj.FrameRate = FrameRate;
open(writerObj);
for i = 1:Para.iter-1
    Img = LabelImage(:,:,:,i);
    Label = Para.LabelCombine(:,:,i);
    
    [dx dy] = gradient(Label);
    Mask = (abs(dx)+abs(dy))>0;
    ZC = zeros(size(Mask));
    Mask2 = cat(3,Mask,ZC);
    Mask2 = cat(3,Mask2,ZC);
    [MeanImg] = LabelToMeanImageRGB(InputImg,Label);
    MeanImgO = MeanImg;
    %         OursSeg = RegionLabelToColor20131029(Label);
    imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) 'MeanNoBoundaryLine.png'])
    MeanImgMask = repmat((1-Mask),[1 1 3]).*InputImg + repmat((Mask),[1 1 3]).*Mask2;
    MeanImg = repmat((1-Mask),[1 1 3]).*MeanImg + repmat((Mask),[1 1 3]).*Mask2;
    imwrite(MeanImg,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) 'Mean.png'])
    imwrite(Img,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) '_' num2str(max(Label(:))) 'Seg.png'])
    imwrite(MeanImgMask,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) 'Mask.png'])
%     imwrite(uint8(Label),[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) '.png'])
    
    CellContract = Para.XYCellCombine(:,:,:,i);
    imwrite(CellContract,[File.OutputPath '/' File.OutputFilename(1:end-4) '_Hier_' num2str(i) 'XY.png'])
   
    %         writeVideo(writerObj,[I_pre dataI1 CellContract;Img MeanImg MeanLabN])
    if (Para.iter-i)<=60
        for k = 1:5
            writeVideo(writerObj,[InputImg CellContract;Img MeanImgO]);
        end
    else
        writeVideo(writerObj,[InputImg CellContract;Img MeanImgO]);
    end
end
close(writerObj)
BuildVideoTime = toc(ImageR)