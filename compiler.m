% compile *.cpp to mex file

cd Source/C++

mex GetNeighborInformation.cpp
mex GetBorderDist.cpp
mex WLDFeature.cpp
mex TextureCascadeWithSlidingWindow.cpp
mex AccumulateRangeofTexture.cpp
mex SparseMatrix.cpp
mex SparseVector1.cpp
mex IndexToMatrix.cpp
mex CountNeiMerge.cpp
mex CountNeiMergeMatrix.cpp
mex LabelSum.cpp
mex MeanToLabel.cpp
mex FreqLabelOfSlidingWindow.cpp

cd ..
cd ..