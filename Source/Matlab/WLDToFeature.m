function [Mask3] = WLDToFeature(Excitation,Orientation,T,D)
% Transform T and D to one dimension value T*D

bins = -pi/2:pi/(T):pi/2;
Mask = zeros(size(Excitation));
count = 0;
for i = 1:length(bins)-1
    MM = Excitation >= bins(i) & Excitation < bins(i+1);
    Mask(MM) = count;
    count = count + 1;
end

bins = 0:360/D:360;
Mask2 = zeros(size(Orientation));
count = 0;
for i = 1:length(bins)-1
    MM = Orientation >= bins(i) & Orientation < bins(i+1);
    Mask2(MM) = count;
    count = count + 1;
end
Mask3 = Mask*D + Mask2;