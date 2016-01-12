function [outImg, outMask, outXDispMap, outYDispMap] = addSeam_wrap(img, mask, xDispMap, yDispMap);

img = permute(img, [3,2,1]);
mask = permute(mask, [2,1]);
xDispMap = permute(xDispMap, [2,1]);
yDispMap = permute(yDispMap, [2,1]);

[outImg, outMask, outXDispMap, outYDispMap] = addSeam(img, mask, xDispMap, yDispMap);

outImg = permute(outImg, [3,2,1]);
outMask = permute(outMask, [2,1]);
outXDispMap = permute(outXDispMap, [2,1]);
outYDispMap = permute(outYDispMap, [2,1]);
