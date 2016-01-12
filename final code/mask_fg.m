function mask = mask_fg(rgbImg, thrs)

[h, w, ~] = size(rgbImg);
grayImg = rgb2gray(rgbImg);

bdLoc = zeros(h+h+w+w, 2);
bdLoc(1:h, 1) = 1:h;
bdLoc(1:h, 2) = 1;
bdLoc(h+(1:h), 1) = 1:h;
bdLoc(h+(1:h), 2) = w;
bdLoc(h+h+(1:w), 1) = 1;
bdLoc(h+h+(1:w), 2) = 1:w;
bdLoc(h+h+w+(1:w), 1) = h;
bdLoc(h+h+w+(1:w), 2) = 1:w;

mask = grayImg < thrs;  %% dirty foreground
mask = imfill(mask, bdLoc, 4) & ~mask;  %% connected background
mask = imfill(mask, floor([h/2, w/2]), 4) & ~mask;  %% connected foreground
mask = imopen(mask, true(3));

%imshow(mask);
