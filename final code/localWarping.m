function [dispMap, outputImg] = localWarping(img,mask)

	mex 'addSeam.cpp';
    img = uint8(img*255);
    mask = int32(mask);

	tstart = tic;
	[h, w, ~] = size(img);
%     thrs = 252;
% 	mask = int32(~mask_fg(img, thrs));
%     maskCheck = mask;

	%% Remove all -background boundary
% 	x = prod(mask,1) == 0;
% 	y = prod(mask,2) == 0;
% 	mask = mask(y,x);
% 	img = img(y,x,:);
% 	[h,w,~] = size(img);

	originImg = img;
    
	%% Local warping
	[xDispMap, yDispMap] = meshgrid(1:w, 1:h);
	xDispMap = int32(xDispMap);
	yDispMap = int32(yDispMap);
	while true %&& false
    	%tic;
    	[d, p0, p1] = findLB(mask);
    	fprintf('Insert seam at side %d from %d to %d\n', d, p0, p1);
    	switch d
        	case 1 % up
        	    [y, x] = meshgrid(h:(-1):1, p0:p1);
        	    axs = 1;
        	case 2 % down
        	    [y, x] = meshgrid(1:h, p1:(-1):p0);
        	    axs = 1;
        	case 3 % left
        	    [x, y] = meshgrid(w:(-1):1, p1:(-1):p0);
        	    axs = 0;
        	case 4 % right
        	    [x, y] = meshgrid(1:w, p0:p1);
        	    axs = 0;
        	otherwise
        	    break;
    	end
    	idx = sub2ind(size(mask), y, x);
    	idx2 = cat(3, idx, idx+h*w);
    	idx3 = cat(3, idx2, idx+2*h*w);
    	subImg = img(idx3);
    	subMask = mask(idx);
    	subXDispMap = xDispMap(idx);
    	subYDispMap = yDispMap(idx);
    	[outImg, outMask, outXDispMap, outYDispMap] = addSeam_wrap(subImg, subMask, subXDispMap, subYDispMap);
    	img(idx3) = outImg;
    	mask(idx) = outMask;
    	xDispMap(idx) = outXDispMap;
    	yDispMap(idx) = outYDispMap;
     	figure(1);imshow(img);
%     	figure; imshow(uint8(mask*255));
    	%toc;
	end

telapsed = toc(tstart);
fprintf('total %f sec used.\n', telapsed);

%% Try recovering image from displacement map
idx = sub2ind(size(mask), yDispMap, xDispMap);
idx3 = cat(3, idx, idx+h*w, idx+h*w*2);
imshow(originImg(idx3));

%Output dispMap
dispMap = zeros(h,w,2);
[xOrigMap, yOrigMap] = meshgrid(1:w, 1:h);

dispMap(:,:,1) = double(xDispMap) - xOrigMap;
dispMap(:,:,2) = double(yDispMap) - yOrigMap;

%Output image
outputImg = img;