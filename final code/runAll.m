clc;
clear all;

nameset = {'1_input','2_input','5_input','6_input','7a_input','7b_input', ...
		'7c_input','7d_input','8a_input','9a_input','9b_input','9c_input',...
		'9d_input','10a_input','10b_input','10c_input'};

name = nameset{1} ;
origImg = double(imread(['sig13pano\' name '.jpg']))/255;
thrs = 253;
mask = int32(~mask_fg(uint8(origImg*255), thrs));

% Downsample
[rows,cols,color] = size(origImg);
megapixel = rows*cols;
scale = sqrt(1e6/megapixel);
origImg1M = imresize(origImg, scale,'bicubic');
mask1M = imresize(mask, scale,'bicubic');

[ dispMap, outImg ] = localWarping( origImg1M, mask1M );
[ Vlocal, Vglobal ] = globalmeshOpt( origImg1M,mask1M,dispMap,name );
outputimg = meshwarp(origImg1M,Vlocal, Vglobal);
%outimg = main(origImg,mask,nameset{1});
upScale = sqrt(megapixel/1e6);
%result = imresize(result,upScale,'bicubic');
finalImg = imresize(outputimg,upScale,'bicubic');