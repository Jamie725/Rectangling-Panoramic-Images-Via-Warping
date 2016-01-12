function [ outimage ] = drawGrid( gridmask, image )
%DRAWGRID Summary of this function goes here
%   Detailed explanation goes here
R= image(:,:,1);
G= image(:,:,2);
B= image(:,:,3);
R(gridmask) = 0;
G(gridmask) = 1;
B(gridmask) = 0;
outimage = cat(3,R,G,B);
end

