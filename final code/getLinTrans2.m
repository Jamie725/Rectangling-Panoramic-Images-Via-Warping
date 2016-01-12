function [ T ] = getLinTrans2( V,p )
%GETLINTRANS Summary of this function goes here
%   Detailed explanation goes here
% V is a 2x4 vector  , p is 2x1 , T is 4x1 
% p = VT 
A = [eye(4) V'; V zeros(2)];
y = [zeros(4,1) ; p];
x = A\y;
T = x(1:4);


end

