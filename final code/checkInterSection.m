function [isInterSect, p]=checkInterSection(lineA, lineB)
%CHECKINTERSECTION Summary of this function goes here
%   Detailed explanation goes here
% A line is represented by starting and ending points
deltaA = lineA(2,:) - lineA(1,:);
deltaB = lineB(2,:) - lineB(1,:);
M = [ -deltaA' deltaB'];
t = (M)\(lineA(1,:) - lineB(1,:))';
isInterSect = false;
if sum(t<=1) == 2 && sum(t>=0) == 2
    isInterSect = true;
end
p = lineA(1,:)+t(1)*deltaA;

end

