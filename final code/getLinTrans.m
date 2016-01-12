function [ T ] = getLinTrans( V,p )
%GETLINTRANS Summary of this function goes here
%   Detailed explanation goes here
% V is a 8x1 vector  , p is 2x1 , T is 2x8

v1 = V(1:2);
v2 = V(3:4);
v3 = V(5:6);
v4 = V(7:8);

v21=v2 - v1;
v31=v3 - v1;
v41=v4 - v1;
p1=p - v1;
a1 = v31(1); % a is x
a2 = v21(1);
a3 = v41(1)- v21(1)-v31(1);
b1 = v31(2);
b2 = v21(2);
b3 = v41(2)- v21(2)-v31(2);
px=p1(1);
py=p1(2);

if a3 == 0 && b3 == 0
    tvec = inv([v31 v21])* p1;
    t1n = tvec(1);
    t2n = tvec(2);
else
    a =(b2*a3-a2*b3);
    b=(-a2*b1+b2*a1+px*b3-a3*py);
    c = px*b1 -py*a1;
    if a == 0
        t2n = -c/b;
    else
        t2n=(-b-sqrt(b^2-4*a*c))/(2*a);
    end
    if abs(a1+t2n*a3)<=0.00000001
        t1n= (py - t2n*b2)/(b1+t2n*b3);
    else
        t1n= (px - t2n*a2)/(a1+t2n*a3);
    end
end

m1 = v1 + t1n*(v3-v1);
m4 = v2 + t1n*(v4-v2);
ptest = m1 + t2n*(m4-m1);

v1w=1-t1n-t2n+t1n*t2n;
v2w=t2n-t1n*t2n;
v3w=t1n-t1n*t2n;
v4w=t1n*t2n;

T = [ v1w 0 v2w 0 v3w 0 v4w 0; ...
         0 v1w 0 v2w 0 v3w 0 v4w ];
assert(norm(T*V - p)<0.0001);
end

