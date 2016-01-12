v1 = [1 5];
v2 = [3 23];
v3 = [37 1];
v4 = [29 31];
V=[v1; v2; v3; v4];
p = (v4+v3)/2;

syms t1 t2

eqn = p== v1 +t1*(v3-v1) + t2*(v2-v1)+t2*t1*(v4-v2+v1-v3);

[t1sol, t2sol] = vpasolve(eqn);
vpa(t2sol)
vpa(t1sol)
%% 
Vv = [v1 v2 v3 v4]';
pv =p';
% v21=v2 - v1;
% v31=v3 - v1;
% v41=v4 - v1;
% p1=p - v1;
% a1 = v31(2);
% a2 = v21(2);
% a3 = v41(2)- v21(2)-v31(2);
% b1 = v31(1);
% b2 = v21(1);
% b3 = v41(1)- v21(1)-v31(1);
% px=p1(2);
% py=p1(1);
% a =(b2*a3-a2*b3);
% b=(-a2*b1+b2*a1+px*b3-a3*py);
% c = px*b1 -py*a1;
% t2n=(-b-sqrt(b^2-4*a*c))/(2*a);
% t1n= (px - t2n*a2)/(a1+t2n*a3);
% 
% pr = v1 +t1n*(v3-v1) + t2n*(v2-v1)+t2n*t1n*(v4-v2+v1-v3);
% 
% v1w=1-t1n-t2n+t1n*t2n;
% v2w=t2n-t1n*t2n;
% v3w=t1n-t1n*t2n;
% v4w=t1n*t2n;
% T = [ v1w 0 v2w 0 v3w 0 v4w 0; ...
%          0 v1w 0 v2w 0 v3w 0 v4w ];

[ T ] = getLinTrans( Vv,pv );
prr= T*Vv;
Vv2=Vv;
Vv2(1:2) = Vv2(1:2)+10;
prr2= T*Vv2;
% syms t1 t2
% syms V2x V3x V4x
% syms V2y V3y V4y
% syms Px Py
% eqnx = Px == t1*(V3x) + t2*(V2x)+t2*t1*(V4x-V2x -V3x);
% eqny = Py == t1*(V3y) + t2*(V2y)+t2*t1*(V4y-V2y -V3y);
% 
% [t1sols, t2sols] = solve([eqnx, eqny], [t1, t2]);

scatter(V(:,2),V(:,1) );
set(gca,'YDir','reverse');
hold on;
scatter(p(2), p(1), '*');
scatter(prr2(2), prr2(1), '^');
