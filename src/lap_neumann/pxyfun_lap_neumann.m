
% proxy function
function [Kpxy,nbr] = pxyfun_lap_neumann(x,slf,nbr,proxy,l,ctr,area)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then

pxy = bsxfun(@plus,proxy*l,ctr');

dx = bsxfun(@minus,pxy(1,:)',x(1,slf));
dy = bsxfun(@minus,pxy(2,:)',x(2,slf));
dz = bsxfun(@minus,pxy(3,:)',x(3,slf));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);
Kpxy1 = 1./(4*pi*dr);
Kpxy2 = -Kpxy1./dr.^2.*dx;
Kpxy3 = -Kpxy1./dr.^2.*dy;
Kpxy4 = -Kpxy1./dr.^2.*dz;
Kpxy1 = bsxfun(@times,Kpxy1,sqrt(area(slf)));
Kpxy2 = bsxfun(@times,Kpxy2,sqrt(area(slf)));
Kpxy3 = bsxfun(@times,Kpxy3,sqrt(area(slf)));
Kpxy4 = bsxfun(@times,Kpxy4,sqrt(area(slf)));

Kpxy = [Kpxy1;Kpxy2;Kpxy3;Kpxy4];
dx = x(1,nbr) - ctr(1);
dy = x(2,nbr) - ctr(2);
dz = x(3,nbr) - ctr(3);
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);
end
