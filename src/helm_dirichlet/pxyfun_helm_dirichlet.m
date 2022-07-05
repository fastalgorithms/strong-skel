
% proxy function
function [Kpxy,nbr] = pxyfun_helm_dirichlet(x,slf,nbr,proxy,l,ctr,zpars,nu,area)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling helm_dirichlet_kernel
zstmp = complex([zpars(1);1.0;0.0]);
pxy = bsxfun(@plus,proxy*l,ctr');
Kpxy1 = helm_dirichlet_kernel(pxy,x(:,slf),zstmp,nu(:,slf));
Kpxy1 = bsxfun(@times,Kpxy1,sqrt(area(slf)));
Kpxy3 = bsxfun(@times,helm_dirichlet_kernel(pxy,x(:,slf),zpars,nu(:,slf)),sqrt(area(slf)));
Kpxy = [Kpxy1;Kpxy3];
dx = x(1,nbr) - ctr(1);
dy = x(2,nbr) - ctr(2);
dz = x(3,nbr) - ctr(3);
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);
end
