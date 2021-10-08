clear
radii = [1.0;2.0;0.25];
scales = [1.2;1.0;1.7];
nu = 20;
nv = 20;
nosc = 5;
norder = 5;
S = wtorus(radii,scales,nosc,nu,nv,norder);
zpars = complex([1.1; 1j*1.0; 1.0]);
eps = 1e-2;
tic, spmat = helm_near_corr(S,zpars,eps); toc;