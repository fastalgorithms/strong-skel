function T = helm_hypersingular_kernel(x, y, z_k, nx, ny)
% HYPERSINGULAR_KERNEL(X,Y,zpars,NU) computes the hypersingular operator
% for the Helmholtz potential evaluated pairwise between points in X and
% points in Y (does not handle the singularity). Z_K is the wavenumber and
% NX is is the normal vector at the targets X.

% r = |x-y|
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);

% e^{ik|x-y|}
zexp = exp(1j*z_k*dr);

rdotny = bsxfun(@times, dx, ny(1, :)) + bsxfun(@times, dy, ny(2, :)) + bsxfun(@times, dz, ny(3, :));
rdotnx = bsxfun(@times, dx, nx(1,:)') + bsxfun(@times, dy, nx(2,:)') + bsxfun(@times, dz, nx(3,:)');
nxdotny = bsxfun(@times, nx(1, :)', ny(1,:)) + bsxfun(@times, nx(2, :)', ny(2,:)) + bsxfun(@times, nx(3, :)', ny(3,:)) ;

T1 = 1./dr.^3.*nxdotny.*zexp.*(1-1j*z_k*dr);

T2 = 1./dr.^5.*rdotny.*zexp.*((dr.^2)*z_k^2+z_k*3*1j*dr -3).*rdotnx;

%T3 = z_k^2./dr.^5.*rdotny.*zexp.*r3dotnx;

T = 1/(4*pi)*(T1+T2);

T(dr == 0) = 0;

end