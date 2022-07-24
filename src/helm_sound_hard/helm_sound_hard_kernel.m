function K_prime = helm_sound_hard_kernel(x, y, z_k, nx)
% KFUN(X,Y,zpars,NU) computes the derivative of the single layer Helmholtz
% potential evaluated pairwise between points in X and points in Y
% (does not handle the singularity). Z_K is the wavenumber and NX is is the
% normal vector at the targets X.

% r = |x-y|
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);

nx = nx';

% e^{ik|x-y|}
zexp = exp(1j*z_k*dr);

% Calculate normal derivative of single layer potential at the targets
rdotnx = bsxfun(@times, dx, nx(:, 1)) + bsxfun(@times, dy, nx(:, 2)) + ...
           bsxfun(@times, dz, nx(:,3));

K_prime = 1/(4*pi).*rdotnx./dr.^3.*zexp.*(1j*z_k*dr-1.0);
K_prime(dr == 0) = 0;
end
