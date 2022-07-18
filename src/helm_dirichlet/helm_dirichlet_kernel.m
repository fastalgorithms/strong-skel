
function K = helm_dirichlet_kernel(x,y,zpuse,nuuse)
% KFUN(X,Y,zpars,NU) computes the Helmholtz potential evaluated
% pairwise between points in X and points in Y (does not handle the
% singularity).  zpars says which layer potential to use
% using the surface normal vectors in NU.
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);

% foo = size(dx)`
% bar = size(nuuse)
zexp = exp(1j*zpuse(1)*dr);
slp = 1/(4*pi)./dr.*zexp;
rdotn = bsxfun(@times,dx,nuuse(1,:)) + bsxfun(@times,dy,nuuse(2,:)) + ...
          bsxfun(@times,dz,nuuse(3,:));
% baz = size(rdotn)
dlp = 1/(4*pi).*rdotn./dr.^3.*zexp.*(1.0-1j*zpuse(1)*dr);
K = zpuse(2)*slp + zpuse(3)*dlp;
K(dr == 0) = 0;

end
