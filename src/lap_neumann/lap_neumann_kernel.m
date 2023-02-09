
function K = lap_neumann_kernel(x,y,nuuse)
% KFUN(X,Y,zpars,NU) computes the Helmholtz potential evaluated
% pairwise between points in X and points in Y (does not handle the
% singularity).  zpars says which layer potential to use
% using the surface normal vectors in NU.
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);
rdotn = bsxfun(@times,nuuse(1,:).',dx) + bsxfun(@times,nuuse(2,:).',dy) + ...
          bsxfun(@times,nuuse(3,:).',dz);
K = -1/(4*pi).*rdotn./dr.^3;
K(dr == 0) = 0;
end
