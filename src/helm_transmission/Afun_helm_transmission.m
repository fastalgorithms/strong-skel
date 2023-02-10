function A = Afun_helm_transmission(i, j, x, z_k, z_k0, nu, area, P, S, mu, mu0)
% AFUN_HELM_TRANSMISSION(I,J,X,Z_K,NU,AREA,P,S) computes entries of the
% matrix A to be factorized at the index sets I and J.  This handles the
% near-field correction. X are the quadrature points, with normal
% derivatives NU, weights AREA. P is the permutation matrix required by the
% skeletonization algorithm and S is a matrix of quadrature corrections for
% the near field.

if isempty(i) || isempty(j)
  A = zeros(2*length(i), 2*length(j));
  return
end

% Assemble (block) system matrix for transmission problem
A = zeros(2*length(i), 2*length(j));

% Compute block operators for transmission problem
targets = x(:,i);
sources = x(:,j);
nx = nu(:,i);
ny = nu(:,j);
zstmp = [z_k; 1; 0];
zdtmp = [z_k; 0; 1];
K_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k, nx), area(j));
K = bsxfun(@times, helm_dirichlet_kernel(targets, sources, zdtmp, ny), area(j));
Sl = bsxfun(@times, helm_dirichlet_kernel(targets, sources, zstmp, ny), area(j));
T = bsxfun(@times, helm_hypersingular_kernel(targets, sources, z_k, nx, ny), area(j));

z0stmp = [z_k0; 1; 0];
z0dtmp = [z_k0; 0; 1];
K0_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k0, nx), area(j));
K0 = bsxfun(@times, helm_dirichlet_kernel(targets, sources, z0dtmp, ny), area(j));
Sl0 = bsxfun(@times, helm_dirichlet_kernel(targets, sources, z0stmp, ny), area(j));
T0 = bsxfun(@times, helm_hypersingular_kernel(targets, sources, z_k0, nx, ny), area(j));

% A11
A(1:2:end, 1:2:end) = -(mu*K-mu0*K0);

% A12
A(1:2:end, 2:2:end) = -(mu^2*Sl-mu0^2*Sl0);

% A21
A(2:2:end, 1:2:end) = T-T0;

% A22
A(2:2:end, 2:2:end) = mu*K_prime-mu0*K0_prime;

% Compute quadrature corrections
M = spget_quadcorr(i, j, P, S);
idx = abs(M) ~= 0;

% Add quadrature corrections to each block
tmp = A(1:2:end, 1:2:end);
tmp(idx) = M(idx);
A(1:2:end, 1:2:end) = tmp;

tmp = A(1:2:end, 2:2:end);
tmp(idx) = M(idx);
A(1:2:end, 2:2:end) = tmp;

tmp = A(2:2:end, 1:2:end);
tmp(idx) = M(idx);
A(2:2:end, 1:2:end) = tmp;

tmp = A(2:2:end, 2:2:end);
tmp(idx) = M(idx);
A(2:2:end, 2:2:end) = tmp;

% Add identity parts to A11 and A22 blocks
tmp = A(1:2:end, 1:2:end);
tmp(I == J) = tmp(I == J) + mu0+mu;
A(1:2:end, 1:2:end) = tmp;

tmp = A(2:2:end, 2:2:end);
tmp(I == J) = tmp(I == J) + mu0 + mu;

% Compute weights for each block
tmp = A(1:2:end, 1:2:end);
tmp = bsxfun(@times, sqrt(area(i)).',tmp);
A(1:2:end, 1:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

tmp = A(1:2:end, 2:2:end);
tmp = bsxfun(@times, sqrt(area(i)).',tmp);
A(1:2:end, 2:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

tmp = A(2:2:end, 1:2:end);
tmp = bsxfun(@times, sqrt(area(i)).',tmp);
A(2:2:end, 1:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

tmp = A(2:2:end, 2:2:end);
tmp = bsxfun(@times, sqrt(area(i)).',tmp);
A(2:2:end, 2:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));
end
