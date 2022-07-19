function A = Afun_helm_sound_hard(i, j, x, z_k, nu, area, P, S)
% AFUN_HELM_SOUND_HARD(I,J,X,Z_K,NU,AREA,P,S) computes entries of the 
% matrix A to be factorized at the index sets I and J.  This handles the 
% near-field correction. X are the quadrature points, with normal
% derivatives NU, weights AREA. P is the permutation matrix required by the
% skeletonization algorithm and S is a matrix of quadrature corrections for
% the near field.

if isempty(i) || isempty(j)
  A = zeros(2*length(i), 2*length(j));
  return
end

[I, J] = ndgrid(i, j);

% Assemble (block) system matrix for sound hard problem
A = zeros(2*length(i), 2*length(j));

% Compute normal derivative of single layer potential wrt to the targets
targets = x(:,i);
sources = x(:,j);
nx = nu(:,i);
K_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k, nx), area(j));

% A11 
A(1:2:end, 1:2:end) = -1j*z_k*K_prime;

% A12
A(1:2:end, 2:2:end) = K_prime;

% A21
A(2:2:end, 1:2:end) = K_prime;

% A22
A(2:2:end, 2:2:end) = 0;

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
tmp(I == J) = tmp(I == J) + (1j*z_k*0.5 - 0.25);
A(1:2:end, 1:2:end) = tmp;

tmp = A(2:2:end, 2:2:end);
tmp(I == J) = complex(-1);
A(2:2:end, 2:2:end) = tmp;

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
