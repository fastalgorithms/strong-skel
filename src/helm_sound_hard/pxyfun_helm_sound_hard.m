% proxy function
function [Kpxy,nbr] = pxyfun_helm_sound_hard(x, slf, nbr, proxy, l, ctr, z_k, nu, area)
% PXYFUN_HELM_SOUND_HARD(X,SLF,NBR,PROXY,L,CTR,ZPARS,NU,AREA) computes 
% interactions between the points X(:,SLF) and the set of proxy points by 
% scaling the proxy sphere to  appropriately contain a box at level L 
% centered at CTR and then calling HELM_SOUND_HARD_PROXY kernel. Z_K is the
% complex wavenumber, NU are the normal derivatives and AREA are the
% weights at each quadrature point.

% Shift and scale precomputed proxy surface
proxy = bsxfun(@plus, proxy*l, ctr');

% Capture interaction over each block, compress the following
targets = proxy;
sources = x(:, slf);

% Allocate a stacked proxy matrix of all kernels to compress
i = length(targets);
j = length(slf);
Kpxy = zeros(4*2*i, 2*j);

% Exterior kernels to compress
K_prime = helm_sound_hard_proxy_kernel(targets, sources, z_k);
K_prime_1 = K_prime{1};
K_prime_1 = bsxfun(@times, K_prime_1, area(j));
K_prime_2 = K_prime{2};
K_prime_2 = bsxfun(@times, K_prime_2, area(j));
K_prime_3 = K_prime{3};
K_prime_3 = bsxfun(@times, K_prime_3, area(j));

% Interior kernel to compress
zdtmp = complex([z_k 0.0 1.0]);
D = helm_dirichlet_kernel(targets, sources, zdtmp, nu(:, slf));
D = bsxfun(@times, D, sqrt(area(j)));

Kpxy(1:8:end, 1:2:end) = K_prime_1;
Kpxy(2:8:end, 1:2:end) = K_prime_1;
Kpxy(1:8:end, 2:2:end) = K_prime_1;

Kpxy(3:8:end, 1:2:end) = K_prime_2;
Kpxy(4:8:end, 1:2:end) = K_prime_2;
Kpxy(3:8:end, 2:2:end) = K_prime_2;

Kpxy(5:8:end, 1:2:end) = K_prime_3;
Kpxy(6:8:end, 1:2:end) = K_prime_3;
Kpxy(5:8:end, 2:2:end) = K_prime_3;

Kpxy(7:8:end, 1:2:end) = D;
Kpxy(8:8:end, 1:2:end) = D;
Kpxy(7:8:end, 2:2:end) = D;

% Compute the neighbours
dx = x(1, nbr) - ctr(1);
dy = x(2, nbr) - ctr(2);
dz = x(3, nbr) - ctr(3);
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);

end
