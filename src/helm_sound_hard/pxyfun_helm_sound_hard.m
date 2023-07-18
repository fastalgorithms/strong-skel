% proxy function
function [Kpxy,nbr] = pxyfun_helm_sound_hard(x, slf, nbr, proxy_dict, l, ctr, z_k, nu, area)
% PXYFUN_HELM_SOUND_HARD(X,SLF,NBR,PROXY,L,CTR,ZPARS,NU,AREA) computes 
% interactions between the points X(:,SLF) and the set of proxy points by 
% scaling the proxy sphere to  appropriately contain a box at level L 
% centered at CTR and then calling HELM_SOUND_HARD_PROXY kernel. Z_K is the
% complex wavenumber, NU are the normal derivatives and AREA are the
% weights at each quadrature point.

proxy = proxy_dict.proxy;
weigt = l*proxy_dict.weigt;
norms = proxy_dict.norms;
% Shift and scale precomputed proxy surface
proxy = bsxfun(@plus, proxy*l, ctr');

% Capture interaction over each block, compress the following
targets = proxy;
sources = x(:, slf);

% Allocate a stacked proxy matrix of all kernels to compress
i = length(targets);
j = length(slf);


% Exterior kernels to compress
K_prime = helm_sound_hard_proxy_kernel(targets, sources, z_k);
K_prime_1 = K_prime{1};
K_prime_1 = bsxfun(@times, K_prime_1, sqrt(area(j)));

K_prime_1 = bsxfun(@times, weigt.', K_prime_1);
K_prime_2 = K_prime{2};
K_prime_2 = bsxfun(@times, K_prime_2, sqrt(area(j)));
K_prime_2 = bsxfun(@times, weigt.', K_prime_2);
K_prime_3 = K_prime{3};
K_prime_3 = bsxfun(@times, K_prime_3, sqrt(area(j)));
K_prime_3 = bsxfun(@times, weigt.', K_prime_3);

Kmat = bsxfun(@times, norms(1,:).', K_prime_1)+...
    bsxfun(@times, norms(2,:).', K_prime_2)+...
    bsxfun(@times, norms(3,:).', K_prime_3);

% Interior kernel to compress
zdtmp = complex([z_k 0.0 1.0]);

D = helm_dirichlet_kernel(targets, sources, zdtmp, nu(:, slf));
D = bsxfun(@times, D, sqrt(area(j)));
D = bsxfun(@times, weigt.', D);


Kpxy = zeros(4*i, 2*j);
Kpxy(2:4:end,1:2:end) = Kmat;
Kpxy(1:4:end,2:2:end) = Kmat;
Kpxy(4:4:end, 1:2:end) = D;
Kpxy(3:4:end, 2:2:end) = D;
ctruse = ctr(:);
dxyz = abs(x(1:3,nbr)-ctruse(1:3))/l;
nbr = nbr(max(dxyz) < 2.5);

end
