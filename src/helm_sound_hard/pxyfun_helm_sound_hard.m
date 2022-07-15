% proxy function
function [Kpxy,nbr] = pxyfun_helm_sound_hard(x, slf, nbr, proxy, l, ctr, zpars, nu, area)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling helm_dirichlet_kernel

% Shift and scale precomputed proxy surface
pxy = bsxfun(@plus, proxy*l, ctr');

% Complex wavenumber
z_k = zpars(1);

% zstmp = complex([zpars(1);1.0;0.0]);
% Kpxy1 = helm_dirichlet_kernel(pxy, x(:,slf), zstmp, nu(:,slf));
% Kpxy1 = bsxfun(@times,Kpxy1,sqrt(area(slf)));
% Kpxy3 = bsxfun(@times,helm_dirichlet_kernel(pxy,x(:,slf),zpars,nu(:,slf)),sqrt(area(slf)));
% Kpxy = [Kpxy1;Kpxy3];

zdtmp = complex([zpars(1) 0.0 1.0]);

% Capture interaction over each block, compress the following
targets = pxy;
sources = x(:, slf);

i = length(targets);
j = length(slf);

Kpxy = zeros(4*2*i, 2*j);

for ii = 1:i
    for jj = 1:j

        % Exterior kernels to compress
%         size(sources)
%         jj
%         size(pxy)
%         ii
%         disp('\n')
        S_x = helm_sound_hard_proxy_kernel(pxy(:,ii), sources(:,jj), z_k);

        S_x_1 = S_x(1,:);
        S_x_2 = S_x(2,:);
        S_x_3 = S_x(3,:);

        % Interior kernel to compress
        D = helm_dirichlet_kernel(pxy(:,ii), sources(:,jj), zdtmp, nu(:,jj));

        Kpxy_block = [
             [S_x_1  S_x_1; S_x_1 0];
             [S_x_2  S_x_2; S_x_2 0];
             [S_x_3  S_x_3; S_x_3 0];
             [D D; D 0]
         ];
        
        r_idx = ii*2*4-7;
        c_idx = jj*2-1;
       
        Kpxy(r_idx:r_idx+7, c_idx:c_idx+1) = Kpxy_block;
    end
end

dx = x(1, nbr) - ctr(1);
dy = x(2, nbr) - ctr(2);
dz = x(3, nbr) - ctr(3);
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);

end
