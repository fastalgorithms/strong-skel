function [Kpxy,nbr] = pxyfun_helm_transmission_wrapper(x, slf, nbr, proxy, l, ctr, x_or, z_k, nu, area)
% PXYFUN_HELM_SOUND_HARD_WRAPPER Wraps the PXYFUN_HELM_SOUND_HARD function
% which only acts on quadrature points so that it works when applied to the
% discretized system that arises out of the sound hard problem. X_OR are
% the physical points, SLF are the indices of the physical points in a box,
% NBR are the indices of the points in its interaction list, L is its level
% CTR is its centre, Z_K is the complex wavenumber, NU the normal derivatives
% and AREA the weights at the quadrature points.

% Convert system source indices back to original indices
jpts = idivide(int64(slf(:)-1), int64(2))+1;
[juni,~,ijuni] = unique(jpts);

% Convert system neighbour indices back to original indices
nbr_pts = idivide(int64(nbr(:)-1), int64(2))+1;
[nbruni, ~, ~] = unique(nbr_pts);

[Kpxy, ~] = pxyfun_helm_transmission(x_or, juni, nbruni, proxy, l, ctr, z_k, nu, area);


% select relevant columns 
ijuni2 = (ijuni-1)*2 + mod(slf(:)-1, 2)+1;

Kpxy = Kpxy(:, ijuni2);


dx = x(1, nbr) - ctr(1);
dy = x(2, nbr) - ctr(2);
dz = x(3, nbr) - ctr(3);
dist = max(abs([dx;dy;dz]));
nbr = nbr(dist/l > 1.5);


end