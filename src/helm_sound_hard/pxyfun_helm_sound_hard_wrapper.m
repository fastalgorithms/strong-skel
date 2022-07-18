function [Kpxy,nbr] = pxyfun_helm_sound_hard_wrapper(x, slf, nbr, proxy, l, ctr, zpars, nu, area)

% (physical) targets
% pxy = bsxfun(@plus, proxy*l, ctr');

% (non-physical) sources
jpts = idivide(int64(slf(:)-1), int64(2))+1;
[juni,~,ijuni] = unique(jpts);

[Kpxy, nbr] = pxyfun_helm_sound_hard(x, juni, nbr, proxy, l, ctr, zpars, nu, area);

% select relevant columns 
ijuni2 = (ijuni-1)*2 + mod(slf(:)-1, 2)+1;

Kpxy = Kpxy(:, ijuni2);

end