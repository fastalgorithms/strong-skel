function K_prime = helm_sound_hard_proxy_kernel(x, y, z_k)
%HELM_SOUND_HARD_PROXY_KERNEL Computes the normal derivative of the single
%layer potential between the targets X and sources Y, where the normal
% derivative is taken with respect to the targets, where Z_K is the complex
% wavenumber.

% r = (x-y), dr = |x-y|
dx = bsxfun(@minus, x(1,:)', y(1,:));
dy = bsxfun(@minus,x(2,:)', y(2,:));
dz = bsxfun(@minus,x(3,:)', y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);

% e^{ik|x-y|}
zexp = exp(1j*z_k*dr);

% Derivative of single layer potential with respect to targets
K_prime = 1/(4*pi)./dr.^3.*zexp.*(1j*z_k*dr-1.0);

% Exclude singularity
K_prime(dr == 0) = 0;

K = 1/(4*pi)*zexp./dr;
K(dr == 0) = 0;

% Return three components wrt to r = x-y
K_prime = {K_prime.*dx, K_prime.*dy, K_prime.*dy, K};
end