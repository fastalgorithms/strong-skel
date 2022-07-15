function S_x = helm_sound_hard_proxy_kernel(x, y, z_k)
%HELM_SOUND_HARD_PROXY_KERNEL Summary of this function goes here
%   Detailed explanation goes here

% r = (x-y), dr = |x-y|
dx = bsxfun(@minus, x(1,:)', y(1,:));
dy = bsxfun(@minus,x(2,:)', y(2,:));
dz = bsxfun(@minus,x(3,:)', y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);
r = x - y;

% e^{ik|x-y|}
zexp = exp(1j*z_k*dr);

% Derivative of single layer potential with respect to targets
S_x = 1/(4*pi)./dr.^3.*zexp.*(1j*z_k*dr-1.0)*r;
S_x(dr == 0) = 0;

end

