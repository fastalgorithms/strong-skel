% Plot geometry
radii = [1.0;2.0;0.25];
scales = [1.0;1.0;1.0];

m = 1000; n = 500;
uu = linspace(0, 2*pi, m);
vv = linspace(0, 2*pi, n);
rr = 1;
[Uu, Vv] = meshgrid(uu, vv);

x = (rr.*cos(Uu) + 2 + 0.25*cos(5*Vv)).*cos(Vv)*1.0;
y = (rr.*cos(Uu) + 2 + 0.25*cos(5*Vv)).*sin(Vv)*1.0;
z = rr.*sin(Uu)*1.0;

figure;
h = surf(x,y,z);
axis off;
h.LineStyle='none';

colormap('copper');
light('Position',[-1 0 0], 'Style','infinite');
lighting gouraud;
material shiny;

% Solve a sample problem
[X, x_or] = ie_fmm3dbie_sound_hard(15, 4, 1000, 1, 5e-7, 50);

% Plot solution
figure;
xflat = reshape(x, 1, []);
yflat = reshape(y, 1, []);
zflat = reshape(z, 1, []);
xyz = [xflat; yflat; zflat];
xyz = xyz';
x_or_t = x_or';
[k, dist]  = dsearchn(x_or_t, xyz);
X_int = X(k);
X_int = reshape(X_int, [n, m]);
h = surf(x, y, z);
axis off;
colorbar;
h.CData = real(X_int);
h.FaceColor = 'flat';
h.LineStyle = 'none';
