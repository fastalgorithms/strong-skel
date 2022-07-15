function A_wrapped = Afun_helm_sound_hard_wrapper(i, j, x, zpars, nu, area, P, S)
% input indices are non-physical, have to work out physical indices.

ipts = idivide(int64(i(:)-1),int64(2))+1;
jpts = idivide(int64(j(:)-1),int64(2))+1;

[iuni,~,iiuni] = unique(ipts);
[juni,~,ijuni] = unique(jpts);

% Get matrix valued entries of system matrix for unique points
A_uni = Afun_helm_sound_hard(iuni, juni, x, zpars, nu, area, P, S);

% relevant rows and columns in A_uni
iiuni2 = (iiuni-1)*2 + mod(i(:)-1, 2)+1;
ijuni2 = (ijuni-1)*2 + mod(j(:)-1, 2)+1;

A_wrapped = A_uni(iiuni2, ijuni2);

end

