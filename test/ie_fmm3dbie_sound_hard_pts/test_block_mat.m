% Create an example block matrix and ordinary matrix
n = 2; m = 2; p = 3; q = 3; % n, m: block size; p, q: number of blocks

% Creating block matrix
blockMatrix = rand(n*p, m*q)

% Creating ordinary matrix
ordMatrix = rand(p,q)

% Expand each element of ordMatrix to form a block
ordMatrix_exp = repelem(ordMatrix, n, m)

% Add blockMatrix and ordMatrix_exp
resultMatrix = blockMatrix + ordMatrix_exp;