function test_shapes()

    % Setup
    addpath('../../fortran_src/')
    addpath('../../../strong-skel/fortran_src/') %git
    addpath('../../')
    addpath('../../sv')
    addpath('../../mv')
    addpath('../../src')
    addpath('../../src/helm_dirichlet/'); %local
    addpath('../../src/helm_sound_hard/'); %local
    run('../../../FLAM/startup.m'); %git

    npu = 5;
    norder = 5;
    occ = 100;
    z_k = 2.0;
    rank_or_tol = 1E-3;
    m = 5;

    % Initialize a wiggly torus
radii = [1.0;2.0;0.25];
scales = [1.0;1.0;1.0];
nnu = npu;
nnv = npu;
nosc = 5;
sinfo = wtorus(radii, scales, nosc, nnu, nnv, norder);

% Pick 'm' points inside object
xyz_in = zeros(3, m);
rng(42);

uu = rand(m,1)*pi/4 -pi/8;
xyz_in(3,:) = radii(1)*sin(uu)*scales(3);
vv = rand(m,1)*2*pi;
rmin = radii(2) + radii(3)*cos(nosc*vv) - radii(1)*abs(cos(uu));
rmax = radii(2) + radii(3)*cos(nosc*vv) + radii(1)*abs(cos(uu));

rr = rand(m,1)*0.1 + 0.45;
rruse = rmin + rr.*(rmax-rmin);
xyz_in(1,:) = rruse.*cos(vv)*scales(1);
xyz_in(2,:) = rruse.*sin(vv)*scales(2);

% Pick 'm' points on parallel surface to the object
xyz_out = zeros(3, m);
uu = rand(m, 1)*2*pi;

zmin = 2.0*radii(1)*scales(3);
zmax = 4.0*radii(1)*scales(3);
xyz_out(3,:) = zmin + rand(m,1)*(zmax-zmin);
rr2 = rand(m,1)*zmax;
xyz_out(1,:) = rr2.*cos(uu);
xyz_out(2,:) = rr2.*sin(uu);

q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);

% Quadrature points
x = sinfo.srcvals(1:3,:);

nu = sinfo.srcvals(10:12,:);
area = sinfo.wts';
N = sinfo.npts;

% Complex parameterizations
zpars = complex([z_k; -1j*z_k; 1.0]);
zstmp = complex([z_k;1.0;0.0]);
zdtmp = complex([z_k;0.0;1.0]);
eps = rank_or_tol;

% Compute the quadrature corrections
tic, Q = helm_neu_near_corr(sinfo, zpars(1:2), eps); tquad=  toc;
S = Q{1};
P = zeros(N,1);
w = whos('S');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)

% Set system matrix, operating on points indices
Afun = @(i, j) Afun_helm_sound_hard(i, j, x, z_k, nu, area, P, S);

% Set proxy function
pxyfun = @(slf, nbr, proxy_dict, l, ctr) pxyfun_helm_sound_hard(x, slf, nbr, proxy_dict, l, ctr, z_k, nu, area);

% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'z_k', z_k);


% Build tree to hold the discretization points
if ~isfield(opts,'ext')
    opts.ext = [];
end % if
if ~isfield(opts,'lvlmax')
    opts.lvlmax = Inf;
end % if
if ~isfield(opts,'symm')
    opts.symm = 'n';
end % if
if ~isfield(opts,'verb')
    opts.verb = 0;
end % if
if ~isfield(opts,'zk')
    opts.zk = 1.0;
end % if
N = size(x,2);
t = shypoct(x,occ,opts.lvlmax,opts.ext);

lvl = 3;
nlvl = 0;
rem = true(N, 1);
nrem1 = sum(rem);

% For each box, pull up information about skeletons from child boxes
for i = t.lvp(lvl)+1:t.lvp(lvl+1)
    t.nodes(i).xi = [t.nodes(i).xi [t.nodes(t.nodes(i).chld).xi]];
end % for

boxsize = t.lrt/2^(lvl - 1);
tol = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nterms = h3dterms(boxsize,opts.zk,tol);
RR = 5/2;
nleg = 2*ceil((nterms+1));
[xleg,wleg] = lege.exps(nleg);
[XL,YL] = meshgrid(xleg,xleg);
[WX,WY] = meshgrid(wleg,wleg);
XL = RR*XL(:);
YL = RR*YL(:);
ZL = ones(size(XL));
WL = RR*sqrt(WX(:).*WY(:));

DO = ones(size(ZL));
DZ = zeros(size(ZL));

ppts = [XL,YL,RR*ZL];
nrml = [DZ,DZ,DO];
wpts = WL;
%
ppts = [ppts; XL,YL,-RR*ZL];
nrml = [nrml;DZ,DZ,-DO];
wpts = [wpts;WL];
%
ppts = [ppts; XL,-RR*ZL,YL];
nrml = [nrml;DZ,-DO,DZ];
wpts = [wpts;WL];
%
ppts = [ppts; XL,RR*ZL,YL];
nrml = [nrml;DZ,DO,DZ];
wpts = [wpts;WL];
%
ppts = [ppts; -RR*ZL,XL,YL];
nrml = [nrml;-DO,DZ,DZ];
wpts = [wpts;WL];
%
ppts = [ppts; RR*ZL,XL,YL];
nrml = [nrml;DO,DZ,DZ];
wpts = [wpts;WL];

proxy_dict = [];
proxy_dict.proxy = ppts';
proxy_dict.weight = wpts';
proxy_dict.norms = nrml';

%%%%%%%%%%%%%%%%%%

ivec = t.lvp(lvl)+1:t.lvp(lvl+1);
i = ivec(1);
slf = t.nodes(i).xi;
nslf = length(slf);
slf = sort(slf);
nbr = [t.nodes(t.nodes(i).nbor).xi];
nnbr = length(nbr);
nbr = sort(nbr);

if lvl == 2
    lst = nbr;
    nbr = [];
    nnbr = 0;
    l = t.lrt/2^(lvl - 1);
else
    lst = [t.nodes(t.nodes(i).ilist).xi];
    l = t.lrt/2^(lvl - 1);
end % if

nlst = length(lst);
nbox = t.lvp(end);
lookup_list = zeros(nbox, 1);

% Compute proxy interactions and subselect neighbors
Kpxy = zeros(0,nslf);
if lvl > 2
    if(~isempty(pxyfun))
        [Kpxy,~] = pxyfun(slf,lst,proxy_dict,l,t.nodes(i).ctr);
    end
end % if

K1 = full(Afun(lst, slf));
% if strcmpi(opts.symm,'n')
%     K1 = [K1; conj(full(Afun(slf,lst)))'];
% end % if

K2 = full(spget('lst', 'slf'));
K2_expanded = repelem(K2, 2, 2);



% % Both Kpxy and K1 need to be reshaped in terms of points so as to do point ID
tic;
K1_reshaped = reshape_pts_idxs(K1, 2, 2);
Kpxy_reshaped = reshape_pts_idxs(Kpxy, 4, 2);
toc;

% size(K1)
% size(K2_expanded)
size(Kpxy)
size(Kpxy_reshaped)
% size(K1)
% % size(K1_reshaped)
% % size(Kpxy_reshaped)

% size(K1_reshaped)
% tol = 1e-2;
% [sk, rd, T, niter] = id(K1, tol);
% % map to physical indices
% [sk, ~, ~] = unique(idivide(int64(sk(:)-1), int64(2))+1);

% [skr, rdr, Tr, niterr] = id(K1_reshaped, tol);

% sort(sk')
% sort(skr)


    pts_idxs = [1];
    mat_idxs = pts_to_mat_idxs(pts_idxs, 3)


    % function A_reshaped = reshape_pts_idxs(A, m, n)

    %     [M, N] = size(A);

    %     M_reshaped = idivide(int64(M), int64(m));
    %     N_reshaped = idivide(int64(N), int64(n))

    %     A_reshaped = zeros(m*n*M_reshaped, N_reshaped);

    %     for i=1:M_reshaped
    %         for j=1:N_reshaped
    %             block = A(m*(i-1)+1:m*i, n*(j-1)+1:n*j);
    %             A_reshaped(m*n*(i-1)+1:m*n*i, (j-1)+1:j) = block(:);
    %         end
    %     end

    % end

    % function mat_idxs = pts_to_mat_idxs(pts_idxs, n)

    %     mat_idxs = zeros(n*length(pts_idxs), 1);
    %     for i = 1:numel(pts_idxs)
    %         pt_idx = pts_idxs(i);

    %         mat_idxs((i-1)*n+1: i*n) = linspace((pt_idx-1)*n+1, pt_idx*n, n)
    %     end
    % end

    function mat_idxs = pts_to_mat_idxs(pts_idxs, n)
        offsets = ((0:n-1).') + 1;
        mat_idxs = bsxfun(@plus, (pts_idxs - 1)*n, offsets);
        mat_idxs = mat_idxs(:);
    end


    function A_reshaped = reshape_pts_idxs(A, m, n)
        [M, N] = size(A);
        M_reshaped = M/m;
        N_reshaped = N/n;
        mnM_reshaped = m*n*M_reshaped;
        A = permute(reshape(A, m, M_reshaped, n, N_reshaped), [1 3 2 4]);
        A_reshaped = reshape(A, mnM_reshaped, N_reshaped);
    end

    function A = inv_reshape_pts_idxs(A_reshaped, m, n)
        [mnM_reshaped, N_reshaped] = size(A_reshaped);
        M = mnM_reshaped/n;
        N = N_reshaped*n;
        M_block = M/m;
        N_block = N/n;

        % Reshape back to original block structure
        B = reshape(A_reshaped, [m, n, M_block, N_block]);

        % Permute back into original order
        B = permute(B, [1, 3, 2, 4]);

        % Reshape to orignal matrix dims
        A = reshape(B, [M, N]);
    end


    function A = spget(Ityp,Jtyp)
        % A = SPGET(ITYP,JTYP) Sparse matrix access function (native MATLAB is
        % slow for large matrices).  We grab the accumulated Schur complement
        % updates to a block of the matrix from previously-skeletonized
        % levels.  Index sets ITYP and JTYP can be 'slf', 'nbr', or 'lst'.

        % Translate input strings to index sets (and their lengths)
        if strcmpi(Ityp,'slf')
          I_ = slf;
          m_ = nslf;
        elseif strcmpi(Ityp,'nbr')
          I_ = nbr;
          m_ = nnbr;
        elseif strcmpi(Ityp,'lst')
          I_ = lst;
          m_ = nlst;
        end % if

        if strcmpi(Jtyp,'slf')
          J_ = slf;
          n_ = nslf;
        elseif strcmpi(Jtyp,'nbr')
          J_ = nbr;
          n_ = nnbr;
        elseif strcmpi(Jtyp,'lst')
          J_ = lst;
          n_ = nlst;
        end % if

        A = zeros(m_,n_);
        update_list = false(nbox,1);

        % TODO: This function is accepting a variable 'i' that is being mutated globally
        get_update_list(i);


        update_list = lookup_list(flip(find(update_list)'));
        update_list = update_list(update_list~=0)';


        foo = 0;
        for jj = update_list
          g = F.factors(jj);
          xj = [g.sk, g.nbr];
          f = length(g.sk);


          if strcmpi(Ityp,Jtyp)
            % For diagonal block
            idxI = ismembc2(xj,I_); % check which skel+neighbour pts in I_ (contains both sk and rd for current box)
            tmp1 = idxI~=0;
            subI = idxI(tmp1); % Pick out skel+nb in I_
            idxI1 = tmp1(1:f); % Pick out skel points in tmp1
            idxI2 = tmp1(f+1:end); % Pick out nbr pts in tmp1


            tmp1 = [g.E(idxI1,:); g.C(idxI2,:)];

            % Different factorization depending on symmetry
            if strcmpi(opts.symm,'p')
              A(subI, subI) = A(subI,subI) - tmp1*tmp1';
            elseif strcmpi(opts.symm,'n')
              tmp2 = [g.F(:,idxI1), g.D(:,idxI2)];
              A(subI, subI) = A(subI,subI) - tmp1*tmp2;
            end % if
          else
            % For off-diagonal block
            idxI = ismembc2(xj,I_);
            idxJ = ismembc2(xj,J_);

            tmp1 = idxI~=0;
            tmp2 = idxJ~=0;

            subI = idxI(tmp1);
            subJ = idxJ(tmp2);
            idxI1 = tmp1(1:f);
            idxI2 = tmp1(f+1:end);
            idxJ1 = tmp2(1:f);
            idxJ2 = tmp2(f+1:end);

            tmp1 = [g.E(idxI1,:); g.C(idxI2,:)];
            % Different factorization depending on symmetry
            if strcmpi(opts.symm,'p')
              tmp2 = [g.E(idxJ1,:); g.C(idxJ2,:)]';
            elseif strcmpi(opts.symm,'n')
              tmp2 = [g.F(:,idxJ1), g.D(:,idxJ2)];
            end % if
            A(subI, subJ) = A(subI,subJ) - tmp1*tmp2;
          end % if
        end % for

        function get_update_list(node_idx)
          % GET_UPDATE_LIST(NODE_IDX) Recursively get the list of all nodes in
          % the tree that could have generated Schur complement updates to
          % points in node NODE_IDX
          update_list(node_idx) = 1;
          update_list(t.nodes(node_idx).snbor) = 1;
          for k = t.nodes(node_idx).chld
            get_update_list(k);
          end % for
        end % get_update_list
      end % spget

end

