function F = srskelf_asym_new_pts(A,x,occ,rank_or_tol,pxyfun,opts)
    % SRSKELF_ASYM   Asymmetric strong recursive skeletonization factorization.
    %
    %    F = SRSKELF_ASYM(A,X,OCC,RANK_OR_TOL,PXYFUN) produces a factorization
    %    F of the interaction matrix A on the points X using tree occupancy
    %    parameter OCC, local precision parameter RANK_OR_TOL, and proxy
    %    function PXYFUN to capture the far field. This is a function of the
    %    form
    %
    %      [KPXY,NBR] = PXYFUN(X,SLF,NBR,proxy,L,CTR)

    %    that is called for every block, where
    %
    %      - KPXY: interaction matrix against artificial proxy points
    %      - NBR:  block neighbor indices (can be modified)
    %      - X:    input points
    %      - SLF:  block indices
    %      - proxy: proxy points on the unit sphere
    %      - L:    block size
    %      - CTR:  block center
    %
    %    See the examples for further details.
    %
    %    F = SRSKELF_ASYM(A,X,OCC,RANK_OR_TOL,PXYFUN,OPTS) also passes various
    %    options to the algorithm. Valid options include:
    %
    %      - EXT: set the root node extent to [EXT(I,1) EXT(I,2)] along
    %             dimension I.  If EXT is empty (default), then the root extent
    %             is calculated from the data.
    %
    %      - LVLMAX: maximum tree depth (default: LVLMAX = Inf).
    %
    %      - SYMM: assume that the matrix is asymmetric if SYMM = 'N' and
    %              Hermitian positive definite if SYMM = 'P' (default: SYMM =
    %              'N'). If SYMM = 'N', then local factors are computed using
    %              the LU decomposition; if SYMM = 'P', the Cholesky
    %              decomposition.
    %
    %      - VERB: display status of the code if VERB = 1 (default: VERB = 0).
    %

      start = tic;
      % Set sane default parameters
      if nargin < 5
        pxyfun = [];
      end % if
      if nargin < 6
        opts = [];
      end % if
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

      if opts.verb
        disp('This is standard asymmetric srskelf (RS-S).');
      end

      % Check inputs are sensible
        assert(strcmpi(opts.symm,'p') || strcmpi(opts.symm,'n'), ...
             'RSS:srskelf_asym:invalidSymm', ...
             'Symmetry parameter must be ''p'' or ''n''.');

      % Build tree to hold the discretization points
      N = size(x,2);
      tic
      t = shypoct(x,occ,opts.lvlmax,opts.ext);

      if opts.verb
        fprintf(['-'*ones(1,80) '\n'])
        fprintf('%3s | %6s | %8s | %8s | %8s | %8s | %10s (s)\n', ...
                'lvl','nblk','nRemIn','nRemOut','inRatio','outRatio','time')
        % Print summary information about tree construction
        fprintf(['-'*ones(1,80) '\n'])
        fprintf('%3s | %63.2e (s)\n','-',toc)

        % Count the nonempty boxes at each level
        pblk = zeros(t.nlvl+1,1);
        for lvl = 1:t.nlvl
          pblk(lvl+1) = pblk(lvl);
          for i = t.lvp(lvl)+1:t.lvp(lvl+1)
            if ~isempty(t.nodes(i).xi)
              pblk(lvl+1) = pblk(lvl+1) + 1;
            end % if
          end % for
        end % for
      end % if

      % Initialize the data structure holding the factorization
      nbox = t.lvp(end);

      e = cell(nbox,1);
      % Each element of F.factors will contain the following data for one box:
      %   - sk: the skeleton DOF indices
      %   - rd: the redundant DOF indices
      %   - nbr: the neighbor (near-field) DOF indices
      %   - T: the interpolation matrix mapping redundant to skeleton
      %   - E: the left factor of the Schur complement update to sk
      %   - F: the right factor of the Schur complement update to sk
      %   - L: the left factor of the diagonal block
      %   - U: the right factor of the diagonal block
      %   - C: the left factor of the Schur complement update to nbr
      %   - D: the right factor of the Schur complement update to nbr
      F = struct('sk',e,'rd',e,'nbr',e,'T',e,'E',e,'F',e,'L',e,'U',e,'C',e,...
                 'D',e);
      F = struct('N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F,...
                 'symm',opts.symm);
      nlvl = 0;
      n = 0;
      % Mark every DOF as "remaining", i.e., not yet eliminated
      rem = true(N,1);
      lookup_list = zeros(nbox,1);

      % Loop over the levels of the tree from bottom to top
      for lvl = t.nlvl:-1:1
        time = tic;
        nlvl = nlvl + 1;
        nrem1 = sum(rem);

        % For each box, pull up information about skeletons from child boxes
        for i = t.lvp(lvl)+1:t.lvp(lvl+1)
          t.nodes(i).xi = [t.nodes(i).xi [t.nodes(t.nodes(i).chld).xi]];
        end % for

        boxsize = t.lrt/2^(lvl - 1);
        tol = rank_or_tol;

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

        % Loop over each box in this level
        for i = t.lvp(lvl)+1:t.lvp(lvl+1)
          slf = t.nodes(i).xi;
          nbr = [t.nodes(t.nodes(i).nbor).xi];

          nslf = length(slf);
          % Sorting not necessary, but makes debugging easier
          slf = sort(slf);

          nnbr = length(nbr);
          % Sorting not necessary, but makes debugging easier
          nbr = sort(nbr);

          % If we are at the second level (i.e., the first level we reach in a
          % bottom-to-top loop in which there do not exist pairs of
          % non-adjacent boxes) then we can do weak skeletonization, so instead
          % of the interaction list we skeletonize against the neighbor set.
          % Currently turned off skeletonization at level 1 in current version
          % also removed selecting subset of indices from interaction list,
          % needs to be fixed..
          if lvl == 2
            lst = [];
            lst = nbr;
            nbr = [];
            nnbr = 0;
            l = t.lrt/2^(lvl - 1);
          else
            lst = [t.nodes(t.nodes(i).ilist).xi];
            l = t.lrt/2^(lvl - 1);
          end % if

          % Compute proxy interactions and subselect neighbors
          Kpxy = zeros(0,nslf);
          if lvl > 2
             if(~isempty(pxyfun))
                [Kpxy,~] = pxyfun(slf,lst,proxy_dict,l,t.nodes(i).ctr);
             end
          end % if


          nlst = length(lst);
          % Sorting not necessary, but makes debugging easier
          lst = sort(lst);

          % Compute interaction matrix between box and far-field (except level
          % 2, where near-field is included).
          K1 = full(A(lst,slf)); % TODO: Will return block matrix
          if strcmpi(opts.symm,'n')
            K1 = [K1; conj(full(A(slf,lst)))'];
          end % if


          K2 = spget('lst','slf'); % TODO will return in terms of pts
          % EXPAND out indices of K2
          K2_expanded = repelem(K2, 2, 2);
          if strcmpi(opts.symm,'n')
              K2 = [K2; conj(spget('slf','lst'))'];
              tmp = repelem(conj(spget('slf','lst')), 2, 2);
              K2_expanded = [K2_expanded; tmp'];
          end % if

          % RESHAPE: to perform ID over point indices
          K1_reshaped = reshape_pts_idxs(K1, 2, 2);
          % Sometimes Kproxy picks out 0 proxy points? hence if statement
          [m_proxy, ~] = size(Kpxy);
          if m_proxy > 0
            Kpxy_reshaped = reshape_pts_idxs(Kpxy, 4, 2);
            % K = [K1+K2_expanded; Kpxy];
            K = [reshape_pts_idxs(K1+K2_expanded, 2, 2); Kpxy_reshaped];
          else
            % K = [K1+K2_expanded];
            K = [reshape_pts_idxs(K1+K2_expanded, 2, 2)];
          end

        %   size(K2_expanded)
        %   size(Kpxy)
          % TODO: in terms of pt indices
        %   K = [K1+K2; Kpxy];
        %   K = [K1+K2_expanded; Kpxy];

          % Compute the skeleton/redundant points and interpolation matrix
          % TODO reshape in terms of pts
          [sk,rd,T] = id(K,rank_or_tol);

          % Move on to next box if no compression for this box
          if isempty(rd)
            continue
          end % if

          % Otherwise, compute the diagonal and off-diagonal blocks for this
          % box
          % TODO: spget will return in terms of points, will need to reshape
          K  = full(A(slf,slf)) + repelem(spget('slf','slf'), 2, 2);

          K2 = full(A(nbr,slf)) + repelem(spget('nbr','slf'), 2, 2);

          if strcmpi(opts.symm,'n')
            K3 = full(A(slf,nbr)) + repelem(spget('slf','nbr'), 2, 2);
          end % if


          % Skeletonize
          % TODO: rd and sk will be interms of points, will need to expand out.
          sk_expanded = pts_to_mat_idxs(sk, 2);
          rd_expanded = pts_to_mat_idxs(rd, 2);
          T_expanded = repelem(T, 2, 2);

        %   K(rd,:) =  K(rd,:) - conj(T)'*K(sk,:);
        %   K(:,rd) = K(:,rd) - K(:,sk)*T;
        %   K2(:,rd) = K2(:,rd) - K2(:,sk)*T;
        %   if strcmpi(opts.symm,'n')
        %     K3(rd,:) = K3(rd,:) - conj(T)'*K3(sk,:);
        %   end % if

          K(rd_expanded,:) =  K(rd_expanded,:) - conj(T_expanded)'*K(sk_expanded,:);
          K(:,rd_expanded) = K(:,rd_expanded) - K(:,sk_expanded)*T_expanded;
          K2(:,rd_expanded) = K2(:,rd_expanded) - K2(:,sk_expanded)*T_expanded;
          if strcmpi(opts.symm,'n')
            K3(rd_expanded,:) = K3(rd_expanded,:) - conj(T_expanded)'*K3(sk_expanded,:);
          end % if


        %   if strcmpi(opts.symm,'p')
        %     % Cholesky for positive definite input
        %     L = chol(K(rd,rd),'lower');
        %     U = [];
        %     E = K(sk,rd)/conj(L)';
        %     G = [];
        %     C = K2(:,rd)/conj(L)';
        %     D = [];
        %   elseif strcmpi(opts.symm,'n')
        %     % Otherwise, LU
        %     [L,U] = lu(K(rd,rd));
        %     E = K(sk,rd)/U;
        %     G = L\K(rd,sk);
        %     C = K2(:,rd)/U;
        %     D = L\K3(rd,:);
        %   end % if

          if strcmpi(opts.symm,'p')
            % Cholesky for positive definite input
            L = chol(K(rd_expanded,rd_expanded),'lower');
            U = [];
            E = K(sk_expanded,rd_expanded)/conj(L)';
            G = [];
            C = K2(:,rd_expanded)/conj(L)';
            D = [];
          elseif strcmpi(opts.symm,'n')
            % Otherwise, LU
            [L,U] = lu(K(rd_expanded,rd_expanded));
            E = K(sk_expanded,rd_expanded)/U;
            G = L\K(rd_expanded,sk_expanded);
            C = K2(:,rd_expanded)/U;
            D = L\K3(rd_expanded,:);
          end % if

          % Store matrix factors for this box
          n = n + 1;
        %   F.factors(n).sk  = slf(sk);
        %   F.factors(n).rd  = slf(sk);
        % sz_sk = size(sk)
        % sz_rd = size(rd)

          F.factors(n).sk  = slf(sk);
          F.factors(n).rd  = slf(sk);
          F.factors(n).nbr = nbr;
          F.factors(n).T = T;
          F.factors(n).E = E;
          F.factors(n).F = G;
          F.factors(n).L = L;
          F.factors(n).U = U;
          F.factors(n).C = C;
          F.factors(n).D = D;
          % Box number i is at index n (more sensible for non-uniform case)
          lookup_list(i) = n;

          t.nodes(i).xi = slf(sk);
          rem(slf(rd)) = 0;


        end % for
        F.lvp(nlvl+1) = n;

        % Print summary for the latest level
        if opts.verb
          nrem2 = sum(rem);
          nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);
          fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
                  lvl,nblk,nrem1,nrem2,nrem1/nblk,nrem2/nblk,toc(time))
        end % if
      end % for

       % Truncate extra storage, and we are done
      F.factors = F.factors(1:n);
      if opts.verb
        fprintf(['-'*ones(1,80) '\n'])
        toc(start)
      end % if

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

      function A_reshaped = reshape_pts_idxs(A, m, n)
        [M, N] = size(A);
        M_reshaped = M/m;
        N_reshaped = N/n;
        A = permute(reshape(A, m, M_reshaped, n, N_reshaped), [1 3 2 4]);
        A_reshaped = reshape(A, m*n*M_reshaped, N_reshaped);
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

    % function mat_idxs = pts_to_mat_idxs(pts_idxs, n)
    %     offsets = ((0:n-1).') + 1;
    %     mat_idxs = bsxfun(@plus, (pts_idxs - 1)*n, offsets);
    %     mat_idxs = mat_idxs(:);
    % end

    function mat_idxs = pts_to_mat_idxs(pts_idxs, n)

        mat_idxs = zeros(n*length(pts_idxs), 1);
        for i = 1:numel(pts_idxs)
            pt_idx = pts_idxs(i);

            mat_idxs((i-1)*n+1: i*n) = linspace((pt_idx-1)*n+1, pt_idx*n, n);
        end
    end
    end % srskelf_asym