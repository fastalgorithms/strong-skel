

      subroutine getnearquad_helm_comb_dir_spmat(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear,irowind,icolind)
c
c       this subroutine generates the near field quadrature
c       for the representation u = (\alpha S_{k}  + \beta D_{k}) ---(1)
c       where the near field is specified by the user 
c       in row sparse compressed format.
c
c       The quadrature is computed by the following strategy
c        targets within a sphere of radius rfac0*rs
c        of a chunk centroid is handled using adaptive integration
c        where rs is the radius of the bounding sphere
c        for the patch
c  
c       All other targets in the near field are handled via
c        oversampled quadrature
c
c       The recommended parameter for rfac0 is 1.25d0
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders - integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            starting location of data on patch i
c  
c         iptype - integer(npatches)
c           type of patch
c           iptype = 1 -> triangular patch discretized with RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c          srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c             srcvals(10:12,i) - normals info
c 
c         ndtarg - integer
c            leading dimension of target array
c        
c         ntarg - integer
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c         ipatch_id - integer(ntarg)
c            id of patch of target i, id = -1, if target is off-surface
c
c         uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if on surface, otherwise
c            set to 0 by default
c            
c          eps - real *8
c             precision requested
c
c          zpars - complex *16 (3)
c              kernel parameters (Referring to formula (1))
c              zpars(1) = k 
c              zpars(2) = alpha
c              zpars(3) = beta
c
c           iquadtype - integer
c              quadrature type
c              iquadtype = 1, use ggq for self + adaptive integration
c                 for rest
c 
c
c           nnz - integer
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer
c               radius parameter for near field
c
c           nquad - integer
c               number of entries in wnear
c
c        output
c            wnear - complex *16(nquad)
c               the desired near field quadrature
c            irowind - integer(nquad)
c              row indices corresponding to near quadrature
c            icolind - integer(nquad)
c              column indices corresponding to near quadrature
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      integer, intent(out) :: irowind(nquad), icolind(nquad)
      complex *16, intent(out) :: wnear(nquad)


      integer ipars
      integer ndd,ndz,ndi
      real *8 dpars

      complex *16 alpha,beta
      integer i,j,jpatch,jstart,jquadstart,l,npols
      integer ipv
      real *8 t1,t2
      
      call cpu_time(t1)
      call getnearquad_helm_comb_dir(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
      call cpu_time(t2)

ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,npols,jquadstart)
ccC$OMP$PRIVATE(jstart,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            irowind(jquadstart+l-1) = i
            icolind(jquadstart+l-1) = jstart+l-1
          enddo
        enddo
      enddo
ccC$OMP END PARALLEL DO      

      return
      end

