C ##################################################################
C ##     Output H, S, SigmaL, SigmaR in selected (iso) region     ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

#ifndef CDF     
      subroutine OutputRegion(USEBULK,ORTHO,joutfile,Voltage,NEn,
     .     ispin,nspin,isobegin,isoend,nou,
     .     Hi,Si,NGL,HLL,NGR,HRR,ZEnergy)

C INPUT:
      logical USEBULK
      logical ORTHO             !if true do Lowdin orthogonalization
c! if true self-energy only is input else
c                               z*S-H-Sigma for bulk is in sfe

      integer joutfile          !outfile unit 
      integer nou               ! no : #used states (minus buffers)
      complex*16 Hi(0:nou*nou-1), Si(0:nou*nou-1)
      complex*16 ZEnergy 
      integer NGL,NGR
      complex*16 HLL(NGL,NGL)   ! Left electrode GF
      complex*16 HRR(NGR,NGR)   ! Right electrode GF
      integer isobegin,isoend   ! beginning/end of contact region (index starting w. NGL+1)
      integer NEn               ! No. energy points to output
      real*8  Voltage           ! The applied voltage
      integer nspin,ispin

      return 
      end

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#else ! Netcdf defined 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine OutputRegion(USEBULK,ORTHO,joutfile,Voltage,NEn,
     .     ispin,nspin,isobegin,isoend,nou,
     .     Hi,Si,NGL,HLL,NGR,HRR,ZEnergy)
      
      use netcdf

      implicit none

C ======================================================================
c     Note: Most indices are starting at 0 in this routine
C ======================================================================
C INPUT:
      logical USEBULK
      logical ORTHO             !if true do Lowdin orthogonalization
c! if true self-energy only is input else
c                               z*S-H-Sigma for bulk is in sfe

      integer joutfile          !outfile unit 
      integer nou               ! no : #used states (minus buffers)
      complex*16 Hi(0:nou*nou-1), Si(0:nou*nou-1)
      complex*16 ZEnergy 
      integer NGL,NGR
      complex*16 HLL(NGL,NGL)   ! Left electrode GF
      complex*16 HRR(NGR,NGR)   ! Right electrode GF
      integer isobegin,isoend   ! beginning/end of contact region (index starting w. NGL+1)
      integer NEn               ! No. energy points to output
      real*8  Voltage           ! The applied voltage
      integer nspin,ispin

C ======================================================================
C OUTPUT: to netcdf file
      complex*16, dimension(:), allocatable ::
     . Hc, Sc,H, invG, SigmaL, SigmaR
      real*8 Energy(2)
      real*8, dimension(:,:), allocatable :: RMat
c================================================================
c NetCDF: Begin
      integer  iret             !error status return
      integer  ncid             ! NetCDF id
* dimension ids
      integer  nc_dim, NEn_dim, Complex_dim, Real_dim
* variable ids
      integer  H_id, S_id
      integer ReSigmaL_id, ImSigmaL_id, ReSigmaR_id, ImSigmaR_id
      integer  En_id, Volt_id
* rank (number of dimensions) for each variable
c      integer  Sigma_rank
c      integer  En_rank
c      parameter (Sigma_rank = 3)
c      parameter (En_rank = 2)
* variable shapes
      integer  En_dims(2),H_dims(2),Sigma_dims(3), Voltage_dim

      save ncid, ReSigmaL_id, ImSigmaL_id
      save ReSigmaR_id, ImSigmaR_id
      save H_id, S_id, En_id,Volt_id
      save nc_dim, NEn_dim, Complex_dim, Real_dim
      save En_dims,H_dims,Sigma_dims, Voltage_dim
* data variables
* attribute vectors
* enter define mode

c NetCDF: End
c================================================================
  
C================================================================
c Helpers, tempos ...
      complex*16, dimension(:), allocatable :: T
      complex*16, dimension(:), allocatable :: U, iU, ipvt
      complex*16 a,b, Ze

      real*8 max

      integer nl,nr ! size of left, right electrode regions which will be 
      integer nc                !no. : states in contact region
      integer nu
      integer icbegin,icend  !begin and end index of contact region (incl. NGL) starting index 0
      integer i,j,io,jo,ii,jj,il,jl,ir,jr

      integer icall
      logical firsttime
      data firsttime /.true./
      save firsttime,icall
      save H, Hc, Sc, T

      real*8 Pi
      parameter(Pi=3.141592653589793238462643383279503d0)
      real*8 eV
      parameter(eV=13.6058d0)


cc used for test:
      complex*16, dimension(:), allocatable :: GF, GFRGF, tt
      complex*16 csum
      real*8 tottrans

c=======================================================================
C BEGIN
c=======================================================================


c note that here we define the contact region as the states in [icend;icbegin]
c these include the first NGL contact states. 
c NOTE: We start with index 0 while isobegin/isoend assumes index 1 


      icbegin = isobegin + NGL  - 1 
      icend   = isoend   + NGL  - 1
      nc = icend - icbegin + 1
      nl = icbegin              ! the first icbegin-1 states are defined as "left electrode"
      nu = nl + nc              ! beginning of nr
      nr = nou - nu             ! the last nr states are defined as "right electrode"

      Ze = ZEnergy

c used for test
c      Ze=dcmplx(dreal(ZEnergy),0d0)



c debug      
c      write(joutfile,*) "nl,nr,nu"
c      write(joutfile,*)  nl,nr,nu
c      write(joutfile,*) "icbegin,icend,nc"
c      write(joutfile,*)  icbegin,icend,nc


      if(nc .le. 0) then
         write(joutfile) "ERROR: OutputRegionData: nc=",nc
         stop "ERROR: nc<0"
      end if

      if(firsttime) icall=0
      icall=icall+1
      if(icall.gt.NEn) 
     .     stop 'ERROR: OutputRegion: NEn data points expected'


cccccccccccccccccccccccccccccccccc
      if(firsttime) then
cccccccccccccccccccccccccccccccccc
         allocate(H(0:nou*nou-1)) ! To be Löwdin transformed
         allocate(T(0:nou*nou-1)) 
         do i=0,nou*nou-1
            H(i)=Hi(i)
            T(i)=Si(i)
         end do


c calculate T=S**(-1/2) for Löwdin transformation and transform H
         if(ORTHO) call Lowdin(.True.,nou, H, T)


ccc pick out interacting region
         allocate(Hc(0:nc*nc-1))
         allocate(Sc(0:nc*nc-1))

         do j=0,nc-1
            jo=j + icbegin
            do i=0,nc-1
               io=i + icbegin
               Hc(i+nc*j) = H(io+nou*jo)
               Sc(i+nc*j) = Si(io+nou*jo)
            end do              !i
         end do                 !j

c check for coupling between nl and nr states:

         if(ORTHO) then
            max=0.
            do jo=0,nl-1
               do io=nou-nr,nou-1
                  a = H(io+nou*jo)
                  if(cdabs(a).gt.max) max=cdabs(a)
               end do           !i
            end do              !j
            if(max.gt.0.0d0) write(joutfile,'(a50,F7.4)') 
     .      "WARNING: ORTHOGONALIZATION: 
     .           L-R coupling: max(eV)=",max*eV
         endif
C     max=0.
C     do jo=0,nl-1
C     do io=nou-nr,nou-1
C     a = Si(io+nou*jo)
C     if(cdabs(a).gt.max) max=cdabs(a)
C     end do              !i
C     end do                 !j
C     if(max.gt.0.0d0) write(joutfile,*) 
C     .        "Warning L-R electrode overlap found: max=",max
            

cccccccccccccccccccccc
c netCDF Begin
cccccccccccccccccccccc

         if(nspin.gt.1) then
            if(ispin.eq.1) 
     . iret = nf90_create('HSSigmaLR.UP.nc', NF90_CLOBBER, ncid)
            if(ispin.eq.2) 
     . iret = nf90_create('HSSigmaLR.DOWN.nc', NF90_CLOBBER, ncid)
         else
           iret = nf90_create('HSSigmaLR.nc', NF90_CLOBBER, ncid)
         end if
         call check_err(iret)
*     define dimensions
         iret = nf90_def_dim(ncid, 'NEn', NF90_UNLIMITED, NEn_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'nc', nc, nc_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'Complex', 2, Complex_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'Real', 1, Real_dim)
         call check_err(iret)

*     define variables
         En_dims(1) = Complex_dim         ! complex number repr. as two reals
         En_dims(2) = NEn_dim

         H_dims(1) = nc_dim
         H_dims(2) = nc_dim

         Sigma_dims(1) = nc_dim
         Sigma_dims(2) = nc_dim
         Sigma_dims(3) = NEn_dim
    
         Voltage_dim = Real_dim

         iret = nf90_def_var(ncid, 'En', NF90_DOUBLE, 
     .        En_dims, En_id)         
         call check_err(iret)

         iret = nf90_def_var(ncid, 'Voltage', NF90_DOUBLE, 
     .        Real_dim, Volt_id)         
         call check_err(iret)

         iret = nf90_def_var(ncid, 'H', NF90_DOUBLE,
     .        H_dims, H_id)
         call check_err(iret)

         if(.not. ORTHO) then
            iret = nf90_def_var(ncid, 'S', NF90_DOUBLE,
     .           H_dims, S_id)
            call check_err(iret)
         endif

         iret = nf90_def_var(ncid, 'ReSigmaL', NF90_DOUBLE,
     .        Sigma_dims, ReSigmaL_id)
         call check_err(iret)

         iret = nf90_def_var(ncid, 'ImSigmaL', NF90_DOUBLE,
     .        Sigma_dims, ImSigmaL_id)
         call check_err(iret)
         

         iret = nf90_def_var(ncid, 'ReSigmaR', NF90_DOUBLE,
     .        Sigma_dims, ReSigmaR_id)
         call check_err(iret)

         iret = nf90_def_var(ncid, 'ImSigmaR', NF90_DOUBLE,
     .        Sigma_dims, ImSigmaR_id)
         call check_err(iret)


*     assign attributes
         iret = nf90_put_att(ncid, En_id, 'long_name',
     .        'Energy Points (Re,Im) in Ryd.')
         call check_err(iret)

         iret = nf90_put_att(ncid, Volt_id, 'long_name',
     .        'Voltage in volts')
         call check_err(iret)


         iret = nf90_put_att(ncid, H_id, 'long_name', 
     .        'Hamiltonian in Ryd.')
         call check_err(iret)
         if(.not.ORTHO) then
            iret = nf90_put_att(ncid, S_id, 'long_name', 
     .           'Overlap matrix')
            call check_err(iret)
         endif

         iret = nf90_put_att(ncid, ReSigmaL_id, 'long_name', 
     .        'Real part of left Sigma in Ryd.')
         call check_err(iret)
         iret = nf90_put_att(ncid, ImSigmaL_id, 'long_name', 
     .        'Imaginary part of left Sigma in Ryd.')
         call check_err(iret)
         iret = nf90_put_att(ncid, ReSigmaR_id, 'long_name', 
     .        'Real part of right Sigma in Ryd.')
         call check_err(iret)
         iret = nf90_put_att(ncid, ImSigmaR_id, 'long_name', 
     .        'Imaginary part of right Sigma in Ryd.')
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 
     .        'source', 'TBTrans')
         call check_err(iret)
*     leave define mode
         iret = nf90_enddef(ncid)
         call check_err(iret)

c write Hc to netCDF:
      allocate(RMat(nc,nc))
      i=0
      do jo = 1,nc
         do io = 1,nc
            RMat(io,jo) = DReal(Hc(i))
            i=i+1
         end do
      end do
      iret = nf90_put_var(ncid, H_id, 
     .     RMat(1:,1:),start=(/1,1/))
      call check_err(iret)
      deallocate(RMat)


      if(.not.ORTHO) then
c     write Sc to netCDF:
         allocate(RMat(nc,nc))
         i=0
         do jo = 1,nc
            do io = 1,nc
               RMat(io,jo) = DReal(Sc(i))
               i=i+1
            end do
         end do
         iret = nf90_put_var(ncid, S_id, 
     .        RMat(1:,1:),start=(/1,1/))
         call check_err(iret)
         deallocate(RMat)
      endif

      iret = nf90_put_var(ncid,Volt_id,Voltage)
      call check_err(iret)


ccccccccccccccccccccccccc
c netCDF End 
ccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccc
         firsttime=.false.
      endif                     !firsttime
ccccccccccccccccccccccccccccccccccccccccccccc


      allocate(invG(0:nou*nou-1))

c     invG without self-energies:
      do i = 0, nou*nou - 1
         invG(i) = ZEnergy*Si(i)-Hi(i) 
      end do                    !i
      
c     Add Left Self-energy
      do j = 0, NGL-1
         jj=j+1
         do i = 0, NGL-1
            ii=i+1
            if(USEBULK) then 
               invG(i + j*nou) = HLL(ii,jj)
            else
               invG(i + j*nou) = invG(i + j*nou) - HLL(ii,jj)
            end if              !USEBULK
         end do                 !i
      end do                    !j

c     Add Right Self-energy
      do j = nou-NGR,nou-1
         jj=j-(nou-NGR)+1
         do i = nou-NGR,nou-1
            ii=i-(nou-NGR)+1
            if(USEBULK) then 
               invG(i + j*nou) = HRR(ii,jj)
            else
               invG(i + j*nou) = invG(i + j*nou) - HRR(ii,jj)
            end if              !USEBULK
         end do                 !i
      end do                    !j
      

c Do Löwdin transform (T is already known and stored)
c      call writezmat(nou,'nonG',invG)

      if(ORTHO) call Lowdin(.False.,nou, invG, T)

c      call writezmat(nou,'invG',invG)
c      stop 'hola'
      allocate(ipvt(nou))
      allocate(U(0:nou*nou-1))
      allocate(iU(0:nou*nou-1))

c-------------------------------------------------
c
c     SigmaL: 
c
      if(nl .lt. NGL) then
         write(joutfile,*) 
     .        "OutputRegionData: ERROR nl < NGL:", nl, "<", NGL
         stop "OutputRegionData: ERROR nl < NGL:"
      end if
      
      do j = 0, nl-1
         do i = 0, nl-1
            U(i + j*nl) = invG(i+j*nou)
         end do                 !i
      end do                    !j
                      
      do j = 0, nc-1
         do i = 0, nl-1
            jj=j + nl
            iU(i+j*nl) = -H(i+jj*nou)
         end do                 !i
      end do                    !j
      
      call csolveg(nl,nc,U,iU,ipvt,joutfile)
      
      do j = 0, nl-1
         do i = 0, nc-1
            ii=i+nl
            U(i + j*nc) = -H(ii+j*nou)
         end do                 !i
      end do                    !j

      a=(0.d0,0.d0)
      b=(1.d0,0.d0)
      
      allocate(sigmaL(0:nc*nc-1))
      
      call zgemm('N','N',nc,nc,nl,b,U,nc,iU,
     &     nl,a,sigmaL,nc)      


c-------------------------------------------------
c
c     SigmaR: 
c
      if(nr .lt. NGR) then
         write(joutfile,*) 
     .        "OutputRegionData: ERROR nr < NGR:", nr, "<", NGR
         stop "OutputRegionData: ERROR nr < NGR:"
      end if

      do j = 0, nr-1
         jr = j + nu
         do i = 0, nr-1
            ir = i + nu
            U(i + j*nr) = invG(ir + jr*nou)
         end do                 !i
      end do                    !j
      
      do j = 0, nc-1
         jj=j+nl
         do i = 0, nr-1
            ii = i + nu
            iU(i+j*nr) = -H(ii+jj*nou)
         end do                 !i
      end do                    !j

      call csolveg(nr,nc,U,iU,ipvt,joutfile)

      do j = 0, nr-1
         jj = j + nu
         do i = 0, nc-1
            ii = i + nl
            U(i + j*nc) = -H(ii+jj*nou)
         end do                 !i
      end do                    !j

      a=(0.d0,0.d0)
      b=(1.d0,0.d0)
      
      allocate(sigmaR(0:nc*nc-1))
      
      call zgemm('N','N',nc,nc,nr,b,U,nc,iU,
     &     nr,a,sigmaR,nc)      


c--------------------------------------------------------------------
c now we have Löwdin transformed Hc, sigmaR, sigmaL for the c-region
c--------------------------------------------------------------------

c write Energy to netCDF:
      Energy(1)=DReal(ZEnergy)
      Energy(2)=DImag(ZEnergy)

      iret = nf90_put_var(ncid,En_id,Energy(1:),start=(/1,icall/))
      call check_err(iret)


c write Sigmas to netCDF:
      allocate(RMat(nc,nc))
c re(SigmaL)
      i=0
      do jo = 1,nc
         do io = 1,nc
            RMat(io,jo) = DReal(SigmaL(i))
            i=i+1
         end do
      end do
      iret = nf90_put_var(ncid, ReSigmaL_id, 
     .     RMat(1:,1:),start=(/1,1,icall/))
      call check_err(iret)

c im(SigmaL)
      i=0
      do jo = 1,nc
         do io = 1,nc
            RMat(io,jo) = DImag(SigmaL(i))
            i=i+1
         end do
      end do
      iret = nf90_put_var(ncid, ImSigmaL_id, 
     .     RMat(1:,1:),start=(/1,1,icall/))
      call check_err(iret)

c re(SigmaR)
      i=0
      do jo = 1,nc
         do io = 1,nc
            RMat(io,jo) = DReal(SigmaR(i))
            i=i+1
         end do
      end do
      iret = nf90_put_var(ncid, ReSigmaR_id, 
     .     RMat(1:,1:),start=(/1,1,icall/))
      call check_err(iret)

c im(SigmaR)
      i=0
      do jo = 1,nc
         do io = 1,nc
            RMat(io,jo) = DImag(SigmaR(i))
            i=i+1
         end do
      end do
      iret = nf90_put_var(ncid, ImSigmaR_id, 
     .     RMat(1:,1:),start=(/1,1,icall/))
      call check_err(iret)

      deallocate(Rmat)
c done writing to netCDF for this energy point

      deallocate(invG)

c =============================================================
c *** TRANSMISSION TEST ...
c =============================================================
c tt  =  GammaL.(GF*.GammaR.GF):
      allocate(tt(0:nc*nc-1))
      allocate(GF(0:nc*nc-1))
      allocate(GFRGF(0:nc*nc-1))
      allocate(invG(0:nc*nc-1))

      do i=0,nc*nc-1
         invG(i) = -(Hc(i) + sigmaL(i) + sigmaR(i))
         GF(i)   = dcmplx(0d0,0d0)
      end do

      do j = 0, nc-1
         invG(j+nc*j) = invG(j+nc*j) + Ze 
         GF(j+nc*j) = dcmplx(1d0,0d0)
      end do                    !j

      call csolve(nc,invG,GF,ipvt,joutfile)


c -----------       
c GammaL -> [U], GammaR -> [iU]
      do i=0,nc*nc-1
         U(i) =  -SigmaL(i)+dconjg(SigmaL(i))
         iU(i) = -SigmaR(i)+dconjg(SigmaR(i))
      end do


c     -GammaR.GF == -[iU].GF -> [Hc]
      CALL zgemm('N','N',nc,nc,nc,dcmplx(-1d0,0d0),
     &     iU,nc,GF,nc,dcmplx(0d0,0d0),invG,nc) 
      
c     GF*.GammaR.GF == GF*.[invG] -> [iU]
      CALL zgemm('C','N',nc,nc,nc,dcmplx(1d0,0d0),
     &     GF,nc,invG,nc,dcmplx(0d0,0d0),iU,nc) 


c     GammaL.GF*.GammaR.GF == [U].[iU] -> tt
      CALL zgemm('N','N',nc,nc,nc,dcmplx(1d0,0d0),
     &     U,nc,iU,nc,dcmplx(0d0,0d0),tt,nc)


      do i=0,nc*nc - 1
         GFRGF(i) = dcmplx(0.5d0,0d0)*iU(i)
      end do

c deallocate
      deallocate(U)
      deallocate(iU)
      deallocate(sigmaL)
      deallocate(sigmaR)
      deallocate(invG)
      deallocate(GF)
      deallocate(GFRGF)


c ==============================================================

c Trace of tt:
      csum=dcmplx(0d0,0d0)
      do i=0,nc-1
         csum = csum + tt(i+nc*i)
       end do
       deallocate(tt)
       
       tottrans = DREAL(csum)
c       open(56,position='append')
c       write(66,'(F9.5,e16.8)') DREAL(Ze)*eV,tottrans
c       write(joutfile,'(a15,F9.5,e16.8)') 'OutputRegion: ',
c     .      DREAL(Ze)*eV,tottrans
c       close(56)

C ====================================================================
       deallocate(ipvt)
       if(icall == NEn) then
          write(joutfile,*) 
     .         "OutputRegionData: Done writing netcdf file: Closing"
          iret = nf90_close(ncid)
c          call check_err(iret)
          return
       end if
      

 
      CONTAINS
       
      

      subroutine check_err(status)
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'netcdf Error: Stopped'
      endif
      end subroutine check_err



c      RETURN
      END subroutine OutputRegion
C ====================================================================

#endif

