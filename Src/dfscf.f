      subroutine dfscf( no, nuo, nuotot, indxuo, ifa, istr, np, 
     .                  maxnd, numd, listdptr, listd, Dscf, Datm,
     .                  dvol, VolCel, Vscf, Vatm, Fscf, Fatm, 
     .                  iaorb, indxua, iphorb, isa, nspin, 
     .                  Fal, Stressl )

C ********************************************************************
C Adds the SCF contribution of a single orbital to atomic forces or
C stress.
C a single orbital.
C Written by P.Ordejon and J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer nuo             : Number of orbitals in unit cell (local)
C integer nuotot          : Number of orbitals in unit cell (global)
C integer indxuo(no)      : Index of equivalent atom in unit cell
C integer np              : Number of columns in C (local)
C integer maxnd           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of row in listd
C integer listd(maxnd)    : List of nonzero elements in each row of Dscf
C real*8  Dscf(maxnd)     : Value of nonzero elemens in each row of Dscf
C real*8  Datm(nuotot)    : Occupations of basis orbitals in free atom
C real*8  dvol            : Volume per mesh point
C real*8  VolCel          : Unit cell volume
C real*4  Vscf(np)        : Value of SCF potential at the mesh points
C real*4  Vatm(np)        : Value of Harris potential (Hartree potential
C                           of sum of atomic desities) at mesh points
C integer MaxDi           : Total size of auxiliary array Di
C *********************** OUTPUT **************************************
C real*8  Fscf(NX)        : Self consistent density contribution, to
C                           atomic force or stress, from orbital iOrb
C real*8  Fatm(NX)        : Free-atom density contribution
C *********************************************************************

C
C  Modules
C
      use precision
      use atmfuncs, only: rcut, phiatm, all_phi, nsmax=>nspecies
      use listsc_module, only: listsc
      use mesh, only: dxa, nsp, xdop, xdsp
      use meshphi, only: endpht, lstpht, listp2
      use meshdscf

      implicit none

C
C  Passed arguments
C
      integer
     .   no, nuo, nuotot, np, maxnd, ifa, istr, indxuo(no),
     .   iaorb(*), indxua(*), iphorb(*), nspin,
     .   isa(*), numd(nuo), listdptr(nuo), listd(maxnd)

      real
     .   Vscf(nsp,np,nspin), Vatm(nsp,np)

      real*8
     .   Dscf(maxnd,nspin), Datm(nuotot), dvol,
     .   VolCel, Fatm(3), Fscf(3), Fal(3,*), Stressl(9)

C Internal variables
      integer i, ia, imp, in, ind, iop, ip, iphi, is, isp, iu, iua, ix,
     .   ic, ispin, iy, ii, iul, j, jc, jmp, lasta, lastop, maxloc, 
     .   maxndl, nphiloc
      real*8
     .   DiiCiV, DijCj, r2o, r2sp, C, dxsp(3,nsp), 
     .   r2cut(nsmax), rvol, phia(100,nsp), FscfX(9), FatmX(9)
      real*8, dimension(:), allocatable, save ::
     .   Di
      real*8, dimension(:,:), allocatable, save ::
     .   VgrCi, VgrCiX
      real*8, dimension(:,:,:), allocatable, save ::
     .   grada, gradCi, xgradCi
      real, dimension(:,:), allocatable, save ::
     .   Clocal
      logical
     .   Parallel
      
      call timer('dfscf',1)

C  Set logical that determines whether we need to use parallel or serial mode
      Parallel = (nuo .ne. nuotot)

C  If parallel, allocate temporary storage for Local Dscf
      if (Parallel) then
        maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        allocate(DscfL(maxndl,nspin))
        call memory('A','D',maxndl*nspin,'meshdscf')

C Redistribute Dscf to DscfL form
        do ispin = 1,nspin
          call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, nuotot,
     .                     Dscf(1,ispin), DscfL(1,ispin) )
        enddo

      endif

C  Find value of maxloc
      maxloc = 0
      do ip = 1,np
        maxloc = max(maxloc,endpht(ip)-endpht(ip-1))
      enddo

C  Allocate local memory
      allocate(Di(no))
      call memory('A','D',no,'dfscf')
      allocate(Clocal(nsp,maxloc))
      call memory('A','S',nsp*maxloc,'dfscf')
      allocate(VgrCi(3,nsp))
      call memory('A','D',3*nsp,'dfscf')
      allocate(VgrCiX(9,nsp))
      call memory('A','D',9*nsp,'dfscf')
      allocate(grada(3,100,nsp))
      call memory('A','D',3*100*nsp,'dfscf')
      allocate(gradCi(3,nsp,maxloc))
      call memory('A','D',3*nsp*maxloc,'dfscf')
      allocate(xgradCi(9,nsp,maxloc))
      call memory('A','D',9*nsp*maxloc,'dfscf')

C  Initialise variables
      Di(1:no) = 0.0d0

C  Find atomic cutoff radiae
      do i = 1,no
        ia = iaorb(i)
        is = isa(ia)
        r2cut(is) = rcut(is,0)**2
      enddo

C  Evaluate constants
      rvol = 1.0d0 / VolCel

C  Loop over grid points
      do ip = 1,np

C  Calculate phi values and derivatives
        ic = 0
        lasta = 0
        lastop = 0
        do imp = 1+endpht(ip-1), endpht(ip)
          ic = ic + 1
          i = lstpht(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iphi = iphorb(i)
          is = isa(ia)
          iua = indxua(ia)
          iop = listp2(imp)

C  Calculate derivatives of phi for orbital imp and point ip

C  Loop over subpoints
          if (ia.ne.lasta .or. iop.ne.lastop) then
            lasta = ia
            lastop = iop
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp = dxsp(1,isp)**2 + dxsp(2,isp)**2 + dxsp(3,isp)**2
              if (r2sp.lt.r2cut(is)) then
                call all_phi(is,+1,dxsp(:,isp),nphiloc,phia(:,isp),
     .            grada(:,:,isp))
              else
                phia(:,isp) = 0.0d0
                grada(1:3,:,isp) = 0.0d0
              endif
            enddo
          endif
          Clocal(:,ic) = phia(iphi,:)
          gradCi(1:3,:,ic) = grada(1:3,iphi,:)

          if (istr.eq.1) then
            do isp = 1,nsp
C  If stresses are to be calculated then generate stress derivatives
              ii = 0
              do ix = 1,3
                do iy = 1,3
                  ii = ii + 1
                  xgradCi(ii,isp,ic)=dxsp(iy,isp)*gradCi(ix,isp,ic)*rvol
                enddo
              enddo
            enddo
          endif

        enddo

C  Loop on first orbital of mesh point
        ic = 0
        do imp = 1+endpht(ip-1), endpht(ip)
          ic = ic + 1
          i = lstpht(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iua = indxua(ia)

C  Loop over spins
          do ispin = 1,nspin

            if (ifa.eq.1) then
C  Initialize Fatm
              Fatm(1:3) = 0.0d0
              Fscf(1:3) = 0.0d0
            endif

            if (istr.eq.1) then
C  Initialize Fatm
              FatmX(1:9) = 0.0d0
              FscfX(1:9) = 0.0d0
            endif

C  Copy full row i of density matrix to Di(j)
            if (Parallel) then
              iul = NeedDscfL(iu)
              if (i.eq.iu) then
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listdl(ind)
                  Di(j) = DscfL(ind,ispin)
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listsc( i, iu, listdl(ind) )
                  Di(j) = DscfL(ind,ispin)
                enddo
              endif
            else
              if (i.eq.iu) then
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listd(ind)
                  Di(j) = Dscf(ind,ispin)
                enddo
              else
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listsc( i, iu, listd(ind) )
                  Di(j) = Dscf(ind,ispin)
                enddo
              endif
            endif

C  Calculate one orbital derivatives
            if (ifa.eq.1.and.istr.eq.1) then
              do isp = 1, nsp
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                Fatm(1:3) = Fatm(1:3) + DiiCiV * gradCi(1:3,isp,ic)
                VgrCi(1:3,isp) = Vscf(isp,ip,ispin)*gradCi(1:3,isp,ic)
                FatmX(1:9) = FatmX(1:9) + DiiCiV * xgradCi(1:9,isp,ic)
                VgrCiX(1:9,isp) = Vscf(isp,ip,ispin)*xgradCi(1:9,isp,ic)
              enddo
            elseif (ifa.eq.1) then
              do isp = 1, nsp
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                Fatm(1:3) = Fatm(1:3) + DiiCiV * gradCi(1:3,isp,ic)
                VgrCi(1:3,isp) = Vscf(isp,ip,ispin)*gradCi(1:3,isp,ic)
              enddo
            else
              do isp = 1, nsp
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                FatmX(1:9) = FatmX(1:9) + DiiCiV * xgradCi(1:9,isp,ic)
                VgrCiX(1:9,isp) = Vscf(isp,ip,ispin)*xgradCi(1:9,isp,ic)
              enddo
            endif

C  Add contribution from Fatm to derivatives of appropriate atom
            if (ispin.le.2) then
              if (ifa.eq.1) then
                Fal(1:3,iua) = Fal(1:3,iua) - Fatm(1:3)
              endif
              if (istr.eq.1) then
                Stressl(1:9) = Stressl(1:9) - FatmX(1:9)
              endif
            endif

C  Loop on second orbital of mesh point
            jc = 0
            do jmp = 1+endpht(ip-1), endpht(ip)
              jc = jc + 1
              j = lstpht(jmp)

C  Loop over sub-points and add 2*Dscf_ij*<Cj|Vscf|gradCi>
C  The following loops are not done using f90 form as this
C  leads to much slower execution on machines with stupid f90
C  compilers at the moment
              if (ifa.eq.1.and.istr.eq.1) then
                do isp = 1, nsp
                  DijCj = 2.d0 * dVol * Di(j) * Clocal(isp,jc)
                  do ix = 1,3
                    Fscf(ix) = Fscf(ix) + DijCj * VgrCi(ix,isp)
                  enddo
                  do ix = 1,9
                    FscfX(ix) = FscfX(ix) + DijCj * VgrCiX(ix,isp)
                  enddo
                enddo
              elseif (ifa.eq.1) then
                do isp = 1, nsp
                  DijCj = 2.d0 * dVol * Di(j) * Clocal(isp,jc)
                  do ix = 1,3
                    Fscf(ix) = Fscf(ix) + DijCj * VgrCi(ix,isp)
                  enddo
                enddo
              else
                do isp = 1, nsp
                  DijCj = 2.d0 * dVol * Di(j) * Clocal(isp,jc)
                  do ix = 1,9
                    FscfX(ix) = FscfX(ix) + DijCj * VgrCiX(ix,isp)
                  enddo
                enddo
              endif

            enddo

C  Add contribution from Fscf to derivatives of appropriate atom
            if (ispin.le.2) then
              if (ifa.eq.1) then
                Fal(1:3,iua) = Fal(1:3,iua) + Fscf(1:3)
              endif
              if (istr.eq.1) then
                Stressl(1:9) = Stressl(1:9) + FscfX(1:9)
              endif
            else
              if (ifa.eq.1) then
                Fal(1:3,iua) = Fal(1:3,iua) + 2.0d0*Fscf(1:3)
              endif
              if (istr.eq.1) then
                Stressl(1:9) = Stressl(1:9) + 2.0d0*FscfX(1:9)
              endif
            endif

C  End of spin loop
          enddo

C  Restore Di for next orbital i
          if (Parallel) then
            if (i.eq.iu) then
              do ii = 1, numdl(iul)
                j = listdl(listdlptr(iul)+ii)
                Di(j) = 0.0d0
              enddo
            else
              do ii = 1, numdl(iul)
                j = listsc( i, iu, listdl(listdlptr(iul)+ii) )
                Di(j) = 0.0d0
              enddo
            endif
          else
            if (i.eq.iu) then
              do ii = 1, numd(iu)
                j = listd(listdptr(iu)+ii)
                Di(j) = 0.0d0
              enddo
            else
              do ii = 1, numd(iu)
                j = listsc( i, iu, listd(listdptr(iu)+ii) )
                Di(j) = 0.0d0
              enddo
            endif
          endif

C  End of orbital imp loop
        enddo

C  End of mesh point loop
      enddo
  
C  Deallocate local memory
      call memory('D','D',size(Di),'dfscf')
      deallocate(Di)
      call memory('D','S',size(Clocal),'dfscf')
      deallocate(Clocal)
      call memory('D','D',size(VgrCi),'dfscf')
      deallocate(VgrCi)
      call memory('D','D',size(VgrCiX),'dfscf')
      deallocate(VgrCiX)
      call memory('D','D',size(grada),'dfscf')
      deallocate(grada)
      call memory('D','D',size(gradCi),'dfscf')
      deallocate(gradCi)
      call memory('D','D',size(xgradCi),'dfscf')
      deallocate(xgradCi)
      if (Parallel) then
        call memory('D','D',size(DscfL),'meshdscf')
        deallocate(DscfL)
      endif

      call timer('dfscf',2)
      end

