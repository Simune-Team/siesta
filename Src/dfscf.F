      subroutine dfscf( NO, indxuo, 
     .                  ifa, istr, NSPmax, 
     .                  NSP, NP, NPF, endCt, listCt, CtoCt,
     .                  NDmax, numDs, listDs, Dscf, Datm,
     .                  VolCel, Vscf, Vatm, 
     .                  Fscf, Fatm, Di, 
     .                  iaorb, indxua, iphorb, isa, 
     .                  maxNO, nspin, Fal, Stressl, dxa, 
     .                  xdop, xdsp, listP2, DscfPtr )

C ********************************************************************
C Adds the SCF contribution of a single orbital to atomic forces or
C stress.
C a single orbital.
C Written by P.Ordejon and J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer indxuo(NO)      : Index of equivalent atom in unit cell
C integer NSPmax          : First dimension of array C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of columns in C (local)
C integer NPF             : Number of columns in C (global)
C integer endCt(NP)       : Accumulated number of nonzero elements in
C                           each column of C
C integer listCt(*)       : List of nonzero elements in each column of C
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C integer NDmax           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numDs(NO)       : Number of nonzero elemts in each row of Dscf
C integer listDs(NDmax,NO): List of nonzero elements in each row of Dscf
C real*8  Dscf(NDmax)     : Value of nonzero elemens in each row of Dscf
C real*8  Datm(NO)        : Occupations of basis orbitals in free atom
C real*8  VolCel          : Unit cell volume
C real*4  Vscf(NP)        : Value of SCF potential at the mesh points
C real*4  Vatm(NP)        : Value of Harris potential (Hartree potential
C                           of sum of atomic desities) at mesh points
C integer MaxDi           : Total size of auxiliary array Di
C integer DscfPtr(NO)     : Pointer from global to locally needed Dscf row
C *********************** OUTPUT **************************************
C real*8  Fscf(NX)        : Self consistent density contribution, to
C                           atomic force or stress, from orbital iOrb
C real*8  Fatm(NX)        : Free-atom density contribution
C ********************* AUXILIARY *************************************
C real*8  Di(NO)          : Array that must be initialized to zero
C                           before call for first orbital and which
C                           MUST NOT BE MODIFIED before call for last
C                           orbital. 
C *********************************************************************

C
C  Modules
C
      use precision
      use atmfuncs, only: rcut, phiatm
      use listsc_module, only: listsc
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .   NO, NP, NPF, NDmax, NSP, NSPmax, 
     .   ifa, istr,
     .   endCt(0:NP), listCt(*), CtoCt(*), 
     .   numDs(NO), listDs(NDmax,NO), indxuo(NO),
     .   iaorb(*), indxua(*), iphorb(*), maxNO, nspin,
     .   isa(*), listP2(*), DscfPtr(*)

      real
     .   Vscf(NSPmax,NP,nspin), Vatm(NSPmax,NP)

      double precision
     .   Dscf(NDmax,MaxNO,nspin), Datm(NO), VolCel,
     .   Fatm(3), Fscf(3), Di(NO), Fal(3,*), dxa(3,*),
     .   xdop(3,*), xdsp(3,*), Stressl(9)

      integer i, in, ip, isp, iu, ix, j, jn, 
     .   Node, iua, ia, iphi, is, imp,
     .   jmp, iop, ispin, iy, ii, iul
#ifdef MPI
      integer MPIerror
#endif
      integer maxloc, ic, jc
      double precision DiiCiV, DijCj, dVol, 
     .   r2o, r2sp, C, dxsp(3), rvol, phip,
     .   FscfX(9), FatmX(9)
      double precision, dimension(:,:), allocatable, save ::
     .   VgrCi, VgrCiX
      double precision, dimension(:,:,:), allocatable, save ::
     .   gradCi, xgradCi
      real, dimension(:,:), allocatable, save ::
     .   Clocal
      
      call timer('dfscf',1)

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C  Find value of maxloc
      maxloc = 0
      do ip = 1,np
        maxloc = max(maxloc,endCt(ip)-endCt(ip-1))
      enddo

C  Allocate local memory
      allocate(Clocal(NSP,maxloc))
      call memory('A','S',nsp*maxloc,'dfscf')
      allocate(VgrCi(3,NSP))
      call memory('A','D',3*nsp,'dfscf')
      allocate(VgrCiX(9,NSP))
      call memory('A','D',9*nsp,'dfscf')
      allocate(gradCi(3,NSP,maxloc))
      call memory('A','D',3*nsp*maxloc,'dfscf')
      allocate(xgradCi(9,NSP,maxloc))
      call memory('A','D',9*nsp*maxloc,'dfscf')

C  Evaluate constants
      dVol = VolCel / (NPF*NSP)
      rvol = 1.0d0 / VolCel

C  Loop over grid points
      do ip = 1,np

C  Calculate phi values and derivatives
        ic = 0
        do imp = 1+endCt(ip-1), endCt(ip)
          ic = ic + 1
          i = listCt(imp)
          in = CtoCt(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iphi = iphorb(i)
          is = isa(ia)
          r2o = rcut(is,iphi)**2
          iua = indxua(ia)
          iop = listP2(in)

C  Calculate derivatives of phi for orbital imp and point ip

C  Loop over subpoints
          do isp = 1,nsp
            do ix = 1,3
              dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
            enddo
            r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
            if (r2sp.lt.r2o) then
              call phiatm(is,iphi,dxsp,phip,gradCi(1,isp,ic))
              Clocal(isp,ic) = phip
            else
              Clocal(isp,ic) = 0.0
              do ix = 1,3
                gradCi(ix,isp,ic) = 0.0d0
              enddo
            endif

C  If stresses are to be calculated then generate stress derivatives
            if (istr.eq.1) then
              ii = 0
              do ix = 1,3
                do iy = 1,3
                  ii = ii + 1
                  xgradCi(ii,isp,ic)=dxsp(iy)*gradCi(ix,isp,ic)*rvol
                enddo
              enddo
            endif

          enddo

        enddo

C  Loop on first orbital of mesh point
        ic = 0
        do imp = 1+endCt(ip-1), endCt(ip)
          ic = ic + 1
          i = listCt(imp)
          in = CtoCt(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iua = indxua(ia)

C  Loop over spins
          do ispin = 1,nspin

            if (ifa.eq.1) then
C     Initialize Fatm
              do ix = 1,3
                Fatm(ix) = 0.d0
                Fscf(ix) = 0.d0
              enddo
            endif

            if (istr.eq.1) then
C     Initialize Fatm
              do ix = 1,9
                FatmX(ix) = 0.d0
                FscfX(ix) = 0.d0
              enddo
            endif

C  Copy full row i of density matrix to Di(j)
            iul = DscfPtr(iu)
            if (i.eq.iu) then
              do ii = 1, numDs(iu)
                j = listDs(ii,i)
                Di(j) = Dscf(ii,iul,ispin)
              enddo
            else
              do ii = 1, numDs(iu)
                j = listsc( i, iu, listDs(ii,iu) )
                Di(j) = Dscf(ii,iul,ispin)
              enddo
            endif

C  Calculate one orbital derivatives
            if (ifa.eq.1.and.istr.eq.1) then
              do isp = 1, NSP
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                do ix = 1,3
                  Fatm(ix) = Fatm(ix) + DiiCiV * gradCi(ix,isp,ic)
                  VgrCi(ix,isp) = Vscf(isp,ip,ispin)*gradCi(ix,isp,ic)
                enddo
                do ix = 1,9
                  FatmX(ix) = FatmX(ix) + DiiCiV * xgradCi(ix,isp,ic)
                  VgrCiX(ix,isp) = Vscf(isp,ip,ispin)*xgradCi(ix,isp,ic)
                enddo
              enddo
            elseif (ifa.eq.1) then
              do isp = 1, NSP
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                do ix = 1,3
                  Fatm(ix) = Fatm(ix) + DiiCiV * gradCi(ix,isp,ic)
                  VgrCi(ix,isp) = Vscf(isp,ip,ispin)*gradCi(ix,isp,ic)
                enddo
              enddo
            else
              do isp = 1, NSP
                DiiCiV = 2.d0*dVol*Datm(iu)*Clocal(isp,ic)*Vatm(isp,ip)
                do ix = 1,9
                  FatmX(ix) = FatmX(ix) + DiiCiV * xgradCi(ix,isp,ic)
                  VgrCiX(ix,isp) = Vscf(isp,ip,ispin)*xgradCi(ix,isp,ic)
                enddo
              enddo
            endif

C  Add contribution from Fatm to derivatives of appropriate atom
            if (ispin.le.2) then
              if (ifa.eq.1) then
                do ix = 1,3
                  Fal(ix,iua) = Fal(ix,iua) - Fatm(ix)
                enddo
              endif
              if (istr.eq.1) then
                do ix = 1,9
                  Stressl(ix) = Stressl(ix) - FatmX(ix)
                enddo
              endif
            endif

C  Loop on second orbital of mesh point
            jc = 0
            do jmp = 1+endCt(ip-1), endCt(ip)
              jc = jc + 1
              j = listCt(jmp)
              jn = CtoCt(jmp)

C         Loop over sub-points and add 2*Dscf_ij*<Cj|Vscf|gradCi>
              if (ifa.eq.1.and.istr.eq.1) then
                do isp = 1, NSP
                  DijCj = 2.d0 * dVol * Di(j) * Clocal(isp,jc)
                  do ix = 1,3
                    Fscf(ix) = Fscf(ix) + DijCj * VgrCi(ix,isp)
                  enddo
                  do ix = 1,9
                    FscfX(ix) = FscfX(ix) + DijCj * VgrCiX(ix,isp)
                  enddo
                enddo
              elseif (ifa.eq.1) then
                do isp = 1, NSP
                  DijCj = 2.d0 * dVol * Di(j) * Clocal(isp,jc)
                  do ix = 1,3
                    Fscf(ix) = Fscf(ix) + DijCj * VgrCi(ix,isp)
                  enddo
                enddo
              else
                do isp = 1, NSP
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
                do ix = 1,3
                  Fal(ix,iua) = Fal(ix,iua) + Fscf(ix)
                enddo
              endif
              if (istr.eq.1) then
                do ix = 1,9
                  Stressl(ix) = Stressl(ix) + FscfX(ix)
                enddo
              endif
            else
              if (ifa.eq.1) then
                do ix = 1,3
                  Fal(ix,iua) = Fal(ix,iua) + 2.0d0*Fscf(ix)
                enddo
              endif
              if (istr.eq.1) then
                do ix = 1,9
                  Stressl(ix) = Stressl(ix) + 2.0d0*FscfX(ix)
                enddo
              endif
            endif

C  End of spin loop
          enddo

C  Restore Di for next orbital i
          if (i.eq.iu) then
            do ii = 1, numDs(iu)
              j = listDs(ii,i)
              Di(j) = 0.d0
            enddo
          else
            do ii = 1, numDs(iu)
              j = listsc( i, iu, listDs(ii,iu) )
              Di(j) = 0.d0
            enddo
          endif

C  End of orbital imp loop
        enddo

C  End of mesh point loop
      enddo
  
C  Deallocate local memory
      call memory('D','S',size(Clocal),'dfscf')
      deallocate(Clocal)
      call memory('D','D',size(VgrCi),'dfscf')
      deallocate(VgrCi)
      call memory('D','D',size(VgrCiX),'dfscf')
      deallocate(VgrCiX)
      call memory('D','D',size(gradCi),'dfscf')
      deallocate(gradCi)
      call memory('D','D',size(xgradCi),'dfscf')
      deallocate(xgradCi)

      call timer('dfscf',2)
      end

