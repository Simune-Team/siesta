! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine cgvc_zmatrix( na, xa, fa, cell, stress, volume, dxmax, 
     .                 tp, ftol, strtol, varcel, relaxd, usesavecg ) 
c ***************************************************************************
c Variable-cell conjugate-gradient geometry optimization
c
c   Energy minimization including atomic coordinates and unit cell vectors.
c   It allows an external target stress:
c              %block MD.TargetStress
c                  3.5  0.0  0.0  0.0  0.0  0.0
c              %endblock MD.TargetStress
c   corresponding to xx, yy, zz, xy, xz, yz.
c   In units of (-MD.TargetPressure)
c   Default: hydrostatic pressure: -1, -1, -1, 0, 0, 0
c
c   Based on E({xa},stress), with {xa} in fractional coor
c   The gradient of the energy given by {cfa} forces (fractional!)
c   The fractional coordinates are multiplied by the initial cell vectors
c   to get them in Bohr for dxmax and preconditioning.
c      
c Written by E. Artacho. November 1998. 
c ******************************** INPUT ************************************
c integer na            : Number of atoms in the simulation cell
c real*8 fa(3,na)       : Atomic forces
c real*8 stress(3,3)    : Stress tensor components
c real*8 volume         : unit cell volume
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
c real*8 strtol         : Maximum stress tolerance
c logical varcel        : true if variable cell optimization
c logical usesavecg     : true if we're using saved CG files
c *************************** INPUT AND OUTPUT ******************************
c real*8 xa(3,na)       : Atomic coordinates
c                         input: current step; output: next step
c real*8 cell(3,3)      : Matrix of the vectors defining the unit cell 
c                         input: current step; output: next step
c                         cell(ix,ja) is the ix-th component of ja-th vector
c real*8 stress(3,3)    : Stress tensor components
c real*8 strtol         : Maximum stress tolerance
c ******************************** OUTPUT ***********************************
c logical relaxd        : True when converged
c ***************************************************************************

C
C  Modules
C
      use precision,   only : dp
      use parallel,    only : Node
      use fdf,         only : fdf_block, fdf_physical
      use zmatrix,     only : VaryZmat, lUseZmatrix, Zmat
      use zmatrix,     only : CartesianForce_to_ZmatForce
      use zmatrix,     only : ZmatForce,ZmatForceVar
      use zmatrix,     only : iZmattoVars,ZmatType
      use zmatrix,     only : Zmat_to_Cartesian
      use zmatrix,     only : coeffA, coeffB, iNextDept
      use Zmatrix,     only : ZmatForceTolLen, ZmatForceTolAng
      use Zmatrix,     only : ZmatMaxDisplLen, ZmatMaxDisplAng
      use m_conjgr,    only : conjgr
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C Subroutine arguments:
      integer,  intent(in)   :: na
      logical,  intent(in)   :: varcel, usesavecg
      logical,  intent(out)  :: relaxd
      real(dp), intent(in)   :: fa(3,na), volume, dxmax, tp, ftol
      real(dp), intent(inout):: xa(3,na), stress(3,3), strtol, cell(3,3)

C Internal variables and arrays
      integer             :: iu, ia, i, j, n, indi,indi1,vi,k
      logical             :: found
      real(dp), parameter :: kBar = 1.d0/1.47108d5
      real(dp)            :: cellin(3,3), sxx, syy, szz, sxy, sxz, syz
      real(dp)            :: celli(3,3)
      real(dp), dimension(:), allocatable :: gxa, gfa
      real(dp)            :: force, force1
#ifdef MPI
      integer           MPIerror
#endif

C Saved internal variables:
      integer,       save :: ndeg, linmin
      logical,       save :: frstme = .true.
      logical,       save :: tarstr = .false.
      real(dp),      save :: cgcntr(0:20) = 0.0_dp
      real(dp),      save :: tstres(3,3) 
      real(dp),      save :: modcel(3) 
      real(dp),      save :: precon 
      real(dp),      save :: strain(3,3)

      real(dp), dimension(:), allocatable, save :: cgaux
      real(dp), dimension(:), allocatable, save :: ftoln, dxmaxn
      real(dp), dimension(:), allocatable, save :: ftoln_tmp


      external          memory

c ---------------------------------------------------------------------------

C Allocate local memory
      if (.not.allocated(cgaux)) then
        allocate(cgaux(2*(na*3 + 6)))
        call memory('A','D',6*na+12,'cgvc_zmatrix')
      endif
      allocate(gfa(na*3 + 6))
      call memory('A','D',3*na+6,'cgvc_zmatrix')
      allocate(gxa(na*3 + 6))
      call memory('A','D',3*na+6,'cgvc_zmatrix')
      allocate(ftoln(na*3 + 6))
      allocate(ftoln_tmp(na*3 + 6))
      call memory('A','D',3*na+6,'cgvc_zmatrix')
      allocate(dxmaxn(na*3 + 6))
      call memory('A','D',3*na+6,'cgvc_zmatrix')

C If first call to cgvc_zmatrix, check dim and get target stress --------------------

      if ( frstme ) then
  
C Look for target stress and read it if found, otherwise generate it --------

        if ( varcel ) then
     
          if (Node.eq.0) then
            tarstr = fdf_block('MD.TargetStress',iu)

            if (tarstr) then
              write(6,'(/a,a)')
     $              'cgvc_zmatrix: Reading %block MD.TargetStress',
     .                       ' (units of MD.TargetPressure).'
              read(iu,*, end=50) sxx, syy, szz, sxy, sxz, syz
              tstres(1,1) = - sxx * tp
              tstres(2,2) = - syy * tp
              tstres(3,3) = - szz * tp
              tstres(1,2) = - sxy * tp
              tstres(2,1) = - sxy * tp
              tstres(1,3) = - sxz * tp
              tstres(3,1) = - sxz * tp
              tstres(2,3) = - syz * tp
              tstres(3,2) = - syz * tp
   50         continue
            else
              write(6,'(/a,a)')
     $              'cgvc_zmatrix: No target stress found, ',
     .              'assuming hydrostatic MD.TargetPressure.'
              do i = 1, 3
                do j = 1, 3
                  tstres(i,j) = 0.0_dp
                enddo
                tstres(i,i) = - tp
              enddo
            endif

C Write target stress down --------------------------------------------------

            write(6,"(/a)") 'cgvc_zmatrix: Target stress (kBar)'
            write(6,"(a,2x,3f12.3)") 
     .        'cgvc_zmatrix:', tstres(1,1)/kBar, tstres(1,2)/kBar, 
     .        tstres(1,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .        'cgvc_zmatrix:', tstres(2,1)/kBar, tstres(2,2)/kBar, 
     .        tstres(2,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .        'cgvc_zmatrix:', tstres(3,1)/kBar, tstres(3,2)/kBar, 
     .        tstres(3,3)/kBar
          endif

#ifdef MPI
C Globalise information
          call MPI_Bcast(tstres(1,1),9,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
#endif

C Moduli of original cell vectors for fractional coor scaling back to au ---
          do n = 1, 3
            modcel(n) = 0.0_dp
            do j = 1, 3
              modcel(n) = modcel(n) + cell(j,n)*cell(j,n)
            enddo
            modcel(n) = dsqrt( modcel(n) )
          enddo

          if (Node.eq.0) then
C Scale factor for strain variables to share magnitude with coordinates -----
C ---- (a length in Bohrs typical of bond lengths ..) -----------------------

            precon = fdf_physical('MD.PreconditionVariableCell',
     .                           9.4486344d0,'Bohr')
          endif
#ifdef MPI
C Globalise information
          call MPI_Bcast(precon,1,MPI_double_precision,0,MPI_Comm_World,
     .      MPIerror)
#endif

C Initialize absolute strain and save initial cell vectors -----------------
C Initialization to 1. for numerical reasons, later substracted ------------
          do i = 1, 3
            do j = 1, 3
              strain(i,j) = 1.0_dp
              cellin(i,j) = cell(i,j)
            enddo
          enddo

C Set number of variables for non-Zmatrix form for iocg call
          ndeg = na*3 + 6

        else

C Set number of variables for non-Zmatrix form for iocg call
          ndeg = na*3

        endif
C End of variable cell commands        

        if (lUseZmatrix) then
C Find number of variables correctly if this is a Z-matrix run for iocg call
          ndeg = 0
          do ia = 1,na
            do n = 1,3
              indi = 3*(ia-1) + n
              if (VaryZmat(indi)) then
                ndeg = ndeg + 1
              endif
            enddo
          enddo
        endif

C Initialize and read cgaux and cgcntr if present and wanted ---------------
        if (usesavecg) then
          if (Node.eq.0) then
            call iocg( 'read', ndeg*2, cgaux, cgcntr, relaxd, found )
          endif
#ifdef MPI
C Globalise information
          call MPI_Bcast(cgaux,2*(na*3+6),MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(cgcntr,21,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(relaxd,1,MPI_logical,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,
     .      MPIerror)
#endif
          if ( found ) then
            linmin = cgcntr(1)
          else
            if (Node.eq.0) then
              write(6,'(/,a)')
     $              'cgvc_zmatrix: WARNING: CG file not found'
            endif
            relaxd = .false.
            cgcntr(0) = 0
            linmin = 1
          endif
        else
          relaxd = .false.
          cgcntr(0) = 0
          linmin = 1
        endif

        frstme = .false.
C End of commands that only get executed on first call of cgvc_zmatrix
      endif


C Variable cell -------------------------------------------------------------

      if ( varcel ) then

C Inverse of matrix of cell vectors  (transverse of) ------------------------

        call reclat( cell, celli, 0 )

C Transform coordinates and forces to fractional ---------------------------- 
C but scale them again to Bohr by using the (fix) moduli of the original ----
C lattice vectors (allows using maximum displacement as before) -------------
C convergence is checked here for input forces as compared with ftol --------

        relaxd = .true.
        ndeg = 0
        if (lUseZmatrix) then
C Loop over degrees of freedom, scaling cartesian coordinates to fractional
          do ia = 1,na
            do n = 1,3
              indi = 3*(ia-1) + n
              if (VaryZmat(indi)) then
                ndeg = ndeg + 1
                if (ZmatType(indi).gt.1) then
                  ftoln(ndeg) = ZmatForceTolLen
                  dxmaxn(ndeg) = ZmatMaxDisplLen
                else
                  ftoln(ndeg) = ZmatForceTolAng
                  dxmaxn(ndeg) = ZmatMaxDisplAng
                endif
                vi = iZmattoVars(indi)
                if (vi.eq.0) then
                  force = ZmatForce(indi)
                else
                  force = ZmatForceVar(vi)
                endif 
                relaxd=relaxd.and.(dabs(force).lt.ftoln(ndeg))
                if (ZmatType(indi).gt.2) then
C Cartesian coordinate                
                  gxa(ndeg) = 0.0_dp
                  gfa(ndeg) = 0.0_dp
                  do i = 1,3
                    indi1 = 3*(ia-1)+i
                    gxa(ndeg) = gxa(ndeg)+
     .                          celli(i,n)*Zmat(indi1)*modcel(n)
                    if (i.eq.n) then
                      force1 = force
                    else
                      force1 = ZmatForce(indi1)
                    endif
                    gfa(ndeg) = gfa(ndeg)+ 
     .                          cell(i,n)*force1/modcel(n)
                  enddo
                else
                  gxa(ndeg) = Zmat(indi)
                  gfa(ndeg) = force
                endif
              endif
            enddo
          enddo
        else
C No Z-Matrix
          do ia = 1,na
            do n = 1,3
              ndeg = ndeg + 1
              ftoln(ndeg) = ftol
              dxmaxn(ndeg) = dxmax
              gxa(ndeg) = 0.0_dp
              gfa(ndeg) = 0.0_dp
              relaxd = relaxd .and. ( dabs(fa(n,ia)) .lt. ftoln(ndeg) )
              do i = 1, 3
                gxa(ndeg) = gxa(ndeg)+celli(i,n)*xa(i,ia)*modcel(n)
                gfa(ndeg) = gfa(ndeg)+ cell(i,n)*fa(i,ia)/modcel(n)
              enddo
            enddo
          enddo
        endif


C Symmetrizing the stress tensor --------------------------------------------
        do i = 1,3
          do j = i+1,3
            stress(i,j) = 0.5_dp*( stress(i,j) + stress(j,i) )
            stress(j,i) = stress(i,j)
          enddo
        enddo

C Append stress (substracting target stress) and strain to gxa and gfa ------ 
C preconditioning: scale stress and strain so as to behave similarly to x,f -
        do i = 1,3
          do j = i,3
            ndeg = ndeg + 1
            gfa(ndeg) = -(stress(i,j) - tstres(i,j))*volume/precon
            gxa(ndeg) = strain(i,j) * precon
            if (lUseZmatrix) then
              dxmaxn(ndeg) = ZmatMaxDisplLen
            else
              dxmaxn(ndeg) = dxmax
            endif
          enddo
        enddo

C  testing

C Check stress convergence --------------------------------------------------
        strtol = dabs(strtol)
        do i = 1,3
          do j = 1,3
            relaxd = relaxd .and. 
     .        ( dabs(stress(i,j)-tstres(i,j)) .lt. strtol )
          enddo
        enddo

C Call conjugate gradient minimization -------------------------------------- 
        if ( .not. relaxd ) then
          do i=1,ndeg 
            ftoln_tmp(i) = 0.d0
          enddo
          call conjgr(ndeg,gxa,gfa,dxmaxn,ftoln_tmp,cgcntr,cgaux)
        endif

C Fixed cell ----------------------------------------------------------------

      else
C Set number of degrees of freedom & variables
        ndeg = 0
        if (lUseZmatrix) then
          do i = 1,na
            do k = 1,3
              indi = 3*(i-1)+k
              if (VaryZmat(indi)) then
                ndeg = ndeg + 1
                gxa(ndeg) = Zmat(indi)
                vi = iZmattoVars(indi)
                if (vi.eq.0) then
                  force = ZmatForce(indi)
                else
                  force = ZmatForceVar(vi)
                endif
                gfa(ndeg) = force
                if (ZmatType(indi).eq.1) then
                  ftoln(ndeg) = ZmatForceTolAng
                  dxmaxn(ndeg) = ZmatMaxDisplAng
                else
                  ftoln(ndeg) = ZmatForceTolLen
                  dxmaxn(ndeg) = ZmatMaxDisplLen
                endif
              endif
            enddo
          enddo
        else
          do i = 1,na
            do k = 1,3
              ndeg = ndeg + 1
              gxa(ndeg) = xa(k,i)
              gfa(ndeg) = fa(k,i)
              dxmaxn(ndeg) = dxmax
              ftoln(ndeg) = ftol
            enddo
          enddo
        endif

C Call conjugant gradient routine
        call conjgr( ndeg, gxa, gfa, dxmaxn, ftoln, cgcntr, cgaux)
        relaxd = (int(cgcntr(0)) .eq. 0)

      endif
C End of fixed cell

C Checking line minimizations and convergence -------------------------------

      if (nint(cgcntr(1)) .ne. linmin) then
        if (Node.eq.0) then
          write(6,'(/a,i4,a,f10.4)')
     .      'cgvc_zmatrix: Finished line minimization ', linmin,
     .      '.  Mean atomic displacement =', cgcntr(18)/sqrt(dble(na))
        endif
        linmin = nint(cgcntr(1))
      endif

C Transform back if variable cell ------------------------------------------- 
      if ( varcel .and. (.not. relaxd) ) then

C New cell ------------------------------------------------------------------
        indi = ndeg-6
        do i = 1,3
          do j = i,3
            indi = indi + 1
            strain(i,j) = gxa(indi) / precon
            strain(j,i) = strain(i,j)
          enddo
        enddo

        do n = 1,3
          do i = 1,3
            cell(i,n) = cellin(i,n)
            do j = 1,3
              cell(i,n) = cell(i,n) + (strain(i,j)-1.0_dp)*cellin(j,n)
            enddo
          enddo
        enddo

C Output fractional coordinates to cartesian Bohr, and copy to xa ----------- 
        ndeg = 0
        if (lUseZmatrix) then
          do ia=1,na
            do k=1,3
              indi = 3*(ia-1)+k
              if (VaryZmat(indi)) then
                ndeg = ndeg + 1
                j = indi
                do while (j.ne.0) 
                  if (ZmatType(j).gt.2) then
                    Zmat(j) = 0.0_dp
                    do n=1,3
                      indi1 = 3*(ia-1)+n
                      Zmat(j)=Zmat(j)+cell(k,n)*gxa(indi1)/modcel(n)
                    enddo
                  else
                    Zmat(j) = gxa(ndeg)
                  endif
                  Zmat(j) = Zmat(j)*coeffA(j) + coeffB(j)
                  j = iNextDept(j)
                enddo
              endif
            enddo
          enddo  
          call Zmat_to_Cartesian(xa)
        else
          do ia = 1,na
            do k = 1,3
              ndeg = ndeg+1
              xa(k,ia) = 0.0_dp
              do n = 1,3
                indi1 = 3*(ia-1) + n
                xa(k,ia) = xa(k,ia) + cell(k,n) * gxa(indi1) / modcel(n)
              enddo
            enddo
          enddo
        endif

      else
C Fixed cell
C Return coordinates to correct arrays 
        ndeg = 0
        if (lUseZmatrix) then
          do ia = 1,na
            do k = 1,3
              indi = 3*(ia-1)+k
              if (VaryZmat(indi)) then
                ndeg = ndeg + 1
                j = indi
                do while (j.ne.0)
                  Zmat(j) = gxa(ndeg)*coeffA(j)+coeffB(j)
                  j = iNextDept(j)
                enddo
              endif
            enddo
          enddo
          call Zmat_to_Cartesian(xa)
        else
          do ia = 1,na
            do k = 1,3
              ndeg = ndeg + 1
              xa(k,ia) = gxa(ndeg)
            enddo
          enddo
        endif

      endif

C Save cgaux ----------------------------------------------------------------

      if (Node.eq.0) 
     .  call iocg( 'write', ndeg*2, cgaux, cgcntr, relaxd, found )

C Deallocate local memory

      call memory('D','D',size(dxmaxn),'cgvc_zmatrix')
      deallocate(dxmaxn)
      call memory('D','D',size(ftoln_tmp),'cgvc_zmatrix')
      deallocate(ftoln_tmp)
      call memory('D','D',size(ftoln),'cgvc_zmatrix')
      deallocate(ftoln)
      call memory('D','D',size(gxa),'cgvc_zmatrix')
      deallocate(gxa)
      call memory('D','D',size(gfa),'cgvc_zmatrix')
      deallocate(gfa)

      return
      end
