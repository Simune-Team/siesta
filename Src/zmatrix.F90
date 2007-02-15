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
!**************************************************************************
!  This module contains subroutines required for Z-matrix based geometry  *
!  optimisation in SIESTA.                                                *
!                                                                         *
!  Written by Rainer Hoft (UTS) & Julian Gale (Curtin), March - May 2005  *
!**************************************************************************
MODULE zmatrix

  use precision, only : dp
  use units, only : Ang, eV, deg
  use sys, only: die

  implicit none

  !  Variables

  !  iZmat holds the pointers to the atoms to which the components
  !  of Zmat relate
  integer,  pointer, save :: iZmat(:)
  !  VaryZmat(3*na) indicates whether the component of Zmat is fixed,
  logical,  pointer, save :: VaryZmat(:)
  !  Zmat(3*na) holds the Z-matrix components for each atom
  real(dp), pointer, save :: Zmat(:)
  !  ZmatForce(3*na) holds the forces on each Zmat coordinate
  real(dp), pointer, save :: ZmatForce(:)
  !  ZmatForceVar(nVar) holds the forces on additional constrained coordinates
  real(dp), pointer, save :: ZmatForceVar(:)
  !  lUseZmatrix indicates whether Zmatrix approach is to be used
  logical,           save :: lUseZmatrix = .false.

  !  nZmol is the number of molecules with Z matrices defined
  integer,           save :: nZmol
  !  nZmolAtoms(nZmol) indicates the number of atoms in each molecule
  integer,  pointer, save :: nZmolAtoms(:)
  !  nZmolStartAtom(nZmol) is a pointer to the first atom of each molecule
  integer,  pointer, save :: nZmolStartAtom(:)
  !  nZcart is the number of cartesian blocks defined in the Zmatrix
  integer,           save :: nZcart
  !  nZcartAtoms(nZcart) indicates the number of atoms in each cartesian block
  integer,  pointer, save :: nZcartAtoms(:)
  !  nZcartStartAtom(nZcart) is a pointer to the first atom of each cartesian block
  integer,  pointer, save :: nZcartStartAtom(:)
  !  Variables to specify the input (and output) units of the Zmatrix
  integer,           save :: ZmatUnitsLen
  integer,           save :: ZmatUnitsAng
  !  Variables to specify the maximum force tolerance for CG optimization
  real(dp),          save :: ZmatForceTolLen
  real(dp),          save :: ZmatForceTolAng
  !  Variables to specify the maximum displacement during CG optimization
  real(dp),          save :: ZmatMaxDisplLen
  real(dp),          save :: ZmatMaxDisplAng
  !  ZmatType(3*na) defines the component type:
  !    1 - angle
  !    2 - bond length
  !    3 - pure cartesian (molecule)
  !    4 - scaled cartesian (molecule)
  !    5 - fractional cartesian (molecule)
  !    6 - pure cartesian
  !    7 - scaled cartesian (scaled by lattice constant)
  !    8 - fractional cartesian (scaled by lattice vectors)
  integer,  pointer, save :: ZmatType(:)
  !  nVars specifies the number of constrained variables
  integer,           save :: nVars
  !  ZmatVarNames(nVars) gives the names of the constrained variables
  character(len=10), allocatable, save :: ZmatVarNames(:)
  !  iZmattoVars(3*na) is the index in the variable array of a symbolic coordinate
  integer, pointer,  save :: iZmattoVars(:) 
  !  iVarstoZmat(nVars) is the Zmat index of the first symbolic coordinate corresponding to the variable
  integer, pointer,  save :: iVarstoZmat(:)
  !  lZmatVarsDef(nVars) tells us whether the variable name has been given a value
  logical, pointer,  save :: lZmatVarsDef(:)
  !  lCalcAllForces specifies whether forces for fixed coorinates should be calculated
  logical,           save :: lCalcAllForces
  !  iNextDept(3*na) gives the linked list of dependencies on symbolic coordinates (note it is set up so that each coordinate can only depend on one other)
  integer, pointer,  save :: iNextDept(:)
  !  coeffA and coeffB contains the linear relationships between variables v2 = a*v1 + B
  real(dp), pointer, save :: coeffA(:)
  real(dp), pointer, save :: coeffB(:)

  real(dp), save  :: scale_length, scale_angle 
  real(dp), save  :: zmatrix_alat
  ! AG: convenience variables

  public :: read_Zmatrix, lUseZmatrix, iofaZmat
  public :: CartesianForce_to_ZmatForce, write_Zmatrix
  public :: VaryZmat, Zmat, ZmatForce, ZmatForceVar
  public :: iZmattoVars, ZmatType, Zmat_to_Cartesian
  public :: coeffA, coeffB, iNextDept
  public :: ZmatForceTolLen, ZmatForceTolAng
  public :: ZmatMaxDisplLen, ZmatMaxDisplAng
  private 

  CONTAINS

  subroutine read_Zmatrix(na,nSpecies,alat,ucell,lOrigin,origin)

   use fdf
   use parsing,     only : parse
   use alloc,       only : re_alloc
   use sys,         only : die
   use parallel,    only : Node
#  ifdef MPI
   use mpi_siesta
#  endif

   ! Passed variables
   integer,  intent(in)          :: na
   integer,  intent(out)         :: nSpecies(na)
   real(dp), intent(in)          :: alat
   real(dp), intent(in)          :: ucell(3,3)
   logical,  intent(in)          :: lOrigin
   real(dp), intent(in)          :: origin(3)

   ! Local variables
   logical                       :: leqi
   character(len=130)            :: line
   character(len= 80)            :: names
   character(len= 10)            :: angleunits
   character(len= 10)            :: lengthunits
   integer                       :: i
   integer                       :: j
   integer                       :: integs(10)
   integer                       :: iu
   integer                       :: lastc
   integer                       :: lc(0:3)
   integer                       :: m
   integer                       :: k
   integer                       :: nStart
   integer                       :: nAtoms
   integer                       :: ni
   integer                       :: nn
   integer                       :: nr
   integer                       :: nv
   logical                       :: eof
   logical,                 save :: firsttime = .true.
   real(dp)                      :: reals(10)
   real(dp)                      :: values(10)
   character(len=100)            :: errormsg
   integer                       :: input_type
   integer                       :: units_type
   character(len=10), parameter  :: unitslen_default = 'Bohr'
   character(len=10), parameter  :: unitsang_default = 'rad' 
   real(dp),          parameter  :: ftollen_default = 1.55574d-3
   real(dp),          parameter  :: ftolang_default = 3.56549d-3
   real(dp),          parameter  :: dxmaxlen_default = 0.2d0
   real(dp),          parameter  :: dxmaxang_default = 0.003d0
   logical,           parameter  :: lCalcAllForces_default=.false.
   integer                       :: order(10)
   logical                       :: found
   integer                       :: nnames
   integer                       :: nreals
   real(dp)                      :: ZmatVal(3)
   real(dp)                      :: coeffBVal(3)
   integer                       :: indi
   integer                       :: indi1
   integer                       :: flagnr
   real(dp)                      :: A
   real(dp)                      :: B
#  ifdef MPI
   integer                       :: integers(2),integers1(4)
   integer                       :: MPIerror
   real(dp)                      :: physicals(4)
#  endif

   zmatrix_alat = alat     ! AG: For future processing

    ! Nullify pointers
    if (firsttime) then
      nullify(iZmat,nZmolAtoms,nZmolStartAtom)
      nullify(nZcartAtoms,nZcartStartAtom)
      nullify(Zmat,VaryZmat,ZmatForce,ZmatForceVar)
      nullify(ZmatType,iZmattoVars,iVarstoZmat)
      nullify(lZmatVarsDef,iNextDept)
      nullify(coeffA,coeffB)

      ! Allocate Zmatrix arrays
      call re_alloc(Zmat,1,3*na)
      call re_alloc(ZmatForce,1,3*na)
      call re_alloc(iZmat,1,3*na)
      call re_alloc(VaryZmat,1,3*na)
      call re_alloc(nZmolAtoms,1,na)
      call re_alloc(nZmolStartAtom,1,na)
      call re_alloc(nZcartAtoms,1,na)
      call re_alloc(nZcartStartAtom,1,na)
      call re_alloc(ZmatType,1,3*na)
      call re_alloc(iZmattoVars,1,3*na)
      call re_alloc(iNextDept,1,3*na)
      call re_alloc(coeffA,1,3*na)

      call re_alloc(coeffB,1,3*na)
      call re_alloc(ZmatForceVar,1,3*na)
      call re_alloc(iVarstoZmat,1,3*na)
      call re_alloc(lZmatVarsDef,1,3*na)
      allocate(ZmatVarNames(3*na))

      ! Read units for lengths and angles
      if (Node.eq.0) then
        lengthunits = fdf_string('ZM.UnitsLength',unitslen_default)
        if (leqi(lengthunits,'Bohr')) then
          ZmatUnitsLen = 0
        elseif (leqi(lengthunits,'Ang').or.leqi(lengthunits,'Angstrom')) then
          ZmatUnitsLen = 1
        else
          call die('Invalid Zmatrix length units')
        endif

        angleunits = fdf_string('ZM.UnitsAngle',unitsang_default)
        if (leqi(angleunits,'rad').or.leqi(angleunits,'radians')) then
          ZmatUnitsAng = 0
        elseif (leqi(angleunits,'deg').or.leqi(angleunits,'degrees')) then
          ZmatUnitsAng = 1
        else
          call die('Invalid Zmatrix angular units')
        endif
      endif
#     ifdef MPI
      if (Node.eq.0) then
        integers(1) = ZmatUnitsLen
        integers(2) = ZmatUnitsAng
      endif
      call MPI_Bcast(integers(1),2,MPI_integer,0,MPI_Comm_World,MPIerror)
      ZmatUnitsLen = integers(1)
      ZmatUnitsAng = integers(2)
#     endif

      ! Read maximum force tolerance for lengths and angles        
      if (Node.eq.0) then
        ZmatForceTolLen = fdf_physical('ZM.ForceTolLength',ftollen_default,'Ry/Bohr')
        ZmatForceTolAng = fdf_physical('ZM.ForceTolAngle',ftolang_default,'Ry/rad')

        ! Read maximum displacement per CG step for lengths and angles
        ZmatMaxDisplLen = fdf_physical('ZM.MaxDisplLength',dxmaxlen_default,'Bohr')
        ZmatMaxDisplAng = fdf_physical('ZM.MaxDisplAngle', dxmaxang_default,'rad')

      endif

#     ifdef MPI
      if (Node.eq.0) then
        physicals(1) = ZmatForceTolLen
        physicals(2) = ZmatForceTolAng
        physicals(3) = ZmatMaxDisplLen
        physicals(4) = ZmatMaxDisplAng
      endif
      call MPI_Bcast(physicals(1),4,MPI_double_precision,0,MPI_Comm_World,MPIerror)
      ZmatForceTolLen = physicals(1)
      ZmatForceTolAng = physicals(2)
      ZmatMaxDisplLen = physicals(3)
      ZmatMaxDisplAng = physicals(4)
#     endif


      ! Check if we should calculate all the forces
      if (Node.eq.0) then
        lCalcAllForces = fdf_boolean('ZM.CalcAllForces',lCalcAllForces_default)
      endif

#     ifdef MPI
      call MPI_Bcast(lCalcAllForces,1,MPI_logical,0, MPI_Comm_World,MPIerror)
#     endif
    endif

      ! Check whether a Z-matrix has been input
    if (Node.eq.0) then
      lUseZmatrix = fdf_block('Zmatrix',iu)
    endif
#    ifdef MPI
    call MPI_Bcast(lUseZmatrix,1,MPI_logical,0, MPI_Comm_World,MPIerror)
#   endif

    ! If not Z matrix return
    if (.not.lUseZmatrix) goto 999

    ! Output information about Z matrix units only for Z matrix case
    if (Node.eq.0) then
      write (6,'(''read_Zmatrix: Length units: '',a10)') lengthunits
      write (6,'(''read_Zmatrix: Angle  units: '',a10)') angleunits
      write (6,'(''read_Zmatrix: Force tolerances:'')')
      write (6,'(''read_Zmatrix:    for lengths = '',f12.6,'' Ry/Bohr'')') ZmatForceTolLen
      write (6,'(''read_Zmatrix:    for angles  = '', f12.6,'' Ry/rad'',/)') ZmatForceTolAng
      write (6,'(''read_Zmatrix: Maximum displacements:'')')
      write (6,'(''read_Zmatrix:    for lengths = '', f12.6,'' Bohr'')') ZmatMaxDisplLen
      write (6,'(''read_Zmatrix:    for angles  = '', f12.6,'' rad'')') ZmatMaxDisplAng
    endif

    ! Initialise number of molecules and Cartesian blocks
    nZmol = 0
    nZcart = 0
    nVars = 0

    if (Node.eq.0) then
      ! Read data block up to endblock statement
      eof = .false.
      nAtoms = 0
      write(*,"(a)") "%block Zmatrix"

      read(iu,'(a)',end=10,err=10) line
      do while (.not.eof.and.index(line,'%endblock').eq.0)
        write(*,"(a)") trim(line)
        lastc = index(line,'#') - 1
        if (lastc .le. 0) lastc = len(line)
        call upper2lower(line,lastc)
        if (index(line(1:lastc),'#').ne.0) then
          continue
        elseif (index(line(1:lastc),'mol').ne.0) then
          ! Found new molecule          
          nZmol = nZmol + 1
          nZmolAtoms(nZmol) = 0
          nZmolStartAtom(nZmol) = nAtoms + 1
          input_type = 0
          if (index(line(1:lastc),'scale').ne.0) then
            units_type = 1
          elseif (index(line(1:lastc),'frac').ne.0) then
            units_type = 2
          else
            units_type = 0
          endif
        elseif (index(line(1:lastc),'cart').ne.0) then
          ! Found new cartesian block
          nZcart = nZcart + 1
          nZcartAtoms(nZcart) = 0
          nZcartStartAtom(nZcart) = nAtoms + 1
          input_type = 1
          units_type = 0
        elseif (index(line(1:lastc),'scale').ne.0) then
          ! Found new cartesian(scaled) block
          nZcart = nZcart + 1
          nZcartAtoms(nZcart) = 0
          nZcartStartAtom(nZcart) = nAtoms + 1
          input_type = 1
          units_type = 1
        elseif (index(line(1:lastc),'frac').ne.0) then
          ! Found new cartesian(fractional) block
          nZcart = nZcart + 1
          nZcartAtoms(nZcart) = 0
          nZcartStartAtom(nZcart) = nAtoms + 1
          input_type = 1
          units_type = 2
        elseif (index(line(1:lastc),'constant').ne.0) then
          input_type = 2
        elseif (index(line(1:lastc),'variable').ne.0) then
          input_type = 3
        elseif (index(line(1:lastc),'constraint').ne.0) then
          input_type = 4

        else 
          ! Parse line
          call parse( line(1:lastc), nn, lc, names, nv, values, ni, integs, nr, reals, order )

          ! Check that the lines contains the correct type of data
          if (input_type.eq.0) then
            ! Molecule            
            do i=1,4
              if (order(i).ne.1) then
                call die('read_Zmatrix: Error in Z-matrix syntax') 
              endif
            enddo
            do i=5,7
              if (order(i).ne.0.and.order(i).ne.2) then
                call die('read_Zmatrix: Error in Z-matrix syntax') 
              endif
            enddo    
          elseif (input_type.eq.1) then
            ! Cartesian            
            if (order(1).ne.1) then
              call die ('read_Zmatrix: Error in Z-matrix syntax')
            endif
            do i=2,4
              if (order(i).ne.0.and.order(i).ne.2) then
                call die ('read_Zmatrix: Error in Z-matrix syntax')
              endif
            enddo
          elseif (input_type.eq.2.or.input_type.eq.3) then
            ! Constant/variable definition
            if (order(1).ne.0.or.order(2).ne.2) then
              call die ('read_Zmatrix: Error in Z-matrix syntax')
            endif
          elseif (input_type.eq.4) then
            ! Constraint definition              
            if (order(1).ne.0.or.order(2).ne.0.or. order(3).ne.2.or.order(4).ne.2) then
              call die ('read_Zmatrix: Error in Z-matrix syntax')
            endif
          endif

          ! Increment molecule/cartesian block count
          if (input_type.eq.0) then
            nZmolAtoms(nZmol) = nZmolAtoms(nZmol) + 1
          else if (input_type.eq.1) then
            nZcartAtoms(nZcart) = nZcartAtoms(nZcart) + 1
          endif

          ! Molecule/cartesian block
          if (input_type.eq.0.or.input_type.eq.1) then

            ! Increment atom number 
            nAtoms = nAtoms + 1
            ! Check number of atoms against input number
            if (nAtoms.gt.na) then
             call die('read_Zmatrix: Too many atoms in Z-matrix')
            endif

            ! Assign values for species and Z matrix
            nSpecies(nAtoms) = integs(1)
          endif

          ! Read dependencies if its a molecule
          if (input_type.eq.0) then
            nStart = nZmolStartAtom(nZmol)
            iZmat(3*nAtoms-2) = integs(2) + nStart-1
            iZmat(3*nAtoms-1) = integs(3) + nStart-1
            iZmat(3*nAtoms)   = integs(4) + nStart-1
          elseif (input_type.eq.1) then
            iZmat(3*nAtoms-2) = 1
            iZmat(3*nAtoms-1) = 1
            iZmat(3*nAtoms) = 1
          endif

          ! Assign type to coordinates
          if (input_type.eq.0) then
            if (nZmolAtoms(nZmol).eq.1) then
              ! Molecule, first atom
              do k=1,3
                ZmatType(3*(nAtoms-1)+k) = 3+units_type
              enddo
            else
              ! Molecule, not first atom
              ZmatType(3*nAtoms-2) = 2
              ZmatType(3*nAtoms-1) = 1
              ZmatType(3*nAtoms)   = 1
            endif
          elseif (input_type.eq.1) then
            ! Cartesian block              
            do k=1,3 
              ZmatType(3*(nAtoms-1)+k) = 6+units_type
            enddo
          endif

          ! Read reals/symbols for values of coordinates
          if (input_type.eq.0.or.input_type.eq.1) then
            nreals = 0
            nnames = 0
            do k=1,3
              flagnr = 4-3*input_type+k
              indi = 3*(nAtoms-1)+k

              if (order(flagnr).eq.0) then
                ! String/symbol
                nnames = nnames + 1
                i = 1
                found = .false.
                do while (i.le.nVars.and..not.found)
                  found = (leqi(names(lc(nnames-1)+1:lc(nnames)), ZmatVarNames(i)))
                  i = i + 1
                enddo
                if (found) then
                  iZmattoVars(indi) = i-1
                  VaryZmat(indi) = .false.
                  j = iVarstoZmat(i-1) 
                  ! Check length/angle dependency incompatibility
                  if (ZmatType(indi).eq.1.and. ZmatType(j).ne.1) then
                      call die('read_Zmatrix: error - angle/length dependency')
                  endif
                  do while (iNextDept(j).ne.0) 
                    j = iNextDept(j)
                  enddo
                  iNextDept(j) = indi
                else
                  nVars = nVars + 1
                  lZmatVarsDef(nVars) = .false.
                  ZmatVarNames(nVars) = names(lc(nnames-1)+1:lc(nnames))
                  iZmattoVars(indi) = nVars
                  iVarstoZmat(nVars) = indi
                endif
                coeffA(indi) = 1.0d0
                coeffB(indi) = 0.0d0
              else
                ! Explicit values                
                nreals = nreals + 1
                Zmat(indi) = reals(nreals)
                iZmattoVars(indi) = 0
                coeffA(indi) = 1.0d0
                coeffB(indi) = 0.0d0
                ! Optional flags for control of optimisation
                if (ni.ge.flagnr) then
                  VaryZmat(indi) = (integs(flagnr).ne.0)
                else
                  VaryZmat(indi) = .true.
                endif
              endif
              iNextDept(indi) = 0
            enddo
          endif

          ! Read variable/constant definitions
          if (input_type.eq.2.or.input_type.eq.3) then
            i = 1
            found = .false.
            do while (i.le.nVars.and..not.found) 
              found =  (leqi(ZmatVarNames(i),names(1:lc(1))))
              i = i+1
            enddo
            if (found) then
              if (lZmatVarsDef(i-1)) then
                call die('read_Zmatrix: Multiple definition of Z-matrix symbol')
              else
                indi = iVarstoZmat(i-1)
                do while (indi.ne.0) 
                  Zmat(indi) = reals(1)
                  indi = iNextDept(indi)
                enddo
                lZmatVarsDef(i-1) = .true.
              endif
              if (input_type.eq.2) then
              ! Fixed, symbolic                  
                VaryZmat(iVarstoZmat(i-1)) = .false.
              else
                ! Vary, symbolic                  
                VaryZmat(iVarstoZmat(i-1)) = .true.
              endif
            else
              call die('read_Zmatrix: Invalid symbol in Z-matrix')
            endif
          endif

          ! Read constraint definitions              
          if (input_type.eq.4) then
            i = 1
            found = .false.
            do while (i.le.nVars.and..not.found)
              found = leqi(ZmatVarNames(i),names(1:lc(1)))
              i = i + 1
            enddo
            if (found) then
              if (lZmatVarsDef(i-1)) then
                call die('read_Zmatrix: Multiple definition of Z-matrix symbol')
              else
                indi = iVarstoZmat(i-1)
                lZmatVarsDef(i-1) = .true.
              endif
            else
              call die('read_Zmatrix: Invalid symbol in Z-matrix')
            endif

            i = 1
            found = .false.
            do while (i.le.nVars.and..not.found)
              found = leqi(ZmatVarNames(i),names(lc(1)+1:lc(2)))
              i = i + 1
            enddo
            if (found) then
              if (.not.lZmatVarsDef(i-1)) then
                call die('read_Zmatrix: Invalid dependency for Z-matrix symbol')
              else
                indi1 = iVarstoZmat(i-1)
              endif
            else
              call die('read_Zmatrix: Invalid symbol in Z-matrix')
            endif
            ! If coefficient A=0 then its really a constant definition
            if (reals(1).eq.0.0d0) then
              call die( 'read_Zmatrix: erroneous constraints definition - '//&
                        'rather put in constants definition ')
            endif

            ! Put values into the newly defined Zmat components
            j = indi
            do while (j.ne.0) 
              Zmat(j) = reals(1)*Zmat(indi1)+reals(2)
              j = iNextDept(j)
            enddo
            ! Decide which is the independent variable
            if (indi.gt.indi1) then
              A = reals(1)
              B = reals(2)
            else
              A = 1.0d0/reals(1)
              B = -reals(2)/reals(1)
              j = indi
              indi = indi1
              indi1 = j
              VaryZmat(indi1) = VaryZmat(indi)
            endif
            VaryZmat(indi) = .false.
            ! Now indi1 is the independent var with indi = A*indi1+B
            ! Link the dependency lists of the two variables and update
            !   coeffiecients and iZmattoVars array
            j = indi1
            do while (iNextDept(j).ne.0) 
              j = iNextDept(j)
            enddo
            iNextDept(j) = indi
            j = indi
            do while (j.ne.0)
              iZmattoVars(j) = iZmattoVars(indi1)
              coeffA(j) = A*coeffA(indi1)
              coeffB(j) = B + A*coeffB(indi1)
              j = iNextDept(j)
            enddo

          endif

        endif
        ! End of parsing input line            
        read(iu,'(a)',end=10,err=10) line

      enddo
      10      continue
      ! Done reading Z-matrix
      write(*,"(a)") "%endblock Zmatrix"

    endif
    ! Done with if(Node.eq.0)        

#   ifdef MPI
    ! Distribute information over processors
    if (Node.eq.0) then
      integers1(1) = nAtoms
      integers1(2) = nZmol
      integers1(3) = nZcart
      integers1(4) = nVars
    endif
    call MPI_Bcast( integers1(1),4,MPI_integer,0,MPI_Comm_World,MPIerror )
    nAtoms = integers1(1)
    nZmol  = integers1(2)
    nZcart = integers1(3)
    nVars  = integers1(4)

    call MPI_Bcast( nSpecies(1), nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( nZmolAtoms(1), nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( nZmolStartAtom(1), nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( nZcartAtoms(1), nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( nZcartStartAtom(1), nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( iZmat(1), 3*nAtoms, MPI_integer, 0, MPI_Comm_World, MPIerror )
    call MPI_Bcast( VaryZmat(1), 3*nAtoms, MPI_logical, 0, MPI_Comm_World, MPIerror )
    call MPI_Bcast( Zmat(1), 3*nAtoms, MPI_double_precision, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( coeffA(1), 3*nAtoms, MPI_double_precision, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( coeffB(1), 3*nAtoms, MPI_double_precision, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( ZmatType(1), 3*nAtoms, MPI_integer, 0, MPI_Comm_World, MPIerror )
    call MPI_Bcast( iZmattoVars(1), 3*nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( iVarstoZmat(1), 3*nAtoms, MPI_integer, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( lZmatVarsDef(1), 3*nAtoms, MPI_logical, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( ZmatVarNames(1), 3*nAtoms*10, MPI_character, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( iNextDept(1), 3*nAtoms, MPI_integer, 0, MPI_Comm_World, &
                    MPIerror )

#   endif

    ! Check that all variables have been defined
    do i = 1,nVars
      if (.not.lZmatVarsDef(i)) then
        call die('read_Zmatrix: Symbol has not been given a value')
      endif
    enddo

    ! Scale all coordinates according to their type            
    do i = 1,3*nAtoms
      ! Angle
      if (ZmatType(i).eq.1) then
        Zmat(i) = Zmat(i) * (1.0d0+(deg-1.0d0)*ZmatUnitsAng)
        coeffB(i) = coeffB(i) * (1.0d0+(deg-1.0d0)*ZmatUnitsAng)
      endif

      ! Bond length/pure cartesian
      if (ZmatType(i).eq.2.or.ZmatType(i).eq.3.or. ZmatType(i).eq.6) then
        Zmat(i) = Zmat(i) * (1.0d0+(Ang-1.0d0)*ZmatUnitsLen)
        coeffB(i) = coeffB(i) * (1.0d0+(Ang-1.0d0)*ZmatUnitsLen)
      endif
         
      ! Scaled cartesian
      if (ZmatType(i).eq.4.or.ZmatType(i).eq.7) then
        Zmat(i) = Zmat(i) * alat
        coeffB(i) = coeffB(i) * alat
      endif
                
    enddo

    ! Translate pure cartesians by the origin
    do i=1,nAtoms
      if ((ZmatType(3*i).eq.3.or.ZmatType(3*i).eq.6).and.lOrigin) then
        Zmat(3*i-2) = Zmat(3*i-2) + origin(1)
        Zmat(3*i-1) = Zmat(3*i-1) + origin(2)
        Zmat(3*i) = Zmat(3*i) + origin(3)
      endif
    enddo

    ! Fractional cartesian 
    do i=1,nAtoms
      if (ZmatType(3*i).eq.5.or.ZmatType(3*i).eq.8) then
        do k=1,3
          ZmatVal(k) = Zmat(3*(i-1)+k)
          coeffBVal(k) = coeffB(3*(i-1)+k)
        enddo
        do k=1,3
          Zmat(3*(i-1)+k) = ucell(k,1)*ZmatVal(1) + ucell(k,2)*ZmatVal(2) +&
                            ucell(k,3)*ZmatVal(3)
          coeffB(3*(i-1)+k) = ucell(k,1)*coeffBVal(1) + ucell(k,2)*coeffBVal(2) +&
                              ucell(k,3)*coeffBVal(3)
        enddo
      endif
    enddo

    ! Set remaining undefined dependencies to 1 (they are not used but if they 
    !   are 0 the program might crash later)
    do m=1,nZmol
      nStart = nZmolStartAtom(m)
      iZmat(3*nStart-2) = 1
      iZmat(3*nStart-1) = 1
      iZmat(3*nStart) = 1
      iZmat(3*nStart+2) = 1
      iZmat(3*nStart+3) = 1
      iZmat(3*nStart+6) = 1
    enddo

    ! Check Z matrix data for consistency
    !   - are all atoms defined in terms of previous atoms?
    !   - are the three dependencies for every atom distinct?

    ! Loop over molecules        
    do m = 1,nZmol
      nStart = nZmolStartAtom(m)
      nAtoms = nZmolAtoms(m)
      ! Second atom
      if (nAtoms.gt.1) then
        if (iZmat(3*nStart+1).ne.nStart) then
          if (Node.eq.0) then
            write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
            write(6,'(''read_Zmatrix: molecule nr '',i7, ''; atom nr 2'')') m
          endif
          call die()
        endif
      endif

      ! Third atom
      if (nAtoms.gt.2) then
        if ( .not. (iZmat(3*nStart+4) .eq. nStart .and. iZmat(3*nStart+5) .eq. &
              nStart+1 .or. iZmat(3*nStart+5) .eq. nStart .and. iZmat(3*nStart+4) &
              .eq.nStart+1)) then
          if (Node.eq.0) then
            write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
            write(6,'(''read_Zmatrix: molecule nr '',i7, ''; atom nr 3'')') m
          endif
          call die()
        endif
      endif

      ! Fourth atom and up -> general case          
      do i = nStart+3,nStart+nAtoms-1
        if (iZmat(3*i-2).gt.i-1.or.iZmat(3*i-2).eq.iZmat(3*i-1)) then
          if (Node.eq.0) then
            write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
            write(6,'(''read_Zmatrix: molecule nr '',i7, ''; atom nr '',i7)') &
                    m,i-nStart
          endif
          call die() 
        endif
        if (iZmat(3*i-1).gt.i-1.or. iZmat(3*i-1).eq.iZmat(3*i)) then
          if (Node.eq.0) then
            write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
            write(6,'(''read_Zmatrix: molecule nr '',i7,''; atom nr '',i7)') &
                      m,i-nStart
          endif
          call die()
        endif
        if (iZmat(3*i).gt.i-1.or.iZmat(3*i).eq.iZmat(3*i-2)) then
          if (Node.eq.0) then
            write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
            write(6,'(''read_Zmatrix: molecule nr '',i7,''; atom nr '',i7)') &
                      m,i-nStart
          endif
          call die() 
        endif
      enddo
    ! End looping over molecules          
    enddo

    if (ZmatUnitsAng == 0) then
       scale_angle = 1.0_dp
    else
       scale_angle = deg
    endif

    if (ZmatUnitsLen == 0) then
       scale_length = 1.0_dp
    else
       scale_length = Ang
    endif

    ! Output Z-matrix coordinates
    if (Node.eq.0) then
      call write_Zmatrix()
    endif

    ! Return point
    999   continue
    firsttime = .false.
    !----------------------------------------------------------------------- END
!1   format((''read_Zmatrix: Length units: '',a10))
!2   format((''read_Zmatrix: Angle  units: '',a10))
!3   format((''read_Zmatrix: Force tolerances:''))
!4   format((''read_Zmatrix:    for lengths = '',f12.6,'' Ry/Bohr''))
!5   format((''read_Zmatrix:    for angles  = '', f12.6,'' Ry/rad'',/))
!6   format((''read_Zmatrix: Maximum displacements:''))

  end subroutine read_Zmatrix

  subroutine Zmat_to_Cartesian(Cartesian)
    !----------------------------------------------------------- Input Variables
    real(dp), intent(inout) ::  Cartesian(:,:)
    !----------------------------------------------------------- Local Variables
    integer                 :: i
    integer                 :: m
    integer                 :: nStart
    integer                 :: nAtoms
    real(dp)                :: RelatedC(3,3)
    real(dp)                :: x
    real(dp)                :: y
    real(dp)                :: z
    !--------------------------------------------------------------------- BEGIN

    !  Loop over molecules
    do m = 1,nZmol

      nStart = nZmolStartAtom(m)
      nAtoms = nZmolAtoms(m)

      !  Loop over atoms within each molecule generating Cartesian coordinates
      !  Use general Z2CGen subroutine which takes special first 3 atoms into account
      do i = nStart,nStart+nAtoms-1
        RelatedC(1:3,1) = Cartesian(1:3,iZmat(3*i-2))
        RelatedC(1:3,2) = Cartesian(1:3,iZmat(3*i-1))
        RelatedC(1:3,3) = Cartesian(1:3,iZmat(3*i))
        call Z2CGen(i-nStart,Zmat(3*i-2),Zmat(3*i-1),Zmat(3*i), RelatedC,x,y,z)
        Cartesian(1,i) = x
        Cartesian(2,i) = y
        Cartesian(3,i) = z
      enddo

    enddo

    !  Loop over cartesian blocks
    do m = 1,nZcart

      nStart = nZcartStartAtom(m)
      nAtoms = nZcartAtoms(m)

      do i = nStart,nStart+nAtoms-1
        Cartesian(1,i) = Zmat(3*i-2)
        Cartesian(2,i) = Zmat(3*i-1)
        Cartesian(3,i) = Zmat(3*i)
      enddo

    enddo
    !----------------------------------------------------------------------- END

  end subroutine Zmat_to_Cartesian

  subroutine iofaZmat()
  
    implicit none
    !----------------------------------------------------------- Local Variables
    integer          :: m
    integer          :: i
    integer          :: k
    integer          :: nStart
    integer          :: nAtoms
    integer          :: jindex
    character(len=4) :: lenstr
    character(len=4) :: angstr
    !--------------------------------------------------------------------- BEGIN

    if (ZmatUnitsLen.eq.0) then
      lenstr = 'Bohr'
    elseif (ZmatUnitsLen.eq.1) then
      lenstr = 'Ang'
    endif
    if (ZmatUnitsAng.eq.0) then
      angstr = 'rad'
    elseif (ZmatUnitsAng.eq.1) then
      angstr = 'deg'
    endif

    write(6,'(/,''zmatrix: Atomic forces (eV/'',a4, ''; eV/'',a4,'')'')') &
              lenstr,angstr
    write(6,'(''zmatrix: (No information if symbols are used)'')')
    !  Print molecule forces in the user's units
    do m = 1,nZmol
      nStart = nZmolStartAtom(m)
      nAtoms = nZmolAtoms(m)
      write(6,'(''molecule'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
      write(6,1) nStart,(ZmatForce(3*(nStart-1)+k)*(1+(Ang-1)*ZmatUnitsLen)/eV,k=1,3)
      write(6,1) (i,ZmatForce(3*i-2)*(1+(Ang-1)*ZmatUnitsLen)/eV,  &
                    ZmatForce(3*i-1)*(1+(deg-1)*ZmatUnitsAng)/eV,  &
                    ZmatForce(3*i)*(1+(deg-1)*ZmatUnitsAng)/eV,    &
                    i=nStart+1,nStart+nAtoms-1)
    enddo

    !  Print cartesian forces in the user's units
    do m = 1,nZcart
      nStart = nZcartStartAtom(m)
      nAtoms = nZcartAtoms(m)
      write(6,'(''cartesian'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
      write(6,1) (i,(ZmatForce(3*(i-1)+k)*(1+(Ang-1)*ZmatUnitsLen)/eV,k=1,3),i=nStart,nStart+nAtoms-1)
    enddo
    write (6,'(/)') 

    !  Print variables' forces in user's units
    write(6,'(''Variable forces (eV/'',a4, ''; eV/'',a4,'')'')') lenstr,angstr
    do i = 1,nVars
       jindex = iVarstoZmat(i)
       if (.not.VaryZmat(jindex)) cycle

      if (ZmatType(iVarstoZmat(i)).eq.1) then
        write(6,*) ZmatVarNames(i), ZmatForceVar(i)* (1+(deg-1)*ZmatUnitsAng)/eV
      else 
        write(6,*) ZmatVarNames(i), ZmatForceVar(i)*(1+(Ang-1)*ZmatUnitsLen)/eV
      endif
    enddo

    return
    !----------------------------------------------------------------------- END

1   format((i6,3f12.6))

  end subroutine iofaZmat

  subroutine CartesianForce_to_ZmatForce(na,Cartesian,CartesianForce)

    use parallel,    only : IONode
#   ifdef MPI
    use mpi_siesta
#   endif

    implicit none

    !----------------------------------------------------------- Input Variables
    integer,  intent(in)  :: na
    real(dp), intent(in)  :: Cartesian(3,na)
    real(dp), intent(in)  :: CartesianForce(3,na)
    !----------------------------------------------------------- Local Variables
    integer                     :: i
    integer                     :: j
    integer                     :: k
    integer                     :: ia
    integer                     :: m
    integer                     :: vi
    real(dp), parameter         :: hfrac = 1e-6
    real(dp), allocatable       :: CartesianB(:,:)
    real(dp), allocatable       :: CartesianF(:,:)
    real(dp), allocatable       :: ZmatF(:)
    real(dp), allocatable       :: ZmatB(:)
    real(dp)                    :: h
    logical,               save :: firsttime = .true.
#   ifdef MPI
    integer                     :: MPIerror
#   endif        
    !--------------------------------------------------------------------- BEGIN

    if (IOnode) then
      !  Allocate local workspace arrays
      allocate(CartesianF(3,na))
      allocate(CartesianB(3,na))
      allocate(ZmatF(3*na))
      allocate(ZmatB(3*na))

      !  Initialise forces & cartesian transformation matrices
      ZmatForce(1:3*na) = 0.0_dp
      ZmatForceVar(1:nVars) = 0.0_dp
      CartesianF(1:3,1:na) = 0.0_dp
      CartesianB(1:3,1:na) = 0.0_dp

      !  Reset the FD and BD Z-Zmatrices
      ZmatF(1:3*na) = Zmat(1:3*na)
      ZmatB(1:3*na) = Zmat(1:3*na)

      ! Zero the Forward & Backward Cartesian matrices:
      CartesianF = 0.0_dp
      CartesianB = 0.0_dp

      !  Loop over all atoms
      do ia = 1,na
        !  Loop over all coordinates of the atom
        do k = 1,3
          !  Index for this coordinate          
          i = 3*(ia-1) + k
          !  Only proceed if this coordinate is allowed to vary or CalcAllForces is true            
          if (VaryZmat(i).or.lCalcAllForces) then
            !  Set the h-increment for the coordinate
            if (Zmat(i).eq.0.0_dp) then 
              h = hfrac
            else 
              h = hfrac*Zmat(i)
            endif
            !  Update the FD and BD Zmatrices for this coordinate
            ZmatF(i) = Zmat(i) + h
            ZmatB(i) = Zmat(i) - h
            !  Calculate the force on this coordinate
            if (iZmattoVars(i).eq.0.or.lCalcAllForces) then
              ZmatForce(i) = GetForce(na,Cartesian,CartesianForce,CartesianF,CartesianB, ZmatF,ZmatB, i,h,.false.)
            endif
            !  Check if variable is independent              
            vi = iZmattoVars(i)
            j = 0
            if (vi.ne.0) then
              j = iVarstoZmat(vi)
            endif
            !  Calculate the force on this variable
            if (j.eq.i) then
              !  Update the FD and BD Z-matrices for all dependents
              do while (iNextDept(j).ne.0) 
                j = iNextDept(j)
                ZmatF(j) = Zmat(j) + h*coeffA(j)
                ZmatB(j) = Zmat(j) - h*coeffA(j)
              enddo
              ZmatForceVar(vi) = GetForce(na,Cartesian,CartesianForce, CartesianF,CartesianB, ZmatF,ZmatB,i,h,.true.)
              !  reset the FD and BD Z-matrices
              j = i
              do while (iNextDept(j).ne.0)
                j = iNextDept(j)
                ZmatF(j) = Zmat(j)
                ZmatB(j) = Zmat(j)
              enddo
            endif
            ZmatF(i) = Zmat(i)
            ZmatB(i) = Zmat(i)
          endif
        enddo
        !  Done looping over 3 coordinate of this atom          
        !  From now on the previous atoms will always be unchanged
        CartesianF(1:3,ia) = Cartesian(1:3,ia)
        CartesianB(1:3,ia) = Cartesian(1:3,ia)
      enddo
      !  Done looping over all atoms

      !  Deallocate workspace arrays
      call memory('D','D',size(ZmatB),'read_Zmatrix')
      deallocate(ZmatB)
      call memory('D','D',size(ZmatF),'read_Zmatrix')
      deallocate(ZmatF)
      call memory('D','D',size(CartesianB),'read_Zmatrix')
      deallocate(CartesianB)
      call memory('D','D',size(CartesianF),'read_Zmatrix')
      deallocate(CartesianF)

    endif !IOnode

#   ifdef MPI
    call MPI_Bcast( ZmatForce, 3*na,MPI_double_precision, 0, MPI_Comm_World,&
                    MPIerror )
    call MPI_Bcast( ZmatForceVar, nVars, MPI_double_precision, 0, MPI_Comm_World,&
                    MPIerror )
#   endif        
    !----------------------------------------------------------------------- END

  end subroutine CartesianForce_to_ZmatForce

  function GetForce(na,Cartesian,CartesianForce, CartesianF,CartesianB, ZmatF,ZmatB,i,h,doAll)

    implicit none

    !----------------------------------------------------------- Input Variables
    integer, intent(in)           :: i
    logical, intent(in)           :: doAll
    integer,  intent(in)          :: na
    real(dp), intent(in)          :: Cartesian(3,na)
    real(dp), intent(in)          :: CartesianForce(3,na)
    real(dp), intent(inout)       :: CartesianF(3,na)
    real(dp), intent(inout)       :: CartesianB(3,na)
    real(dp), intent(inout)       :: ZmatF(3*na)
    real(dp), intent(inout)       :: ZmatB(3*na)
    real(dp), intent(in)          :: h

    !----------------------------------------------------------- Local Variables
    real(dp)        :: GetForce
    integer         :: m
    integer         :: prevm
    integer         :: j
    integer         :: l
    integer         :: ca
    real(dp)        :: RelatedC(3,3)
    real(dp)        :: X,Y,Z
    !--------------------------------------------------------------------- BEGIN

    !  Initialize force
    GetForce = 0.0d0
    !  Add force contributions from all the affected cartesian coordinates 
    m = 1
    prevm = 0
    !  Loop over affected Zmat Coords            
    j = i
    do while (j.ne.0) 
    !  l is the atom number for index j            
      l = (j+2)/3
      if (ZmatType(j).lt.6) then
        !  Molecule              
        !  Find correct molecule
        do while (l.gt.nZmolStartAtom(m)+nZmolAtoms(m)-1)
          m = m + 1
        enddo
       !  Only add force contributions if this is the first affected 
        !     coordinate for a specific molecule (otherwise it has been done already)
        if (m.ne.prevm) then
          prevm = m
          !  Initialise FD and BD cartesians that come before changed atom
          do ca = nZmolStartAtom(m),l-1
            CartesianF(1:3,ca) = Cartesian(1:3,ca)
            CartesianB(1:3,ca) = Cartesian(1:3,ca)
          enddo
          !  Add force contributions from all atoms later than and
          !     including the changed one in this molecule
          do ca = l,nZmolStartAtom(m)+nZmolAtoms(m)-1
            !  Forward cartesian                
            RelatedC(1:3,1) = CartesianF(1:3,iZmat(3*ca-2))
            RelatedC(1:3,2) = CartesianF(1:3,iZmat(3*ca-1))
            RelatedC(1:3,3) = CartesianF(1:3,iZmat(3*ca))
            call Z2CGen( ca-nZmolStartAtom(m), ZmatF(3*ca-2), ZmatF(3*ca-1),       &
                         ZmatF(3*ca), RelatedC, CartesianF(1,ca), CartesianF(2,ca),&
                         CartesianF(3,ca) )
            !  Backward cartesian               
            RelatedC(1:3,1) = CartesianB(1:3,iZmat(3*ca-2))
            RelatedC(1:3,2) = CartesianB(1:3,iZmat(3*ca-1))
            RelatedC(1:3,3) = CartesianB(1:3,iZmat(3*ca))
            call Z2CGen( ca-nZmolStartAtom(m), ZmatB(3*ca-2), ZmatB(3*ca-1),       &
                         ZmatB(3*ca), RelatedC, CartesianB(1,ca), CartesianB(2,ca),&
                         CartesianB(3,ca) )
            !  Add in force constribution     
            X = CartesianForce(1,ca)*(CartesianF(1,ca)-CartesianB(1,ca))/(2*h)
            Y = CartesianForce(2,ca)*(CartesianF(2,ca)-CartesianB(2,ca))/(2*h)
            Z = CartesianForce(3,ca)*(CartesianF(3,ca)-CartesianB(3,ca))/(2*h)
            GetForce = GetForce + X+Y+Z
            ! testing
          enddo
        endif
      else  
        !  Cartesian              
        GetForce = GetForce+CartesianForce(j-3*(l-1),l)*coeffA(j)
      endif  
      if (doAll) then
        j = iNextDept(j)
      else 
        j=0
      endif
    enddo
    !----------------------------------------------------------------------- END

  end function GetForce


  subroutine write_Zmatrix()

    implicit none
    !----------------------------------------------------------- Local Variables
    integer             :: nStart
    integer             :: nAtoms
    integer             :: k
    integer             :: i
    integer             :: m
    character(len=4)    :: lenstr
    character(len=4)    :: angstr
    !--------------------------------------------------------------------- BEGIN

    !  Output Z matrix information
    if (ZmatUnitsLen.eq.0) then
      lenstr = 'Bohr'
    elseif (ZmatUnitsLen.eq.1) then
      lenstr = 'Ang'
    endif
    if (ZmatUnitsAng.eq.0) then
      angstr = 'rad'
    elseif (ZmatUnitsAng.eq.1) then
      angstr = 'deg'
    endif


    !  Write molecule coordinates in user's units
    write(6,'(/,''zmatrix: Z-matrix coordinates: ('',a4, ''; '',a4,'')'')') lenstr,angstr
    write(6,'(''zmatrix: '', ''(Fractional coordinates have been converted '', ''to cartesian)'')')
    do m = 1,nZmol
      nStart = nZmolStartAtom(m)
      nAtoms = nZmolAtoms(m)
      write(6,'(''molecule'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
      write(6,'(3f16.8)') (Zmat(3*(nStart-1)+k)/(1.0d0+(Ang-1.0d0)* ZmatUnitsLen),k=1,3)
      do i=nStart+1,nStart+nAtoms-1
        write(6,'(3f16.8)') Zmat(3*i-2)/(1.0d0+(Ang-1.0d0)*ZmatUnitsLen), &
                            Zmat(3*i-1)/(1.0d0+(deg-1.0d0)*ZmatUnitsAng), &
                            Zmat(3*i)/(1.0d0+(deg-1.0d0)*ZmatUnitsAng)
      enddo
    enddo

    !  Write cartesian coordinates in user's units        
    do m = 1,nZcart
      nStart = nZcartStartAtom(m)
      nAtoms = nZcartAtoms(m)
      write(6,'(''cartesian block'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
      write(6,'(3f16.8)')((Zmat(3*(i-1)+k)/(1+(Ang-1)*ZmatUnitsLen), k=1,3),i=nStart,nStart+nAtoms-1)
    enddo
    write(6,'(/)')

    call print_variables()
    !----------------------------------------------------------------------- END

  end subroutine write_Zmatrix

  subroutine Z2CGen(atomNr,r,theta,phi,RelatedC,x,y,z)
    !  Subroutine for generation of cartesian coordinates from Z-matrix
    !  given the following special cases:
    !    - the first entry is pure cartesian coordinates
    !    - the second entry is spherical coordinates relative to the first atom
    !    - the third entry is general Z-matrix coordinates, but with the torsion atom 
    !      a dummy atom 1 unit in the z direction above the second
    !    - the fourth entry onwards are general Z-matrix coordinates
    implicit none
    !----------------------------------------------------------- Input Variables
    integer,  intent(in)     :: atomNr
    real(dp), intent(in)     :: r
    real(dp), intent(in)     :: theta
    real(dp), intent(in)     :: phi
    real(dp), intent(inout)  :: RelatedC(3,3)
    real(dp), intent(out)    :: x
    real(dp), intent(out)    :: y
    real(dp), intent(out)    :: z
    !--------------------------------------------------------------------- BEGIN

    if (atomNr.eq.0) then
      x = r
      y = theta
      z = phi
    elseif (atomNr.eq.1) then
      x = r*sin(theta)*cos(phi) + RelatedC(1,1)
      y = r*sin(theta)*sin(phi) + RelatedC(2,1)
      z = r*cos(theta) + RelatedC(3,1)
    elseif (atomNr.eq.2) then
      RelatedC(1,3) = RelatedC(1,1)
      RelatedC(2,3) = RelatedC(2,1)
      RelatedC(3,3) = RelatedC(3,1) + 1.0d0
      call Z2C(r,theta,phi,RelatedC,x,y,z)
    elseif(atomNr.gt.2) then
      call Z2C(r,theta,phi,RelatedC,x,y,z)
    endif
    !----------------------------------------------------------------------- END

  end subroutine Z2CGen

  subroutine Z2C(r,theta,phi,RelatedC,x,y,z)
    !  Subroutine for generation of Cartesian coordinates from Z-Matrix
    !  Julian Gale, NRI, Curtin University, March 2004
    implicit none
    !----------------------------------------------------------- Input Variables
    real(dp), intent(in)  :: r
    real(dp), intent(in)  :: theta
    real(dp), intent(in)  :: phi
    real(dp), intent(in)  :: RelatedC(3,3)
    real(dp), intent(out) :: x
    real(dp), intent(out) :: y
    real(dp), intent(out) :: z
    !----------------------------------------------------------- Local Variables
    integer           :: i
    integer           :: j
    integer           :: k
    real(dp)          :: rji
    real(dp)          :: rn
    real(dp)          :: rp
    real(dp)          :: xi
    real(dp)          :: yi
    real(dp)          :: zi
    real(dp)          :: xj
    real(dp)          :: yj
    real(dp)          :: zj
    real(dp)          :: xji
    real(dp)          :: yji
    real(dp)          :: zji
    real(dp)          :: xk
    real(dp)          :: yk
    real(dp)          :: zk
    real(dp)          :: xki
    real(dp)          :: yki
    real(dp)          :: zki
    real(dp)          :: xn
    real(dp)          :: yn
    real(dp)          :: zn
    real(dp)          :: xp
    real(dp)          :: yp
    real(dp)          :: zp
    !--------------------------------------------------------------------- BEGIN

    !  Find coordinates for related atoms
    xi = RelatedC(1,1)
    yi = RelatedC(2,1)
    zi = RelatedC(3,1)
    xj = RelatedC(1,2)
    yj = RelatedC(2,2)
    zj = RelatedC(3,2)
    xk = RelatedC(1,3)
    yk = RelatedC(2,3)
    zk = RelatedC(3,3)
    !  Find unit vector along j->i vector
    xji = xi - xj
    yji = yi - yj
    zji = zi - zj
    rji = xji*xji + yji*yji + zji*zji
    rji = sqrt(rji)
    xji = xji/rji
    yji = yji/rji
    zji = zji/rji
    !  Find j->k vector
    xki = xk - xj
    yki = yk - yj
    zki = zk - zj
    !  Find unit vector normal to the i-j-k plane
    xn = yji*zki - yki*zji
    yn = zji*xki - zki*xji
    zn = xji*yki - xki*yji
    rn = xn*xn + yn*yn + zn*zn
    rn = sqrt(rn)
    xn = xn/rn
    yn = yn/rn
    zn = zn/rn
    !  Find unit vector normal to the other 2 directions already found
    !  Since original vectors are normalised the result should be likewise
    xp = yn*zji - yji*zn
    yp = zn*xji - zji*xn
    zp = xn*yji - xji*yn
    !  Find distances along each unit vector
    rji = r*cos(theta)
    rn  = r*sin(theta)*sin(phi)
    rp  = r*sin(theta)*cos(phi)
    !  Multiply unit vector by distances and add to origin to get position
    x = xi - rji*xji + rn*xn + rp*xp
    y = yi - rji*yji + rn*yn + rp*yp
    z = zi - rji*zji + rn*zn + rp*zp
    !----------------------------------------------------------------------- END

  end subroutine z2c

  function user_print(i,offset) result(x)
    !     Converts variables and constants to the
    !     format originally supplied by the user.
    use m_cell, only  : cart2frac

    ! Index of Zmatrix variable
    integer, intent(in)  :: i
    real(dp), intent(in), optional :: offset
    real(dp)             :: x  ! final value

    real(dp) :: xin
    real(dp) :: r(3), rfrac(3)
    integer  :: natom, k
    !     Offset: prepared to deal with origin shift... *****
    !
    !  ZmatType(3*na) defines the component type:
    !    1 - angle
    !    2 - bond length
    !    3 - pure cartesian (molecule)
    !    4 - scaled cartesian (molecule)
    !    5 - fractional cartesian (molecule)
    !    6 - pure cartesian
    !    7 - scaled cartesian (scaled by lattice constant)
    !    8 - fractional cartesian (scaled by lattice vectors)
    !--------------------------------------------------------------------- BEGIN

    xin = Zmat(i)
    if (present(offset)) then
       xin = xin + offset
    endif

    select case (ZmatType(i))

    case(1)
      x =  xin/scale_angle
    case(2,3,6)
      x=   xin/scale_length
    case(4,7)
      x =  xin / zmatrix_alat
    case(5,8)

      ! Get other cart coords and compute fractional coords
      ! we assume that any use of fractional coordinates is
      ! on a whole-atom basis, i.e., there are no bond lenghts
      ! or angles involved, so Zmat(base+1:base+3) is homogeneous

      natom = 1 + (i-1)/3
      r(1:3) = Zmat(3*(natom-1)+1:3*(natom-1)+3)
      call cart2frac(r,rfrac)
      k = i - 3*(natom-1)
      x = rfrac(k)

    case default
      call die("Wrong type for Zmatrix coordinate")
    end select
    !----------------------------------------------------------------------- END

  end function user_print

  subroutine print_variables()
    use siesta_geom, only: ucell
    use m_cell, only  : celli
    implicit none

    integer  :: i
    integer  :: jindex
    integer  :: ivars

    external :: io_assign, io_close
    !--------------------------------------------------------------------- BEGIN

    !Update celli, just in case
    call reclat( ucell,celli,0 )
    call io_assign( ivars )
    open( ivars, file="ZMATRIX_VARS", form='formatted', position='rewind', &
                status='unknown')

    write(6,"(a)") "Z-matrix Symbol Section -------"
    write(6,"(a)") "Variables"
    do i = 1, nVars
      jindex = iVarstoZmat(i)
      if (.not.VaryZmat(jindex)) cycle
      print *, ZmatVarNames(i), user_print(jindex)
      write(unit=ivars, fmt="(a20,g25.15)") ZmatVarNames(i), user_print(jindex)
    enddo
    call io_close(ivars)

    write(6,"(a)") "Constants"
    do i = 1, nVars
       jindex = iVarstoZmat(i)
       if (VaryZmat(jindex)) cycle
       print *, ZmatVarNames(i), user_print(jindex)
    enddo
    write(6,"(a,/)") "------------ End of Z-matrix Information"
    !----------------------------------------------------------------------- END

  end subroutine print_variables

  subroutine upper2lower(string,nchar)
    !  upper2lower accepts a string of nchar characters and replaces
    !  any lowercase letters by uppercase ones.
    character :: string*(*)
    character :: char*1
    integer   :: nchar
    integer   :: ncopy
    integer   :: i
    integer   :: itemp
    !--------------------------------------------------------------------- BEGIN

    ncopy = nchar
    if (ncopy.le.0) then
      ncopy = len(string)
    endif
    do i = 1,ncopy
      if (lge(string(I:I),'A').and.lle(string(I:I),'Z')) then
        itemp = ichar(string(i:i))-ichar('A')+ichar('a')
        string(I:I) = char(itemp)
      endif
    enddo
    !----------------------------------------------------------------------- END

  end subroutine upper2lower

END MODULE zmatrix
