      module atminit

      use precision
      use atmparams

      use atmfuncs, only:species

      implicit none

      CONTAINS

      subroutine initatom(ns, na, isa, no, nokb)

C Routine to initialize the Pseudopotentials and Atomic Orbitals.
C
C Writen by J.Soler and P.Ordejon, August-October'96 
C Strongly modified by D. Sanchez-Portal, Oct. 1998 
C Modified by DSP, July 1999.
C Modified by JDG, October 1999 for dynamic memory usage
C **************************** INPUT **********************************
C integer ns                : Number of species
C integer na                : Number of atoms
C integer isa(na)           : Species index of each atom
C **************************** OUTPUT *********************************
C integer no                : Total number of orbitals
C integer nokb              : Total number of KB projectors
C *********************************************************************
C
#ifdef MPI
      use mpi
#endif

      use fdf

      integer, intent(in)  :: ns, na, isa(na)
      integer, intent(out) :: no, nokb

C Internal variables ...................................................
      integer ia, is, js, nkba, noa

      logical new_user_basis, atmonly

      integer Node

#ifdef MPI
      integer MPIerror
#endif
!
!     Automatic arrays
!
      logical           :: semic(ns)
      integer           :: lmxkbs(ns), lmaxs(ns)
      integer           :: nsemic(0:lmaxd,ns), nkbl(0:lmaxd,ns)
      integer           :: polorb(0:lmaxd,nsemx,ns)
      integer           :: nzls(0:lmaxd,nsemx,ns)
      real*8            :: erefkb(nkbmx,0:lmaxd,ns)
      real*8            :: charge(ns)
      real*8            :: contrf(nzetmx,0:lmaxd,nsemx,ns)
      real*8            :: rcls(nzetmx,0:lmaxd,nsemx,ns)
      integer           :: izs(ns), nkbs(ns), nos(ns)
      real*8            :: smass(ns)
      real*8            ::  qos(maxos,ns)

      character(len=10) :: basistype(ns)
      character(len=20) :: atomlabel(ns)
      
      external transfer, read_user_basis, broadcast_basis

C Get Node number  --- Put it in atmfuncs
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

c Reading input for the pseudopotentials and atomic orbitals 

      if (Node.eq.0) then
        write(6,'(/2a)') 
     .    'initatom: Reading input for the pseudopotentials ',
     .    'and atomic orbitals'
      endif

      if (Node.eq.0) 
     $    new_user_basis = fdf_boolean('new-user-basis',.false.)
#ifdef MPI
      call MPI_Bcast(new_user_basis,1,MPI_logical,
     $               0,MPI_Comm_World,MPIerror)
#endif
      if (new_user_basis) then
         if (Node.eq.0) call read_user_basis
         GOTO 1000    !! Do not use the old atomic routines
      endif

!
!     Old atomic routines
!
      if (Node.eq.0) then
        call read_basis(ns,izs,lmxkbs,nkbl,erefkb,lmaxs,
     .    nzls,rcls,contrf,atomlabel,polorb,semic,nsemic,
     .    charge,smass,basistype)
      endif

C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(izs,ns,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(lmxkbs,ns,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(nkbl(0,1),(lmaxd+1)*ns,MPI_integer,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(erefkb(1,0,1),nkbmx*(lmaxd+1)*ns,DAT_double,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(lmaxs,ns,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(nzls(0,1,1),(lmaxd+1)*nsemx*ns,MPI_integer,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(rcls(1,0,1,1),nzetmx*(lmaxd+1)*nsemx*ns,
     .  DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(contrf(1,0,1,1),nzetmx*(lmaxd+1)*nsemx*ns,
     .  DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(atomlabel,20*ns,MPI_character,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(polorb(0,1,1),(lmaxd+1)*nsemx*ns,MPI_integer,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(semic,ns,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(nsemic(0,1),(lmaxd+1)*ns,MPI_integer,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(charge,ns,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(smass,ns,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(basistype,10*ns,MPI_character,0,
     .  MPI_Comm_World,MPIerror)
#endif

C .....................
c Initialize pseudopotentials and atomic orbitals 
      do is = 1,ns
        call atom( ns, izs(is), lmxkbs(is), nkbl(0,is), 
     .           erefkb(1,0,is),lmaxs(is), nzls(0,1,is), rcls(1,0,1,is), 
     .           contrf(1,0,1,is),
     .           atomlabel(is), polorb(0,1,is), semic(is), nsemic(0,is),
     .           charge(is), smass(is), basistype(is), 
     .           js, nkbs(is), nos(is), qos(1,is))
        if (js .ne. is) then
          if (Node.eq.0) then
            stop 'initatom: Unexpected species index returned by atom'
          else
            stop
          endif
        endif
      enddo 

!     Create the new data structures

      if (Node.eq.0)   call transfer

 1000 CONTINUE                 !! Temporary point 

!
      call broadcast_basis

c Initialize atomic lists

      no = 0
      nokb = 0
      do ia = 1,na
        is = isa(ia)
        noa  = species(is)%norbs
        nkba = species(is)%nprojs
        no = no + noa
        nokb = nokb + nkba
      enddo

      if (Node.eq.0) then
          call dump2cdf
          call dump_ascii
       endif

      if (Node.eq.0) atmonly = fdf_defined('Atom-Setup-Only')
#ifdef MPI
      call MPI_Bcast(atmonly,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (atmonly) call bye("End of atom setup")
!
!      if (Node.eq.0) then
!        call dump_atom_old
!        call dump_atom_new
!      endif

      end subroutine initatom

      end module atminit




