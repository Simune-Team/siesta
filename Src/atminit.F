      module atminit

      use precision
      use atmfuncs, only: lmaxd, nsemx, nkbmx, nzetmx,
     $                    maxos, Node,
     $                    rcut

      implicit none

      private
      public initatom, initatomlists
!
!
      integer, allocatable, save :: izs(:), nkbs(:), nos(:)
      real*8, allocatable, save  :: smass(:)
      real*8, allocatable, save  :: qos(:,:)

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

      integer, intent(in)  :: ns, na, isa(na)
      integer, intent(out) :: no, nokb

C Internal variables ...................................................
      integer ia, is, js, nkba, noa
c     integer io

#ifdef MPI
      integer MPIerror
#endif
!
!     Automatic arrays
!
      logical           :: semic(ns)
      integer           :: lmxkbs(ns), lmaxs(ns)
      integer           :: lsemic(0:lmaxd,ns), nkbl(0:lmaxd,ns)
      integer           :: polorb(0:lmaxd,nsemx,ns)
      integer           :: nzls(0:lmaxd,nsemx,ns)
      real*8            :: erefkb(nkbmx,0:lmaxd,ns)
      real*8            :: charge(ns)
      real*8            :: contrf(nzetmx,0:lmaxd,nsemx,ns)
      real*8            :: rcls(nzetmx,0:lmaxd,nsemx,ns)

      character(len=10) :: basistype(ns)
      character(len=20) :: atomlabel(ns)

C Get Node number  --- Put it in atmfuncs
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Allocate memory --- deallocated in initatomlists

      allocate(qos(maxos,ns))
      call memory('A','D',maxos*ns,'atminit')
      allocate(izs(ns))
      call memory('A','I',ns,'atminit')
      allocate(nkbs(ns))
      call memory('A','I',ns,'atminit')
      allocate(nos(ns))
      call memory('A','I',ns,'atminit')
      allocate(smass(ns))
      call memory('A','D',ns,'atminit')

c Reading input for the pseudopotentials and atomic orbitals 

      if (Node.eq.0) then
        write(6,'(/2a)') 
     .    'initatom: Reading input for the pseudopotentials ',
     .    'and atomic orbitals'
      endif

      if (Node.eq.0) then
        call read_basis(ns,izs,lmxkbs,nkbl,erefkb,lmaxs,
     .    nzls,rcls,contrf,atomlabel,polorb,semic,lsemic,
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
      call MPI_Bcast(lsemic(0,1),(lmaxd+1)*ns,MPI_integer,0,
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
     .           atomlabel(is), polorb(0,1,is), semic(is), lsemic(0,is),
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

c Print some information on atomic orbitals
c       do is = 1,ns
c         write(6,'(/,a,2i4)') 'initatom: is,iz =', is, izs(is)
c         do io = 1,nos(is)
c           write(6,'(a,i4,2f12.6)') 'initatom: io,rcut,q   =',
c    .        io, rcut(is,io), qos(io,is)
c         enddo
c         write(6,'(a,i4,2f12.6)')   'initatom: io,rcut     =',
c    .       0, rcut(is,0)
c         do io = 1,nkbs(is)
c           write(6,'(a,i4,2f12.6)') 'initatom: io,rcut,eps =',
c    .       -io, rcut(is,-io), epskb(is,-io)
c         enddo
c       enddo

c Initialize atomic lists
      no = 0
      nokb = 0
      do ia = 1,na
        is = isa(ia)
        noa  = nos(is)
        nkba = nkbs(is)
        no = no + noa
        nokb = nokb + nkba
      enddo

      end subroutine initatom

      subroutine initatomlists(na, no, nokb, isa, qtot, 
     .                    rmaxv, rmaxo, rmaxkb, lasto, 
     .                    lastkb, iza, amass, iaorb,
     .                    iphorb, Datm, qa, iaKB, iphKB)
C *********************************************************************
C Routine to initialize the atomic lists.
C
C Created from initatom.F
C Writen by J.Soler and P.Ordejon, August-October'96 
C Strongly modified by D. Sanchez-Portal, Oct. 1998 
C Modified by DSP, July 1999.
C Modified by JDG, October 1999 for dynamic memory usage
C **************************** INPUT **********************************
C integer na                : Number of atoms
C integer no                : Total number of orbitals
C integer nokb              : Total number of KB projectors
C integer isa(na)           : Species index of each atom
C **************************** OUTPUT *********************************
C real*8 qtot               : Total number of electrons
C real*8 rmaxv              : Maximum cutoff for local potential Vna
C real*8 rmaxo              : Maximum cutoff for atomic orbitals
C real*8 rmaxkb             : Maximum cutoff for KB projectors
C integer lasto(0:na)       : Position of last orbital of each atom
C integer lastkb(0:na)      : Position of last KB proj. of each atom
C integer iza(na)           : Atomic number of each atom
C real*8 amass(na)          : Atomic mass of each atom
C integer iaorb(no)         : Atomic index of each orbital
C integer iphorb(no)        : Orbital index of each orbital in its atom
C real*8 Datm(no)           : Neutral atom charge of each orbital
C real*8 qa(na)             : Neutral atom charge of each atom
C integer iaKB(nokb)        : Atomic index of each KB projector
C integer iphKB(nokb)       : Index of each KB projector in its atom 
C                             (negative)
C *********************************************************************

      integer, intent(in)  ::  na, no, nokb, isa(na)

      integer, intent(out) ::
     .  iaKB(nokb), iaorb(no), iphKB(nokb), iphorb(no), 
     .  iza(na), lastkb(0:na), lasto(0:na)

      double precision, intent(out) ::
     .  amass(na), Datm(no), 
     .  qa(na), qtot,
     .  rmaxv, rmaxo, rmaxkb 

C Internal variables 
      integer  ia, io, is, nkba, noa, nol, nokbl

c Initialize atomic lists
      nol = 0
      nokbl = 0
      qtot = 0.d0
      rmaxv  = 0.d0
      rmaxo  = 0.d0
      rmaxkb = 0.d0
      lasto(0) = 0
      lastkb(0) = 0
      do ia = 1,na
        is = isa(ia)
        noa  = nos(is)
        nkba = nkbs(is)
        lasto(ia)  = lasto(ia-1)  + nos(is)
        lastkb(ia) = lastkb(ia-1) + nkbs(is)
        rmaxv = max( rmaxv, rcut(is,0) )
        iza(ia) = izs(is)
        amass(ia) = smass(is)
        qa(ia) = 0.0d0
        do io = 1,noa
          nol = nol + 1
          rmaxo = max( rmaxo, rcut(is,io) )
          iaorb(nol) = ia
          iphorb(nol) = io
          Datm(nol) = qos(io,is)
          qa(ia) = qa(ia) + qos(io,is)
          qtot = qtot + qos(io,is)
        enddo
        do io = 1,nkba
          nokbl = nokbl + 1
          rmaxkb = max( rmaxkb, rcut(is,-io) ) 
          iaKB(nokbl) = ia
          iphKB(nokbl) = -io
        enddo
      enddo

C Deallocate local memory

      call memory('D','D',size(qos),'atminit')
      deallocate(qos)
      call memory('D','I',size(izs),'atminit')
      deallocate(izs)
      call memory('D','I',size(nkbs),'atminit')
      deallocate(nkbs)
      call memory('D','I',size(nos),'atminit')
      deallocate(nos)
      call memory('D','D',size(smass),'atminit')
      deallocate(smass)

      end subroutine initatomlists

      end module atminit


