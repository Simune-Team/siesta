C $Id: initatom.f,v 1.11 1999/01/31 11:14:55 emilio Exp $

      subroutine initatom(ns, na, maxo, maxkb, isa, overflow,
     .                    no, nokb, qtot, rmaxv, rmaxo, rmaxkb,
     .                    lasto, lastkb, iza, amass,
     .                    iaorb, iphorb, Datm, qa, iaKB, iphKB)
C *********************************************************************
C Routine to initialize the Pseudopotentials and Atomic Orbitals,
C and the atomic lists.
C
C Writen by J.Soler and P.Ordejon, August-October'96 
C Strongly modified by D. Sanchez-Portal, Oct. 1998 
C **************************** INPUT **********************************
C integer ns                : Number of species
C integer na                : Number of atoms
C integer maxo              : Maximum number of orbitals
C integer maxkb             : Maximum number of KB projectors
C integer isa(na)           : Species index of each atom
C **************************** OUTPUT *********************************
C logical overflow          : True if dimensions are too small
C integer no                : Total number of orbitals
C integer nokb              : Total number of KB projectors
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
      implicit none

      integer
     .  maxkb, maxo, na, no, nokb, ns

      integer
     .  iaKB(*), iaorb(*), iphKB(*), iphorb(*), 
     .  isa(na), iza(na),
     .  lastkb(0:na), lasto(0:na)

      double precision
     .  amass(na), Datm(*), epskb,
     .  qa(na), qtot,
     .  rcut, rmaxv, rmaxo, rmaxkb 

      logical
     .  overflow

      external
     .  atom, chkdime, epskb, rcut

C Internal variables ...................................................
      include 'atom.h'
      integer
     .  ia, io, is, js, kbmax, nkba, noa, omax, maxos 
         parameter (maxos=2*nzetmx*lmx2)

      integer
     . lmxkbs(nsmax), lmaxs(nsmax), nzls(0:lmaxd,nsmax),
     . izs(nsmax), polorb(lmaxd,nsmax),
     . lsemic(nsmax),nkbs(nsmax),nos(nsmax)

      double precision
     . contrf(nzetmx,0:lmaxd,nsmax), rcls(nzetmx,0:lmaxd,nsmax),
     . smass(nsmax), charge(nsmax), qos(maxos,nsmax)

      logical
     .  semic(nsmax)

      character
     .   atm_label(nsmax)*20, basistype(nsmax)*10

C ...................

C ...................

      overflow = .false.
      kbmax = 1
      omax = 1

c Reading input for the pseudopotentials and atomic orbitals 

       write(6,'(/2a)') 
     .   'initatom: Reading input for the pseudopotentials ',
     .   'and atomic orbitals'

       call read_basis(ns,izs,lmxkbs,lmaxs,nzls,rcls,contrf,
     .   atm_label,polorb,semic,lsemic,charge,smass,basistype)

C .....................
c Initialize pseudopotentials and atomic orbitals 
      do is = 1,ns
        call atom( ns, izs(is), lmxkbs(is), lmaxs(is), 
     .             nzls(0,is), rcls(1,0,is), contrf(1,0,is), 
     .             atm_label(is), polorb(1,is), semic(is), lsemic(is),
     .             charge(is), smass(is), basistype(is), 
     .             js, nkbs(is), nos(is), qos(1,is))
        if (js .ne. is)
     .    stop 'initatom: Unexpected species index returned by atom'
      enddo 

c Print some information on atomic orbitals
        do is = 1,ns
c         write(6,'(/,a,2i4)') 'initatom: is,iz =', is, izs(is)
          do io = 1,nos(is)
c           write(6,'(a,i4,2f12.6)') 'initatom: io,rcut,q   =',
c    .        io, rcut(is,io), qos(io,is) 
          enddo
c         write(6,'(a,i4,2f12.6)')   'initatom: io,rcut     =',
c    .       0, rcut(is,0)
          do io = 1,nkbs(is)
c           write(6,'(a,i4,2f12.6)') 'initatom: io,rcut,eps =',
c    .       -io, rcut(is,-io), epskb(is,-io)
          enddo
        enddo

c Initialize atomic lists
      no = 0
      nokb = 0
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
        qa(ia) = 0.d0
c       if (ia .eq. 1) write(6,*) ' '
c       write(6,'(a,3i4,2i6)') 
c    .     'initatom: ia, is, iz, lasto, lastkb =',
c    .      ia, is, iza(ia), lasto(ia), lastkb(ia)
        do io = 1,noa
          no = no + 1
          call chkdime(maxo, no, overflow, omax )
          rmaxo = max( rmaxo, rcut(is,io) )
          if (.not.overflow) then
            iaorb(no) = ia
            iphorb(no) = io
            Datm(no) = qos(io,is)
            qa(ia) = qa(ia) + qos(io,is)
            qtot = qtot + qos(io,is)
          endif
        enddo
        do io = 1,nkba
          nokb = nokb + 1
          call chkdime(maxkb, nokb, overflow, kbmax )
            rmaxkb = max( rmaxkb, rcut(is,-io) ) 
          if (.not.overflow) then
            iaKB(nokb) = ia
            iphKB(nokb) = -io
          endif
        enddo
      enddo


      return
      end

