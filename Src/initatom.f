      subroutine initatom(ns, na, izs, smass, lmxkbs, 
     .                    lmaxs, maxl, lmax, 
     .                    nzls, maxos, maxzet, zetmax, maxo, maxkb,
     .                    rcls, contrf, isa, atm_label,
     .                    overflow, omax, osmax, kbmax,
     .                    nkbs, nos, qos,
     .                    no, nokb, qtot, rmaxv, rmaxo, rmaxkb,
     .                    lasto, lastkb, iza, amass,
     .                    iaorb, iphorb, Datm, qa, iaKB, iphKB)
C *********************************************************************
C Routine to initialize the Pseudopotentials and Atomic Orbitals,
C and the atomic lists.
C
C Writen by J.Soler and P.Ordejon, August-October'96
C **************************** INPUT **********************************
C integer ns                : Number of species
C integer na                : Number of atoms
C integer izs(ns)           : Atomic number of each species
C real*8 smass(ns)          : Atomic mass of each species
C integer maxl              : Maximum angular momentum for orbitals
C integer maxos             : Maximum num. of basis orbs of any atom
C integer maxzet            : Maximum number of orbital zetas
C integer maxo              : Maximum number of orbitals
C integer maxkb             : Maximum number of KB projectors
C integer isa(na)           : Species index of each atom
C character*20 atm_label(ns): Label for the atomic-files (i.e. pseudopot.
C                              files, input/output basis-sets files, etc..)
C **************************** OUTPUT *********************************
C integer lmxkbs(ns)        : L cutoff for KB projectors for each species
C integer lmaxs(ns)         : L cutoff for orbitals for each species
C integer nzls(0:maxl,ns)   : Number of zetas per L, for each species
C real*8 rcls(maxzet,0:maxl,ns) 
C                           : Cutoff radius for each atomic orbital
C real*8 contrf(maxzet,0:maxl,ns)
C                           : Scaling factor for each atomic orbital
C logical overflow          : True if dimensions are too small
C integer omax              : Actual value needed for maxo
C integer osmax             : Actual value needed for maxos
C integer kbmax             : Actual value needed for maxkb
C integer nkbs(ns)          : Tot num of KB projectors for each species
C integer nos(ns)           : Tot num of orbitals for each species
C real*8 qos(maxos,ns)      : Neutral atom occup. of each orbital
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
     .  maxkb, maxl, maxo, maxos, maxzet, na, no, nokb, ns,
     .  kbmax, omax, osmax, lmax, zetmax

      integer
     .  iaKB(*), iaorb(*), iphKB(*), iphorb(*), 
     .  isa(na), iza(na), izs(ns), 
     .  lastkb(0:na), lasto(0:na), lmaxs(ns), 
     .  lmxkbs(ns), nkbs(ns), nos(ns),
     .  nzls(0:maxl,ns)

      double precision
     .  amass(na), contrf(maxzet,0:maxl,ns), Datm(*), epskb,
     .  qa(na), qos(maxos,ns), qtot,
     .  rcls(maxzet,0:maxl,ns), rcut, rmaxv, rmaxo, rmaxkb, 
     .  smass(ns)

      logical
     .  overflow

      character 
     .   atm_label(ns)*20

      external
     .  atom, chkdime, epskb, rcut

C Internal variables ...................................................

      integer
     .  ia, io, is, js, nkba, noa, maxzetout
C ...................

      overflow = .false.
      kbmax = 1
      omax = 1
      osmax = 1

c Initialize pseudopotentials and atomic orbitals
      do is = 1,ns
        call atom( ns, izs(is), lmxkbs(is), lmaxs(is), maxl,
     .             nzls(0,is), maxzet, maxzetout, rcls(1,0,is), 
     .             contrf(1,0,is),
     .             js, atm_label(is),  nkbs(is), nos(is), 
     .             qos(1,is), maxos)
        if (js .ne. is)
     .    stop 'initatom: Unexpected species index returned by atom'
        call chkdime(maxl, lmaxs(is), overflow, lmax)
        call chkdime(maxos, nos(is), overflow, osmax )
        call chkdime(maxzet,maxzetout, overflow, zetmax)
      enddo 
      if (overflow) return

c Print some information on atomic orbitals
      if (.not.overflow) then
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
      endif

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

