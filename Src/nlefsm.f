C $Id: nlefsm.f,v 1.10 1999/05/05 17:25:35 emilio Exp $

      subroutine nlefsm( scell, nua, na, isa, xa, indxua,
     .                   nomax, nnomax, ndmax,
     .                   lasto, lastkb, iphorb, iphKB,
     .                   numd, listd, numh, listh, nspin, Dscf, 
     .                   Enl, fa, stress, H )
C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, June 1997.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C integer nomax            : Maximum number of atomic orbitals
C integer nnomax           : Maximum number of basis orbitals interacting
C                            either directly or thru a KB projector
C integer ndmax            : Maximum number of elements of each row of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(no)         : Number of nonzero elements of each row of the
C                            density matrix
C integer listd(ndmax,no)  : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(no)         : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listh(nnomax,no) : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin components
C real*8  Dscf(ndmax,nomax,nspin) : Density matrix
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(nnomax,nomax,nspin) : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
      implicit none

      integer
     . na, ndmax, nnomax, nomax, nspin, nua

      integer
     . indxua(na), iphKB(*), iphorb(*), isa(na),  
     . lasto(0:na), lastkb(0:na), listd(ndmax,*), listh(nnomax,*),
     . numd(*), numh(*)

      double precision
     . scell(3,3), Dscf(ndmax,nomax,nspin), Enl, epskb,
     . fa(3,nua), H(nnomax,nomax,nspin),
     . rcut, stress(3,3), xa(3,na), volcel

      external
     . chkdim, epskb, matel, neighb, rcut, timer, volcel

C Internal variables ................................................
C maxkba = maximum number of KB projectors of one atom
C maxna  = maximum number of neighbour atoms of one atom
C maxno  = maximum number of basis orbitals overlapping a KB projector
C maxo   = maximum total number of basis orbitals
      integer maxkba, maxna, maxno, maxo
      parameter ( maxkba =    16 )
      parameter ( maxna  =  1000 )
      parameter ( maxno  =  2000 )
      parameter ( maxo   = 20000 )
  
      integer
     .  ia, iano(maxno), ikb, ina, iana(maxna), ino,
     .  io, ioa, iono(maxno), is, ispin, ix, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no

      double precision
     .  Cijk, Di(maxo), epsk, fik, grSki(3,maxkba,maxno), 
     .  r2ki(maxna), rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, Ski(maxkba,maxno),
     .  Vi(maxo), volume, xki(3,maxna), xno(3,maxno)

      logical
     .  frstme, listed(maxo), within
     
      data frstme /.true./
C ......................

C     Start time counter
      call timer( 'nlefsm', 1 )

C     Print array sizes
      if (frstme) then
        call prmem( 0, 'nlefsm', 'iano',  'i', maxno          )
        call prmem( 0, 'nlefsm', 'iana',  'i', maxna          )
        call prmem( 0, 'nlefsm', 'grSki', 'd', 3*maxkba*maxno )
        call prmem( 0, 'nlefsm', 'r2ki',  'd', maxna          )
        call prmem( 0, 'nlefsm', 'Ski',   'd', maxkba*maxo    )
        call prmem( 0, 'nlefsm', 'Vi',    'd', maxo           )
        call prmem( 0, 'nlefsm', 'xki',   'd', 3*maxna        )
        call prmem( 0, 'nlefsm', 'xno',   'd', maxno          )
        call prmem( 0, 'nlefsm', ' ',     ' ', 0              )
        frstme = .false.
      endif
      
C     Find unit cell volume
      volume = volcel( scell ) * nua / na

C     Find maximum range
      rmaxo = 0.d0
      rmaxkb = 0.d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxkb

C     Initialize neighb subroutine
      nna = maxna
      call neighb( scell, rmax, na, xa, 0, 0,
     .             nna, iana, xki, r2ki )

C     Initialize arrays Di and Vi only once
      no = lasto(na)
      call chkdim( 'nlefsm', 'maxo', maxo, no, 1 )
      Enl = 0.d0
      do jo = 1,no
        Di(jo) = 0.d0
        Vi(jo) = 0.d0
        listed(jo) = .false.
      enddo

C     Loop on atoms with KB projectors      
      do ka = 1,na
        kua = indxua(ka)
        ks = isa(ka)
        nkb = lastkb(ka) - lastkb(ka-1)
        call chkdim( 'nlefsm', 'maxkba', maxkba, nkb, 1 )

C       Find neighbour atoms
        nna = maxna
        call neighb( scell, rmax, na, xa, ka, 0,
     .               nna, iana, xki, r2ki )
        call chkdim( 'nlefsm', 'maxna', maxna, nna, 1 )

C       Find neighbour orbitals
        nno = 0
        do ina = 1,nna
          ia = iana(ina)
          is = isa(ia)
          rki = sqrt(r2ki(ina))
          do io = lasto(ia-1)+1,lasto(ia)
            ioa = iphorb(io)
          
C           Find if orbital is within range
            within = .false.
            do ko = lastkb(ka-1)+1,lastkb(ka)
              koa = iphKB(ko)
              if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) within = .true.
            enddo

C           Find overlap between neighbour orbitals and KB projectors
            if (within) then
              nno = nno + 1
              call chkdim( 'nlefsm', 'maxno', maxno, nno, 1 )
              iono(nno) = io
              iano(nno) = ia
              do ix = 1,3
                xno(ix,nno) = xki(ix,ina)
              enddo
              ikb = 0
              do ko = lastkb(ka-1)+1,lastkb(ka)
                ikb = ikb + 1
                ioa = iphorb(io)
                koa = iphKB(ko)
                call matel( 'S', ks, is, koa, ioa, xki(1,ina),
     .                      Ski(ikb,nno), grSki(1,ikb,nno) )
              enddo
            endif
          enddo
        enddo
        
C       Loop on neighbour orbitals
        do ino = 1,nno
          io = iono(ino)
          ia = iano(ino)
          if (ia .le. nua) then

C           Scatter density matrix row of orbital io
            do j = 1,numd(io)
              jo = listd(j,io)
              do ispin = 1,nspin
                Di(jo) = Di(jo) + Dscf(j,io,ispin)
              enddo
            enddo
          
C           Scatter filter of desired matrix elements
            do j = 1,numh(io)
              jo = listh(j,io)
              listed(jo) = .true.
            enddo

C           Find matrix elements with other neighbour orbitals
            do jno = 1,nno
              jo = iono(jno)
              if (listed(jo)) then

C               Loop on KB projectors
                ikb = 0
                do ko = lastkb(ka-1)+1,lastkb(ka)
                  ikb = ikb + 1
                  koa = iphKB(ko)
                  epsk = epskb(ks,koa)
                  Sik = Ski(ikb,ino)
                  Sjk = Ski(ikb,jno)
                  Vi(jo) = Vi(jo) + epsk * Sik * Sjk
                  Cijk = Di(jo) * epsk
                  Enl = Enl + Cijk * Sik * Sjk
                  do ix = 1,3
                    fik = 2.d0 * Cijk * Sjk * grSki(ix,ikb,ino)
                    fa(ix,ia)  = fa(ix,ia)  - fik
                    fa(ix,kua) = fa(ix,kua) + fik
                    do jx = 1,3
                      stress(jx,ix) = stress(jx,ix) +
     .                                xno(jx,ino) * fik / volume
                    enddo
                  enddo
                enddo

              endif
            enddo

C           Pick up contributions to H and restore Di and Vi
            do j = 1,numh(io)
              jo = listh(j,io)
              do ispin = 1,nspin
                H(j,io,ispin) = H(j,io,ispin) + Vi(jo)
              enddo
              Vi(jo) = 0.d0
              listed(jo) = .false.
            enddo
            do j = 1,numd(io)
              jo = listd(j,io)
              Di(jo) = 0.d0
            enddo

          endif
        enddo
      enddo

      call timer( 'nlefsm', 2 )
      end

