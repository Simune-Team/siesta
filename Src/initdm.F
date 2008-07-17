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
      subroutine initdm(Datm, Dscf, Dold, lasto, maxa,
     .                  maxnh, maxo, maxuo, maxspn, nua, nuo, nspin, 
     .                  numh, numhold, listhptr, listhptrold,
     .                  listh, listhold, iaorb, found, inspn, 
     .                  usesavedm )

c *******************************************************************
c Density matrix initialization
c
c    If UseSaveDM is true, it is read from file if present.
c    Otherwise it is generated assuming atomic charging 
c      (filling up atomic orbitals). The DM originated that way is
c      not a good DM due to overlaps, but the SCF cycling corrects
c      that for the next cycle.
c    Spin polarized calculations starting from atoms:
c      Default: All atoms with maximum polarization compatible with 
c               atomic configuration. In Ferromagnetic ordering (up).
c      If DM.InitSpinAF is true, as default but in Antiferro order:
c               even atoms have spin down, odd up.
c      If fdf %block DM.InitSpin is present it overwrites previous 
c         schemes: magnetic moments are explicitly given for some atoms.
c         Atoms not mentioned in the block are initialized non polarized.
c      
c Written by E. Artacho. December 1997. Taken from the original piece
c of siesta.f written by P. Ordejon.
c Non-collinear spin added by J.M.Soler, May 1998.
c ********* INPUT ***************************************************
c logical usesavedm     : whether DM has to be read from files or not
c logical found         : whether DM was found in files
c logical inspn         : true : AF ordering according to atom ordering
c                                if no DM files, no DM.InitSpin, ispin=2
c                         false: Ferro ordering  (fdf DM.InitSpinAF)
c integer nua           : Number of atoms in the unit cell
c integer nuo           : Number of orbitals in the unit cell
c integer nspin         : Number of spin components
c integer maxa          : Max num. atoms for dimension
c integer maxo          : Max. number of orbitals (globally)
c integer maxuo         : Max. number of orbitals (locally)
c integer maxnh         : Max number of nonzero interactions
c integer maxspn        : Max number of spin components in dimensions
c integer lasto(0:maxa) : List with last orbital of each atom
c integer numh(*)       : Dscf matrix sparse information
c integer listhptr(*)   :  "
c integer listh(maxnh)  :  "
c integer numhold(*)    : Same for Dold
c integer listhptrold(*):  "
c integer listhold(maxnh) :  "
c integer iaorb(maxo)   : List saying to what atom an orbital belongs 
c double Datm(no)       : Occupations of basis orbitals in free atom
c ********* OUTPUT **************************************************
c double Dscf(maxnh,maxspn) : Density matrix in sparse form
c double Dold(maxnh,maxspn) : same Dscf (for previous cycle)
c *******************************************************************

C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use fdf
      use sys,          only : die
      use units, only        : pi
      use alloc, only : re_alloc, de_alloc

      implicit          none

      logical           found, inspn, usesavedm
      integer           nuo, nua, maxo, maxuo, maxnh, nspin, maxa, 
     .                  maxspn
      integer           lasto(0:maxa), numh(maxuo), numhold(maxuo),
     .                  listhptr(maxuo), listhptrold(maxuo),
     .                  listh(maxnh), listhold(maxnh),
     .                  iaorb(maxo)
      real(dp)          Dscf(maxnh,maxspn), 
     .                  Dold(maxnh,maxspn), 
     .                  Datm(maxo)

c ---------------------------------------------------------------------

C Internal variables and arrays
 
      character         updo*1, msg*80
      logical           noncol, peratm, badsyntax
      integer           nh, ni, nn, nr, iat, nat, ia,
     .                  i1, i2, in, ind, ispin, jo, io, 
     .                  iio, maxatnew

      integer, save ::  maxat

      integer, pointer, save :: atom(:)
      real(dp)          aspin, cosph, costh, epsilon, 
     .                  qio, rate, sinph, sinth, spinat, spio

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp), pointer, save :: phi(:), spin(:), theta(:)

      data maxat / 1000 /
      data epsilon / 1.d-8 /

C Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------

      if (usesavedm) then
        call iodm( 'read', maxnh, nuo, nspin, numhold, 
     .             listhptrold, listhold, Dscf, found )
      else
        found = .false.
      endif

C If found update Dold, otherwise initialize with neutral atoms
      if (found) then

        if (nuo.gt.0) then
          nh = listhptrold(nuo) + numhold(nuo)
          Dold(1:nh,1:nspin) = Dscf(1:nh,1:nspin)
        endif

      else

C See whether specific initial spins are given in a DM.InitSpin block
C and read them in a loop on atoms where lines are read and parsed
C   integer nat       : how many atoms to polarize
C   integer atom(nat) : which atoms
C   double  spin(nat) : what polarization -----------------------------

        noncol = .false.
        peratm = fdf_block('DM.InitSpin',bfdf)
        if ((Node .eq. 0) .and. peratm .and. (nspin .lt. 2)) then
          write(6,'(/,a)')
     .      'initdm: WARNING: DM.InitSpin not used because nspin < 2'
        endif

        if (peratm .and. nspin.ge.2) then

C Allocate local memory
          nullify(atom,phi,spin,theta)
          call re_alloc( atom, 1, maxat, 'atom', 'initdm' )
          call re_alloc( phi, 1, maxat, 'phi', 'initdm' )
          call re_alloc( spin, 1, maxat, 'spin', 'initdm' )
          call re_alloc( theta, 1, maxat, 'theta', 'initdm' )

          nat = 0
          badsyntax = .FALSE.
          do while(fdf_bline(bfdf,pline) .and. (nat .lt. nua) .and.
     .            (.not. badsyntax))

            nn = fdf_bnnames(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

            if (ni .eq. 1) then
              if (nat .eq. maxat) then
                maxatnew = nat + nint(0.1*nat)
C
                call re_alloc( atom, 1, maxatnew, 'atom', 'initdm',
     &                         copy=.true. )
C
                call re_alloc( phi, 1, maxatnew, 'phi', 'initdm',
     &                         copy=.true. )
                call re_alloc( spin, 1, maxatnew, 'spin', 'initdm',
     &                         copy=.true. )
                call re_alloc( theta, 1, maxatnew, 'theta', 'initdm',
     &                         copy=.true. )
C
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = fdf_bintegers(pline,1)

              if (nn .eq. 0) then
C Read value of spin
                if (nr .eq. 3) then
C Read spin value and direction
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = fdf_breals(pline,2) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,3) * pi/180.0d0
                elseif (nr .eq. 1) then
C Read spin value. Default direction.
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              elseif (nn .eq. 1) then
C Read spin as + or - (maximun value)
                updo = fdf_bnames(pline,1)
                if (updo .eq. '+') then
                  spin(nat) =  100.d0
                elseif (updo .eq. '-') then
                  spin(nat) = -100.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
                if (nr .eq. 2) then
                  theta(nat) = fdf_breals(pline,1) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,2) * pi/180.0d0
                elseif (nr .eq. 0) then
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              else
C Print bad-syntax error and stop
                badsyntax = .TRUE.
              endif

              if ((atom(nat) .lt. 1) .or. (atom(nat) .gt. nua)) then
                write(msg,'(a,a,i4)') 'intdm: ERROR: Bad atom ' //
     .            'index in DM.InitSpin, line', nat+1
                call die(TRIM(msg))
              endif
              if (abs(theta(nat)) .gt. 1.d-12) noncol = .true.
            else
C Print bad-syntax error and stop
              badsyntax = .TRUE.
            endif
          enddo

          if (badsyntax) then
            write(msg,'(a,i4)')
     .        'initdm: ERROR: bad syntax in DM.InitSpin, line', nat+1
            call die(msg)
          endif

          if (nat .eq. nua) then
            call die('initdm: ERROR: Too many' //
     $               ' atom entries in DM.InitSpin')
          endif

          if (noncol .and. nspin.lt.4) then
            if (Node.eq.0) then
              write(6,'(/,2a)') 'initdm: WARNING: noncolinear spins ',
     .                 'in DM.InitSpin not used because nspin < 4'
            endif
            noncol = .false.
          endif

C Initialize to 0

          Dscf(1:maxnh,1:nspin) = 0.0d0

C Initialize all paramagnetic 

          do ia = 1, nua
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    Dscf(ind,1) = 0.5d0 * Datm(io)
                    Dscf(ind,2) = Dscf(ind,1)
                    Dold(ind,1) = Dscf(ind,1)
                    Dold(ind,2) = Dscf(ind,2)
                  endif
                enddo
              endif
            enddo
          enddo

C Loop on atoms with spin

          do iat = 1, nat
            ia = atom(iat)

C Find maximum atomic moment that the atoms involved can carry
          
            spinat = 0.d0
            do io = lasto(ia-1) + 1, lasto(ia)
              spinat = spinat + min( Datm(io), 2.d0 - Datm(io) )
            enddo
            if (spinat.lt.epsilon .and. Node.eq.0) print'(a,i6,a)', 
     .        'initdm: WARNING: atom ', atom(iat),
     .        ' has a closed-shell and cannot be polarized'

C If given spin is larger than possible, make it to max atomic

            aspin = abs(spin(iat))
            if ((aspin .gt. spinat) .and. (aspin .gt. epsilon)) 
     .         spin(iat) = spinat*spin(iat)/aspin 

C Initialize orbitals with same rate as atom

            rate = spin(iat) / (spinat+epsilon)
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                qio = Datm(io)
                spio = rate * min( Datm(io), 2.d0 - Datm(io) )
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    if (noncol) then
C Store non-collinear-spin density matrix as
C   ispin=1 => D11, ispin=2 => D22;
C   ispin=3 => Real(D12); ispin=4 => Imag(D12)
                      costh = cos(theta(iat))
                      sinth = sin(theta(iat))
                      cosph = cos(phi(iat))
                      sinph = sin(phi(iat))
                      Dscf(ind,1) = (qio + spio * costh) / 2
                      Dscf(ind,2) = (qio - spio * costh) / 2
                      Dscf(ind,3) =   spio * sinth * cosph / 2
                      Dscf(ind,4) = - spio * sinth * sinph / 2
                    else
                      Dscf(ind,1) = (qio + spio) / 2
                      Dscf(ind,2) = (qio - spio) / 2
                    endif
                    do ispin = 1,nspin
                      Dold(ind,ispin) = Dscf(ind,ispin)
                    enddo
                  endif
                enddo
              endif
            enddo

          enddo

C Deallocate local memory
          call de_alloc( atom, 'atom', 'initdm' )
          call de_alloc( phi, 'phi', 'initdm' )
          call de_alloc( spin, 'spin', 'initdm' )
          call de_alloc( theta, 'theta', 'initdm' )
c ---------------------------------------------------------------------

        else

C Initialize to 0
          Dscf(1:maxnh,1:nspin) = 0.0d0

C Automatic, for non magnetic (nspin=1) or for Ferro or Antiferro -----
          do io = 1, nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
              ind = listhptr(io)+in
              jo = listh(ind)
              if (iio .eq. jo) then
                if (nspin .eq. 1) then

C No spin polarization

                  Dscf(ind,1) = Datm(iio)
                  Dold(ind,1) = Datm(iio)
                else

C Spin polarization

                  i1 = 1
                  i2 = 2

C Ferro or antiferro according to DM.InitSpinAF (inspn)

                  if (inspn) then
                    if (mod(iaorb(iio),2).eq.0) then
                      i1 = 2
                      i2 = 1
                    endif
                  endif
                  Dscf(ind,i1) = min( Datm(iio), 1.d0 )
                  Dscf(ind,i2) = Datm(iio) - Dscf(ind,i1)
                  Dold(ind,i1) = Dscf(ind,i1)
                  Dold(ind,i2) = Dscf(ind,i2)
                endif
              endif
            enddo
          enddo

        endif

      endif

      return
      end
