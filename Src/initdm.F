      subroutine initdm(Datm, Dscf, Dold, lasto, maxa,
     .                  maxnh, maxo, maxuo, maxspn, nua, nuo, nspin, 
     .                  numh, numhold, listhptr, listhptrold,
     .                  listh, listhold, iaorb, found, inspn, 
     .                  usesavedm, nuotot)

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
      use parsing
#ifdef MPI
      use mpi_siesta
#endif
      use units, only : pi

      implicit          none

      logical           found, inspn, usesavedm
      integer           nuo, nua, maxo, maxuo, maxnh, nspin, maxa, 
     .                  maxspn
      integer           lasto(0:maxa), numh(maxuo), numhold(maxuo),
     .                  listhptr(maxuo), listhptrold(maxuo),
     .                  listh(maxnh), listhold(maxnh),
     .                  iaorb(maxo), nuotot
      real(dp)          Dscf(maxnh,maxspn), 
     .                  Dold(maxnh,maxspn), 
     .                  Datm(maxo)
      external          memory

c ---------------------------------------------------------------------

C Internal variables and arrays
 
      character         updo*1, line*130, names*80
      logical           noncol, peratm
      integer           nh, ni, nn, nr, nv, iat, nat, ia, iu,
     .                  i1, i2, in, ind, ispin, jo, io, 
     .                  iio, maxatnew

      integer, save ::  maxat

      integer           integs(4), lastc, lc(0:3)

      integer, dimension(:), allocatable, save ::
     .                  atom, ibuffer
#ifdef MPI
      integer
     .  MPIerror
#endif
      real(dp)          aspin, cosph, costh, epsilon, 
     .                  qio, rate, reals(4),
     .                  sinph, sinth, spinat, spio, values(4)

      real(dp), dimension(:), allocatable, save ::
     .                  phi, spin, theta, dpbuffer1

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
        if (Node.eq.0) then
          peratm = fdf_block('DM.InitSpin',iu)
          if (peratm .and. nspin.lt.2) write(6,'(/,a)')
     .    'initdm: WARNING: DM.InitSpin not used because nspin < 2'
        endif
#ifdef MPI
        call MPI_Bcast(peratm,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

        if (peratm .and. nspin.ge.2) then

C Allocate local memory
          allocate(atom(maxat))
          call memory('A','I',maxat,'initdm')
          allocate(phi(maxat))
          call memory('A','D',maxat,'initdm')
          allocate(spin(maxat))
          call memory('A','D',maxat,'initdm')
          allocate(theta(maxat))
          call memory('A','D',maxat,'initdm')

          nat = 0
          do iat = 1, nua+1
C Read and parse a line of the data block
            if (Node.eq.0) then
              read(iu,'(a)', end=50) line
            endif
#ifdef MPI
            call MPI_Bcast(line,130,MPI_character,0,MPI_Comm_World,
     .        MPIerror)
#endif
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (nn.ge.1 .and. names(lc(0)+1:lc(1)).eq.'%endblock') then
C End data reading
              goto 50
            elseif (ni .eq. 1) then
              if (nat.eq.maxat) then
                maxatnew = nat + nint(0.1*nat)
C
                allocate(ibuffer(maxat))
                call memory('A','I',maxat,'initdm')
                ibuffer(1:maxat) = atom(1:maxat)
                call memory('D','I',size(atom),'initdm')
                deallocate(atom)
                allocate(atom(maxatnew))
                call memory('A','I',maxatnew,'initdm')
                atom(1:maxat) = ibuffer(1:maxat)
                call memory('D','I',size(ibuffer),'initdm')
                deallocate(ibuffer)
C
                allocate(dpbuffer1(maxat))
                call memory('A','D',maxat,'initdm')
                dpbuffer1(1:maxat) = phi(1:maxat)
                call memory('D','D',size(phi),'initdm')
                deallocate(phi)
                allocate(phi(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                phi(1:maxat) = dpbuffer1(1:maxat)
                dpbuffer1(1:maxat) = spin(1:maxat)
                call memory('D','D',size(spin),'initdm')
                deallocate(spin)
                allocate(spin(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                spin(1:maxat) = dpbuffer1(1:maxat)
                dpbuffer1(1:maxat) = theta(1:maxat)
                call memory('D','D',size(theta),'initdm')
                deallocate(theta)
                allocate(theta(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                theta(1:maxat) = dpbuffer1(1:maxat)
                call memory('D','D',size(dpbuffer1),'initdm')
                deallocate(dpbuffer1)
C
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = integs(1)
            else
C Print bad-syntax error and stop
              goto 40
            endif
            if (nn .eq. 0) then
C Read value of spin
              if (nr .eq. 3) then
C Read spin value and direction
                spin(nat)  = reals(1)
                theta(nat) = reals(2) * pi/180.0d0
                phi(nat)   = reals(3) * pi/180.0d0
              elseif (nr .eq. 1) then
C Read spin value. Default direction.
                spin(nat)  = reals(1)
                theta(nat) = 0.d0
                phi(nat)   = 0.d0
              else
                goto 40
              endif
            else if (nn .eq. 1) then
C Read spin as + or - (maximun value)
              updo = names(lc(0)+1:lc(1))
              if (updo .eq. '+') then
                spin(nat) =  100.d0
              elseif (updo .eq. '-') then
                spin(nat) = -100.d0
              else
                goto 40
              endif
              if (nr .eq. 2) then
                theta(nat) = reals(1) * pi/180.0d0
                phi(nat)   = reals(2) * pi/180.0d0
              elseif (nr .eq. 0) then
                theta(nat) = 0.d0
                phi(nat)   = 0.d0
              else
                goto 40
              endif
            else
              goto 40
            endif
            if (atom(nat).lt.1 .or. atom(nat).gt.nua) then
              if (Node.eq.0) then
                write(6,'(/,a)')
     .          'intdm: ERROR: Bad atom index in DM.InitSpin, line', iat
                stop 'intdm: ERROR: Bad atom index in DM.InitSpin'
              else
                stop
              endif
            endif
            if (abs(theta(nat)).gt.1.d-12) noncol = .true.
          enddo

          if (Node.eq.0) then
            write(6,'(a)') 
     .         'initdm: ERROR: Too many atom entries in DM.InitSpin'
            stop 'initdm: ERROR: Too many atom entries in DM.InitSpin'
          else
            stop
          endif

   40     continue
          if (Node.eq.0) then
            write(6,*)
     .       'initdm: ERROR: bad syntax in DM.InitSpin, line', iat
            stop 'initdm: ERROR: bad syntax in DM.InitSpin'
          else
            stop
          endif

   50     continue

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
          call memory('D','I',size(atom),'initdm')
          deallocate(atom)
          call memory('D','D',size(phi),'initdm')
          deallocate(phi)
          call memory('D','D',size(spin),'initdm')
          deallocate(spin)
          call memory('D','D',size(theta),'initdm')
          deallocate(theta)

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
