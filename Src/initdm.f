      subroutine initdm(Datm, Dscf, Dold, Escf, lasto, maxa,
     .                  maxno, maxo, maxspn, na, no, nspin, 
     .                  numh, numhold, listh, listhold, iaorb, 
     .                  found, inspn, usesavedm)

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
c ********* INPUT ***************************************************
c logical usesavedm     : whether DM has to be read from files or not
c logical found         : whether DM was found in files
c logical inspn         : true : AF ordering according to atom ordering
c                                if no DM files, no DM.InitSpin, ispin=2
c                         false: Ferro ordering  (fdf DM.InitSpinAF)
c integer na            : Number of atoms
c integer no            : Number of orbitals
c integer nspin         : Number of spins (1 or 2)
c integer maxa          : Max num. atoms for dimension
c integer maxo          : Max. number of orbitals
c integer maxno         : Max number of nonzero interactions
c integer maxspn        : Max number of spins (1 or 2) for dimensions
c integer lasto(0:maxa) : List with last orbital of each atom
c integer numh(maxo)        : Dscf matrix sparse information
c integer listh(maxno,maxo) :  "
c integer numhold(maxo)        : Same for Dold
c integer listhold(maxno,maxo) :  "
c integer iaorb(maxo)   : List saying to what atom an orbital belongs 
c double Datm(NO)       : Occupations of basis orbitals in free atom
c ********* OUTPUT **************************************************
c double Dscf(maxno,maxo,maxspn) : Density matrix in sparse form
c double Dold(maxno,maxo,maxspn) : same Dscf (for previous cycle)
c double Escf(maxno,maxo,maxspn) : Energy-density matrix in sparse form
c *******************************************************************

      implicit          none
      logical           found, inspn, usesavedm
      integer           no, na, maxo, maxno, nspin, maxa, maxspn
      integer           lasto(0:maxa), numh(maxo), numhold(maxo),
     .                  listh(maxno,maxo), listhold(maxno,maxo),
     .                  iaorb(maxo)
      double precision  Dscf(maxno,maxo,maxspn), 
     .                  Dold(maxno,maxo,maxspn), 
     .                  Escf(maxno,maxo,maxspn),
     .                  Datm(maxo)
      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         updo*1, line*130, names*80
c     character         fname*24
      logical           peratm
      integer           maxat, ni, nn, nr, nv, iat, nat, ia, iu,
     .                  i1, i2, in, ispin, jo, io

      parameter         (maxat = 500)
      integer           integs(4), lastc, lc(0:3), atom(maxat)
      double precision  spin(maxat), spinat, reals(4), values(4),
     .                  epsilon, qio, spio, aspin, rate


c enable FDF input/output

      include 'fdf/fdfdefs.h'

      data epsilon / 1.d-8 /

c try to read DM from disk if wanted (DM.UseSaveDM true) ---------------------

      if (usesavedm) then
        call iodm( 'read', maxno, maxo, no, nspin,
     .             numhold, listhold, Dscf, found )
      else
        found = .false.
      endif

c if found update Dold, otherwise initialize with neutral atoms

      if (found) then

        do ispin = 1,nspin
          do io = 1,no
            do in = 1,numhold(io)
              Dold(in,io,ispin) = Dscf(in,io,ispin)
            enddo
          enddo
        enddo

      else

c see whether specific initial spins are given through a DM.InitSpin block
c and read them eventually in a loop on atoms where lines are read and parsed
c   integer nat       : how many atoms to polarize
c   integer atom(nat) : which atoms
c   double  spin(nat) : what polarization ------------------------------------
 
        peratm = fdf_block('DM.InitSpin',iu)

        if ( peratm .and. (nspin .eq. 2) ) then

          nat = 0
          do iat = 1, na+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if ((ni .ge. 1) .and. (nr .ge. 1)) then
              nat = nat + 1
              if (nat .gt. maxat) then
                 write(6,"(/,a)") 
     .                'initdm: ERROR: BAD DIMENSIONS. Too small maxat'
                 stop 'initdm: ERROR: BAD DIMENSIONS. Too small maxat'
              endif
              atom(nat) = integs(1)
              spin(nat) = reals(1)
              if ( atom(nat) .gt. na ) then
                write(6,'(/,a)')
     .               'intdm: ERROR: Wrong entry in DM.InitSpin'
                stop 'intdm: ERROR: Wrong entry in DM.InitSpin'
              endif
            else if ( (ni .ge. 1) .and. (nn .ge. 1) ) then
              updo = names(lc(0)+1:lc(1))
              if ((updo .eq. '+') .or. (updo .eq. '-')) then
                nat = nat +1
                if (nat .gt. maxat) then
                   write(6,"(/,a)") 
     .                  'initdm: ERROR: BAD DIMENSIONS. Too small maxat'
                   stop 'initdm: ERROR: BAD DIMENSIONS. Too small maxat'
                endif
                atom(nat) = integs(1)
                if ( atom(nat) .gt. na ) then
                  write(6,'(/,a)')
     .                 'intdm: ERROR: Wrong entry in DM.InitSpin'
                  stop 'intdm: ERROR: Wrong entry in DM.InitSpin'
                endif
              else
                write(6,'(/,a)') 
     .               'intdm: ERROR: Wrong entry in DM.InitSpin'
                stop 'intdm: ERROR: Wrong entry in DM.InitSpin'
              endif
              if (updo .eq. '+') spin(nat) =  100.d0
              if (updo .eq. '-') spin(nat) = -100.d0
            else
              goto 50
            endif
          enddo

          write(6,'(a)') 
     .         'initdm: ERROR: Too many atom entries in DM.InitSpin'
          stop 'initdm: ERROR: Too many atom entries in DM.InitSpin'

   50     continue

c initialize to 0

          do io = 1,no
            do in = 1,numh(io)
              do ispin = 1,nspin
                Dscf(in,io,ispin) = 0.d0
                Escf(in,io,ispin) = 0.d0
              enddo
            enddo
          enddo

c initialize all paramagnetic 

          do ia = 1, na
            do io = lasto(ia-1) + 1, lasto(ia)
              do in = 1, numh(io)
                jo = listh(in,io)
                if (io .eq. jo) then
                  Dscf(in,io,1) = 0.5d0 * Datm(io)
                  Dscf(in,io,2) = Dscf(in,io,1)
                  Dold(in,io,1) = Dscf(in,io,1)
                  Dold(in,io,2) = Dscf(in,io,2)
                endif
              enddo
            enddo
          enddo

c loop on atoms with spin

          do iat = 1, nat
            ia = atom(iat)

c find maximum atomic moment that the atoms involved can carry
          
            spinat = 0.d0
            do io = lasto(ia-1) + 1, lasto(ia)
              spinat = spinat + min( Datm(io), 2.d0 - Datm(io) )
            enddo

c if given spin is larger than possible, make it to max atomic

            aspin = abs(spin(iat))
            if ((aspin .gt. spinat) .and. (aspin .gt. epsilon)) 
     .         spin(iat) = spinat*spin(iat)/aspin 

c initialize orbitals with same rate as atom

            rate = spin(iat) / spinat
            do io = lasto(ia-1) + 1, lasto(ia)
              qio = Datm(io)
              spio = rate * min( Datm(io), 2.d0 - Datm(io) )
              do in = 1, numh(io)
                jo = listh(in,io)
                if (io .eq. jo) then
                  Dscf(in,io,1) = 0.5d0 * (qio + spio)
                  Dscf(in,io,2) = 0.5d0 * (qio - spio)
                  Dold(in,io,1) = Dscf(in,io,1)
                  Dold(in,io,2) = Dscf(in,io,2)
                endif
              enddo
            enddo

          enddo

c ----------------------------------------------------------------------------

        else

c automatic, for non magnetic (nspin=1) or for Ferro or Antiferro ------------

          do io = 1, no
            do in = 1, numh(io)
              do ispin = 1, nspin
                Dscf(in,io,ispin) = 0.d0
                Escf(in,io,ispin) = 0.d0
              enddo
            enddo
            do in = 1,numh(io)
              jo = listh(in,io)
              if (io .eq. jo) then
                if (nspin .eq. 1) then

c No spin polarization

                  Dscf(in,io,1) = Datm(io)
                  Dold(in,io,1) = Datm(io)
                else

c Spin polarization

                  i1 = 1
                  i2 = 2

C ferro or antiferro according to DM.InitSpinAF (inspn)

                  if(inspn) then
                    if (mod(iaorb(io),2).eq.0) then
                      i1 = 2
                      i2 = 1
                    endif
                  endif
                  Dscf(in,io,i1) = min( Datm(io), 1.d0 )
                  Dscf(in,io,i2) = Datm(io) - Dscf(in,io,i1)
                  Dold(in,io,i1) = Dscf(in,io,i1)
                  Dold(in,io,i2) = Dscf(in,io,i2)
                endif
              endif
            enddo
          enddo

        endif

      endif

      return
      end

