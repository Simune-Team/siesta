      subroutine remass(smass,maxs,ns)

c *******************************************************************
c Reading atomic masses of different species.
c
c Reads fdf block. Not necessarily all species have to be given. The 
c ones not given at input will be assumed to have their natural mass
c (according to atmass subroutine).
c      
c Written by E. Artacho. March 1998.97. 
c ********* INPUT ***************************************************
c integer maxs             : Maximum number of species
c integer ns               : Number of species
c double smass(maxs)       : Atomic masses of each species (amu)
c                            as given by default by atmass
c ********* OUTPUT **************************************************
c double smass(maxs)       : Atomic masses of each species (amu)
c                            modified by explicit numbers at input
c *******************************************************************

      implicit          none
      integer           maxs, ns
      double precision  smass(maxs)

      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80
      integer           ni, nn, nr, nv, ist, nst, iu
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)


c enable FDF input/output

      include 'fdf/fdfdefs.h'

c check for block and read

      if ( fdf_block('AtomicMass',iu) ) then

         nst = 0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if ((ni .ge. 1) .and. (nr .ge. 1)) then
              nst = nst + 1
              if (nst .gt. maxs) then
                 write(6,"(/,a)") 
     .                'remass: ERROR: BAD DIMENSIONS. Too small maxs'
                 stop 'remass: ERROR: BAD DIMENSIONS. Too small maxs'
              endif
              if ((integs(1).gt.ns) .or. (integs(1).lt.1) ) then
                 write(6,"(/,a,i4)") 
     .           'remass: WARNING: Species out of range. Ignored is =',
     .           integs(1)
              else
                 smass(integs(1)) = reals(1)
                 write(6,"(a, i4, a, f12.5)")
     .            'remass: Read atomic mass for species ', integs(1),
     .            ' as ', reals(1)
              endif
            else
              return
            endif
         enddo

         write(6,'(a)') 
     .        'remass: ERROR: Too many entries in AtomicMass'
         stop 'remass: ERROR: Too many entries in AtomicMass'

      endif

 50   continue

      return
      end

