c
      subroutine compat_params(str)
c
c     Some internal parameters tend to change over time... This routine
c     assigns values to them based on a "compatibility string".
c  
c     This is a temporary kludge. Ideally, this kind of thing should be
c     done with fdf. However, the input file can contain references to
c     several concatenated calculations...
c
c     Currently, the only meaningful strings are:
c
c     'mons'
c     '  '    :  Use a denser grid up to larger radii.
c                Use a larger value for the ps's ecuts.
c                Use the Soler-Balbas XC package
c       
c     'ucb'   :  Revert to the standard circa 1990 UCB values. Note
c                that these correspond to the first released version
c                of Jose Luis Martins code, not to the old Froyen
c                version (that would be 'froyen' --to be implemented)
c
      character*(*) str

      include 'compat.h'

      logical leqi
      external leqi
c
      if (str .eq. ' ' .or. leqi(str,'mons')) then

         write(6,'(/,a,/)') ' *** MONS compatibility mode ***'
         aa_def = 6.d0
         bb_def = 80.d0
         rmax_def = 120.d0 
         ecuts = 1.d-3
         use_excorr = .false.

      else if (leqi(str,'ucb')) then

         write(6,'(/,a,/)') ' *** UCB compatibility mode ***'
         aa_def = 6.d0
         bb_def = 40.d0
         rmax_def = 80.d0
         ecuts = 1.2d-4
         use_excorr = .true.
c

      else

         stop 'COMPAT'

      endif

      return
      end

