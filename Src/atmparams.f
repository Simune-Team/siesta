      module atmparams

      implicit none 
!
!    Hard-wired parameters to be eliminated in the future
!
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
C                      with the same angular  momentum and 
C                      for the same specie.       
C INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
C                      for each angular momentum
C INTEGER  NSMX    : Maximum number of semicore shells for each angular
C                    momentum present in the atom ( for normal atom nsmx=0)
C INTEGER  NRMAX    : Maximum number of points in the functions read
C                      from file '.vps' or '.psatom.data' (this number is
C                      determined by the parameter nrmax in the
C                      program atm, which generates the files with
C                      the pseudopotential information).

       integer  ntbmax, lmaxd, nzetmx, nkbmx, nrmax, lmx2, nsemx, nsmx
         parameter ( ntbmax =  500 )
         parameter ( lmaxd  =    4 )
         parameter ( nzetmx =    3 )
         parameter ( nkbmx  =    2 )
         parameter ( nsmx  =    1 )
         parameter ( nsemx = 1 + nsmx) 
         parameter ( nrmax  = 2000 )
         parameter ( lmx2   = (lmaxd+1)*(lmaxd+1) )

      integer, parameter      :: maxos=2*nzetmx*lmx2*nsemx

      end module atmparams
