      module atmparams

      implicit none 
!
!    Hard-wired parameters to be eliminated in the future
!

C INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
C                     with the same angular  momentum and 
C                     for the same species.       

         integer, parameter  :: nzetmx =    3  

C INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
C                     for each angular momentum

         integer, parameter  :: nkbmx  =    2  

C INTEGER  NSMX    : Maximum number of semicore shells for each angular
C                    momentum present in the atom ( for normal atom nsmx=0)

         integer, parameter  :: nsmx  =    2  
         integer, parameter  :: nsemx = 1 + nsmx  

C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                     orbitals,projectors and local neutral-atom 
C                     pseudopotential.

         integer, parameter  :: ntbmax =  500  

C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.

         integer, parameter  :: lmaxd  =    4  
         integer, parameter  :: lmx2   = (lmaxd+1)*(lmaxd+1)  

C INTEGER  NRMAX    : Maximum number of points in the functions read
C                     from file '.vps' or '.psatom.data' (this number is
C                     determined by the parameter nrmax in the
C                     program atm, which generates the files with
C                     the pseudopotential information).

         integer, parameter  :: nrmax  = 2000  

         integer, parameter      :: maxos=2*nzetmx*lmx2*nsemx

      end module atmparams
