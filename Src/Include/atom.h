C Dimension parameters for routine atom:
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of orbitals with the same angular
C                      momentum and for the same specie.       
C INTEGER  NRMAX    : Maximum number of points in the functions read 
C                      from file 'psatom.data' (this number is
C                      determined by the parameter nrmax in the 
C                      program atm, which generates the files with 
C                      the pseudopotential information).

         parameter ( nsmax  =    7 )
         parameter ( ntbmax =  100 )
         parameter ( lmaxd  =    4 )
         parameter ( nzetmx =    2 )
         parameter ( nrmax  = 2000 )

         parameter ( lmx2   = (lmaxd+1)*(lmaxd+1) )

