C $Id: atom.h,v 1.3 1999/01/31 11:56:50 emilio Exp $

C Dimension parameters for routine atom:
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
C                      with the same angular  momentum and 
C                      for the same specie.       
C INTEGER  NRMAX    : Maximum number of points in the functions read 
C                      from file 'psatom.data' (this number is
C                      determined by the parameter nrmax in the 
C                      program atm, which generates the files with 
C                      the pseudopotential information).

       integer nsmax, ntbmax, lmaxd, nzetmx, nrmax, lmx2
         parameter ( nsmax  =    7 )
         parameter ( ntbmax =  500 )
         parameter ( lmaxd  =    4 )
         parameter ( nzetmx =    2 )
         parameter ( nrmax  = 2000 )
         parameter ( lmx2   = (lmaxd+1)*(lmaxd+1) )

