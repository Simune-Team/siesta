      subroutine mixing(iscf,nbasis,maxo,ndmax,numd,listd,alpha,
     .                 dmnew,dmold,dmax)
C ***********************************************************************
C Subroutine to perform the mixing of the density matrix (D.M.) between
C SCF iterations. It also checks if the tolerance in the density
C matrix has been achieved, to stop the selfconsistency cicle.
C The density matrix is stored in sparse form. 
C
C Written by P.Ordejon, Noviembre'96
C ************************** INPUT **************************************
C integer iscf               : Current SCF iteration
C integer nbasis             : Number of atomic orbitals 
C integer maxo               : Maximum number of atomic orbitals (dimensio)
C integer ndmax              : First dimension of listd and D.M., and 
C                              maximum number of nonzero elements of 
C                              each row of the D.M.
C integer numd(maxo)         : Control vector of D.M.
C                              (number of nonzero elements of each row)
C integer listd(ndmax,maxo)  : Control vector of D.M.
C                              (list of nonzero elements of each row)
C real*8 alpha               : Mixing parameter
C ********************* INPUT AND OUTPUT*********************************
C real*8 dmnew(ndmax,maxo)   : Density Matrix
C                           Input: d.m. output in current SCF iter
C                           Output: d.m. input for next SCF interation
C real*8 dmold(ndmax,maxo)   : Density matrix
C                           Input: d.m. input in current SCF step
C                           Output: d.m. output in current SCF iter
C ************************** OUTPUT *************************************
C real*8 dmax                : Maximum change of a DM element between 
C                              input and output
C ************************ BEHAVIOUR ************************************
C Mixes the density matrix of the current and previous iterations
C like  n_new = (1 - alpha) * n_old + alpha * n_new
C For the first SCF cycle (iscf = 1), no mixing is done (alpha = 0.0),
C to erase density matrix from atom or former step, which derives from 
C nonorthogonal wavefunctions.
C ***********************************************************************
      implicit none

      integer
     .  iscf,maxo,nbasis,ndmax

      integer
     .  listd(ndmax,maxo),numd(maxo)

      double precision
     .  alpha,dmax,dmnew(ndmax,maxo),dmold(ndmax,maxo)

C Internal variables ........................................................
      integer
     .  i,in
C ........................

      dmax = 0.0
      do i = 1,nbasis
        do in = 1,numd(i)
          dmax = max( dmax, abs( dmnew(in,i) - dmold(in,i) ) )
          if (iscf .ne. 1) then
            dmnew(in,i) = (1.0-alpha)*dmold(in,i) + alpha*dmnew(in,i)
          endif
          dmold(in,i) = dmnew(in,i)
        enddo
      enddo
C .......................
	    
      return
      end



