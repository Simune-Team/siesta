      subroutine vlist( iOpt, NO, endC, listC, NP, endCt, listCt,
     .                  NVmax, numVs, listVs, MaxInd, ind )

C ********************************************************************
C Finds the list of non-zero matrix elements of the potential.
C Written by P.Ordejon.
C Name and interfase modified by J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer iOpt      : Option switch: iOpt=0 => Check list
C                                    iOpt=1 => Create list
C integer NO        : Number of rows in C (num. of basis orbitals)
C integer endC(NO)  : Accumulated umber of nonzero elements in
C                     each row of C
C integer listC(*)  : List of nonzero elements in each row of C
C integer NP        : Number of columns in C (num. of mesh points)
C integer endCt(NP) : Accumulated number of nonzero elements in
C                     each column of C
C integer listCt(*) : List of nonzero elements in each column of C
C integer NVmax     : First dimension of listV and Vs, and maxim.
C                     number of nonzero elements in any row of Vs
C integer MaxInd    : Dimension of auxiliary array ind
C ********************* INPUT or OUTPUT (depending of iOpt) ***********
C integer numVs(NO)       : Number of nonzero elements in each row of Vs
C integer listVs(NVmax,NO): List of nonzero elements in each row of Vs
C ********************* AUXILIARY *************************************
C integer ind(MaxInd)   : Space that can be used freely between calls
C                         MaxInd must be at least NO
C *********************************************************************

      implicit none

      integer
     .   iOpt, NO, NP, NVmax, MaxInd,
     .   endC(0:NO),  listC(*),
     .   endCt(0:NP), listCt(*),
     .   numVs(NO), listVs(NVmax,NO), ind(NO)

      integer i, in, j, K, kn, n, nmax
      logical baddim

C     Check size of auxiliary array
      if (NO .gt. MaxInd) stop 'VLIST: dimension MaxInd too small'

C     Full initialization of array ind done only once
      do 10 i = 1, NO
        ind(i) = 0
   10 continue

C     Loop on orbitals
      baddim = .false.
      nmax = 0
      do 60 i = 1, NO

C       Copy listVs to array 'ind' in full-line format
        if (iOpt .eq. 0) then
          do 20 in = 1, numVs(i)
            j = listVs(in,i)
            ind(j) = 1
   20     enddo
        endif
        
C       Loop on mesh points of orbital i
        n = 0
        do 40 in = 1+endC(i-1), endC(i)
          K = listC(in)

C         Loop on other non-zero orbitals in mesh point
          do 30 kn = 1+endCt(K-1), endCt(K)
            j = listCt(kn)
C           If overlapping orbital j is not yet included in listVs
            if (ind(j) .eq. 0) then
              if (iOpt .eq. 0) then
                write(6,*) 
     .             'VLIST: List inconsistency. io,jo,kp =', i, j, K
                stop
              else
C               Switch-on orbital j in full-row vector ind
              	  ind(j) = 1
                n = n+1
C               Add orbital j to listVs, in sparse format
                if (n .le. NVmax) then
                  listVs(n,i) = j
                else
                  baddim = .true.
                endif
              endif
            endif
   30     enddo

   40   enddo
        if (iOpt .ne. 0) numVs(i) = n
        nmax = max( n, nmax )

C       Restore ind(j) for next orbital i
        do 50 in = 1, numVs(i)
          j = listVs(in,i)
          ind(j) = 0
   50   enddo

   60 enddo

C     Write correct value for dimension NVmax, if it was too small.
      if (baddim) then
        write(6,*) 'VLIST: Dimension too small. Make NVmax =', nmax
        stop
      endif

      end

