      subroutine transp( NO, endC, listC,
     .                   NP, endCt, listCt, CtoCt,
     .                   MaxNCt, numCt )

C ********************************************************************
C Finds the 'transpose' of a sparse matrix ordered by rows, or rather
C it constructs information to obtain elements of Ct (transpose of C)
C from the elements of C in sparse form (control vector C to Ct)
C Written by P.Ordejon.
C Name and interfase modified by J.M.Soler. April'95.
C *********************** INPUT **************************************
C integer NO              : Number of rows
C integer endC(0:NO)      : Accumulated umber of nonzero elements in
C                           each row
C integer listC(NCmax,NO) : List of nonzero elements in each row
C integer NP              : Number of columns
C integer MaxNCt          : Dimension of auxiliary array numCt
C *********************** OUTPUT **************************************
C integer endCt(0:NP)     : Accumulated number of nonzero elements in
C                           each column
C integer listCt(*)       : List of nonzero elements in each column
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C *********************** AUXILIARY **********************************
C integer numCt(MaxNCt)  : Space that can be used freely between calls
C                          MaxNCt must be at least NP
C ********************************************************************

      implicit none
      integer  NO, NP, MaxNCt, endC(0:NO),  listC(*), endCt(0:NP),
     .         listCt(*), CtoCt(*), numCt(NP)
      integer  i, in, J, n

C     Check size of auxiliary array
      if (NP .gt. MaxNCT) stop 'TRANSP: dimension MaxNCT too small'

C     Find number of nonzero orbitals at each mesh point
      do J = 1, NP
        numCt(J) = 0
      enddo
      do i = 1, NO
        do in = endC(i-1)+1, endC(i)
          J = listC(in)
          numCt(J) = numCt(J) + 1
        enddo
      enddo

C     Find limit of each point in arrays listCt and CtoCt
      endCt(0) = 0
      do J = 1, NP
        endCt(J) = endCt(J-1) + numCt(J)
        numCt(J) = 0
      enddo

C     Make listCt and CtoCt lists
      do i = 1, NO
        do in = endC(i-1)+1, endC(i)
          J = listC(in)
          numCt(J) = numCt(J) + 1
          n = endCt(J-1) + numCt(J)
          listCt(n) = i
          CtoCt(n) = in
*         write(6,'(a,i3,5i6)')
*    .     'transp: i,endC(i-1),in,ip,endCt(ip-1),kn=',
*    .      i,endC(i-1),in,J,endCt(J-1),n
        enddo
      enddo

      end

