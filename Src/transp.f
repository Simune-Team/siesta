      subroutine transp( NO, endC, listC,
     .                   NP, NPF, NPxyz, endCt, listCt, CtoCt,
     .                   Node, Nodes )

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
C integer NP              : Number of columns stored locally
C integer NPF             : Number of columns in full system
C integer NPxyz(3)        : Dimensions of the grid
C *********************** OUTPUT **************************************
C integer endCt(0:NP)     : Accumulated number of nonzero elements in
C                           each column
C integer listCt(*)       : List of nonzero elements in each column
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C ********************************************************************

      use parallel

      implicit none
      integer  NO, NP, endC(0:NO),  listC(*), endCt(0:NP),
     .         listCt(*), CtoCt(*), NPF, NPxyz(3)
      integer  i, in, j, jj, n, Node, Nodes
      integer, dimension(:), allocatable, save :: numCt
      external memory

C     Allocate local memory
      allocate(numCt(NPF))
      call memory('A','I',NPF,'transp')

C     Find number of nonzero orbitals at each mesh point
      do J = 1, NPF
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
      do J = 1, NPF
        call GlobalToLocalMesh(j,NPxyz,Node,Nodes,jj)
        if (jj.gt.0) then
          endCt(jj) = endCt(jj-1) + numCt(J)
        endif
        numCt(J) = 0
      enddo

C     Make listCt and CtoCt lists
      do i = 1, NO
        do in = endC(i-1)+1, endC(i)
          J = listC(in)
          call GlobalToLocalMesh(j,NPxyz,Node,Nodes,jj)
          if (jj.gt.0) then
            numCt(jj) = numCt(jj) + 1
            n = endCt(jj-1) + numCt(jj)
            listCt(n) = i
            CtoCt(n) = in
*         write(6,'(a,i3,5i6)')
*    .     'transp: i,endC(i-1),in,ip,endCt(ip-1),kn=',
*    .      i,endC(i-1),in,J,endCt(J-1),n
          endif
        enddo
      enddo

C     Free local memory
      call memory('D','I',size(numCt),'transp')
      deallocate(numCt)

      end

