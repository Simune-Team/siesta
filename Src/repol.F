! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine repol(nmeshk,displ)

c *******************************************************************
c Reading k_space grids for the calculation of the macroscopic 
c polarization
c
c Reads fdf block. All three grids must be specified, therefore 
c a 3x3 matrix of integer numbers must be given:
C    For example: 
C        %block PolarizationGrids
C              100 11   6  yes 
C              10  120  7  no 
C              12  15   85 yes
C        %endblock PolarizationGrids  
C
c   *  First grid will be used to calculate the polarization
c      along the direction of the first lattice vector.
c   *  Second grid  will be used for the calculation along the 
C      the direction of the second lattice vector.
C   *  Third grid  will be used for the calculation along the 
C      the direction of the third lattice vector. 
C   *  The last column specifies if the bidimensional grids
C      are going to be diplaced from the origin or not as in the
C      Monkhorst-Pack algorithm (PHYS.REV. B13, 5188 (1976)).
C      This last column is optional
C
C   The number in the diagonal of the matrix specified the number
C   of points to be used in the one dimensional line integrals along 
C   the different directions. The other number specified the mesh 
c   used in the surface integrals. For more details see, 
c   R.D. King-Smith, and D. Vanderbilt, PRB 47, 1651 (1993).
c   If the number of point in one of the grids is zero, the calculation
C   will not be done for this particular direction. 
c      
c Written by D. Sanchez-Portal. March 1999
c ********* OUTPUT **************************************************
c integer   nmeshk(3,3)     : Grids for the integration in reciprocal 
c                             space
c real*8    displ(3)        : Displacement of the bidimensional 
c                             integration grid from the origin
c *******************************************************************

C
C  Modules
C
      use precision
      use parallel,  only : Node
      use fdf
      use parsing
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

C Passed variables and arrays
      integer           nmeshk(3,3)
      real(dp)          displ(3)

C Internal variables and arrays
 
      character         line*130, names*80, yiorno
      integer           ni, nn, nr, nv, igrd, iu, ix
      integer           integs(4), lastc, lc(0:3)
      real(dp)          reals(4), values(4)
#ifdef MPI
      integer           MPIerror
#endif

C Set outputs to zero
      do igrd = 1,3
        do ix = 1,3
          nmeshk(ix,igrd) = 0
        enddo 
        displ(igrd) = 0.0d0
      enddo

C Check for block and read

      if (Node.eq.0) then
        if ( fdf_block('PolarizationGrids',iu) ) then

          do igrd = 1, 3
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (ni .eq. 3) then
              nmeshk(1,igrd) = integs(1)
              nmeshk(2,igrd) = integs(2)
              nmeshk(3,igrd) = integs(3)
                 
              if (nn.ge.1) then 
                yiorno = names(lc(0)+1:lc(1)) 
                if ( yiorno.eq.'yes'  .or.
     .               yiorno.eq.'YES'  .or.
     .               yiorno.eq.'y'    .or.
     .               yiorno.eq.'Y'    .or.
     .               yiorno.eq.'true' .or.
     .               yiorno.eq.'TRUE' .or.
     .               yiorno.eq.'t'    .or.
     .               yiorno.eq.'T' ) then 
                        
                  displ(igrd) = 0.5d0
                             
                endif 
              endif         
            else
              goto 10
            endif

          enddo

        endif

  10    continue

      endif

  50  continue 

#ifdef MPI
C Globalise output
      call MPI_Bcast(displ,3,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(nmeshk(1,1),9,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif

      return
      end
