      SUBROUTINE PDOSK( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH, 
     .                  MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .                  E1, E2, NHIST, SIGMA, 
     .                  XIJ, INDXUO, NK, KPOINT, WK, EO, 
     .                  HAUX, SAUX, PSI, AUX, DTOT, DPR, NUOTOT,
     .                  Node, Nodes )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C ****  INPUT  *********************************************************
C INTEGER NSPIN             : Number of spin components (1 or 2)
C INTEGER NUO               : Number of atomic orbitals in the unit cell
C INTEGER NO                : Number of atomic orbitals in the supercell
C INTEGER MAXSPN            : Second dimension of eo and qo 
C                             (maximum number of differents spin polarizations)
C INTEGER MAXUO             : Maximum number of atomic orbitals in the unit cell
C INTEGER MAXNH             : Maximum number of orbitals interacting
C                             with any orbital
C INTEGER MAXO              : First dimension of eo
C INTEGER NUMH(NUO)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C INTEGER LISTHPTR(NUO)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C INTEGER LISTH(MAXNH)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C REAL*8  H(MAXNH,NSPIN)    : Hamiltonian in sparse format
C REAL*8  S(MAXNH)          : Overlap in sparse format
C REAL*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C INTEGER NHIST             : Number of the subdivisions of the histogram
C REAL*8  SIGMA             : Width of the gaussian to expand the eigenvectors
C REAL*8  XIJ(3,MAXNH)      : Vectors between orbital centers (sparse)
C                             (not used if only gamma point)
C INTEGER INDXUO(NO)        : Index of equivalent orbital in unit cell
C INTEGER NK                : Number of k points
C REAL*8  KPOINT(3,NK)      : k point vectors
C REAL*8  WK(NK)            : Weights for k points
C REAL*8  EO(MAXO,MAXSPN,NK): Eigenvalues
C INTEGER NUOTOT            : Total number of orbitals per unit cell
C integer Node              : Local node number
C integer Nodes             : Total number of nodes
C ****  AUXILIARY  *****************************************************
C REAL*8  HAUX(2,NUO,NUO)   : Auxiliary space for the hamiltonian matrix
C REAL*8  SAUX(2,NUO,NUO)   : Auxiliary space for the overlap matrix
C REAL*8  PSI(2,NUO,NUO)    : Auxiliary space for the eigenvectors
C REAL*8  AUX(2*NUO*5)      : Extra auxiliary space
C ****  OUTPUT  ********************************************************
C REAL*8  DTOT(NHIST,2)   : Total density of states
C REAL*8  DPR(NHIST,NUO,2): Proyected density of states
C **********************************************************************

      use precision
      use parallel
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT NONE

      INTEGER
     .  NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH, NK, 
     .  MAXO, NHIST, NUOTOT, Node, Nodes

      INTEGER
     .  NUMH(NUO), LISTHPTR(NUO), LISTH(MAXNH),
     .  INDXUO(NO)

      DOUBLE PRECISION
     .  H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA, 
     .  XIJ(3,MAXNH), KPOINT(3,NK), EO(MAXO,MAXSPN,NK),
     .  HAUX(2,NUOTOT,NUO), SAUX(2,NUOTOT,NUO), PSI(2,NUOTOT,NUO),
     .  AUX(2*NUOTOT*5), DTOT(NHIST,2), DPR(NHIST,NUOTOT,2), WK(NK)

C Internal variables ---------------------------------------------------
      INTEGER
     .  IK, ISPIN, IUO, JUO, IO, J, JO, IHIST, IBAND, JOTY, IND,
     .  IERROR, IBANDG

      DOUBLE PRECISION
     .  KXIJ, CKXIJ, SKXIJ, DELTA, ENER, DIFF, PI, PIPJ1, PIPJ2, 
     .  PIPJS1, PIPJS2, GAUSS, NORM, WKSUM

#ifdef MPI
      integer ::
     .  BNode, Bnuo, maxnuo, MPIerror
      real*8, dimension(:,:,:), allocatable, save ::
     .  Sloc
#endif

      EXTERNAL
     .  CDIAG

C Initialize some variables
      DELTA = (E2 - E1)/NHIST
      PI = 4.0D0 * ATAN(1.0D0)

C Solve eigenvalue problem for each k-point
      DO ISPIN = 1, NSPIN

        DO IK = 1, NK

C Initialize auxiliary variables 
          DO IUO = 1, NUO
            DO JUO = 1, NUOTOT
              SAUX(1,JUO,IUO) = 0.D0
              SAUX(2,JUO,IUO) = 0.D0
              HAUX(1,JUO,IUO) = 0.D0
              HAUX(2,JUO,IUO) = 0.D0
            ENDDO
          ENDDO

          DO IUO = 1, NUO
            DO J = 1, NUMH(IUO)
              IND = LISTHPTR(IUO) + J
              JO = LISTH(IND)
              JUO= INDXUO(JO)
C Calculate the phases k*r_ij
              KXIJ = KPOINT(1,IK) * XIJ(1,IND) +
     .               KPOINT(2,IK) * XIJ(2,IND) +
     .               KPOINT(3,IK) * XIJ(3,IND) 
              CKXIJ = COS(KXIJ)
              SKXIJ = SIN(KXIJ)
C Calculate the Hamiltonian and the overlap in k space
C H(k) = Sum(R) exp(i*k*R) * H(R)
              SAUX(1,JUO,IUO) = SAUX(1,JUO,IUO) + S(IND) * CKXIJ
              SAUX(2,JUO,IUO) = SAUX(2,JUO,IUO) - S(IND) * SKXIJ
              HAUX(1,JUO,IUO) = HAUX(1,JUO,IUO) + H(IND,ISPIN) * CKXIJ
              HAUX(2,JUO,IUO) = HAUX(2,JUO,IUO) - H(IND,ISPIN) * SKXIJ
            ENDDO
          ENDDO

C Diagonalize for each k point
          CALL CDIAG( HAUX, NUOTOT, SAUX, NUOTOT, NUO,
     .                EO(1,ISPIN,IK), PSI, NUOTOT, NUOTOT, IERROR )

C Recalculate again the overlap matrix in k-space
          DO IUO = 1, NUO
            DO JUO = 1, NUOTOT
              SAUX(1,JUO,IUO) = 0.D0
              SAUX(2,JUO,IUO) = 0.D0
            ENDDO
          ENDDO

          DO IUO = 1, NUO
            DO  J = 1, NUMH(IUO)
              IND = LISTHPTR(IUO) + J
              JO = LISTH(IND)
              JUO= INDXUO(JO)
C Calculates the phases k*r_ij
              KXIJ = KPOINT(1,IK) * XIJ(1,IND) +
     .               KPOINT(2,IK) * XIJ(2,IND) +
     .               KPOINT(3,IK) * XIJ(3,IND) 
              CKXIJ = COS(KXIJ)
              SKXIJ = SIN(KXIJ)
C Calculates the hamiltonian and the overlap in k space
C H(k) = Sum(R) exp(i*k*R) * H(R)
              SAUX(1,JUO,IUO) = SAUX(1,JUO,IUO) + S(IND) * CKXIJ
              SAUX(2,JUO,IUO) = SAUX(2,JUO,IUO) - S(IND) * SKXIJ
            ENDDO
          ENDDO

#ifdef MPI
C Find maximum number of orbitals per node
          call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max,
     .      MPI_Comm_World,MPIerror)

C Allocate workspace array for broadcast overlap matrix
          allocate(Sloc(2,nuotot,maxnuo))
          call memory('A','D',2*nuotot*maxnuo,'pdosk')

C Loop over nodes broadcasting overlap matrix
          do BNode = 0,Nodes-1

C Find out how many orbitals there are on the broadcast node
            call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

C Transfer data
            if (Node.eq.BNode) then
              Sloc(1:2,1:nuotot,1:Bnuo) = Saux(1:2,1:nuotot,1:Bnuo)
            endif
#ifdef NODAT
            call MPI_Bcast(Sloc(1,1,1),2*nuotot*Bnuo,
     .        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#else
            call MPI_Bcast(Sloc(1,1,1),2*nuotot*Bnuo,
     .        DAT_double,BNode,MPI_Comm_World,MPIerror)
#endif

C Loop over all the energy range
            DO IHIST = 1, NHIST
              ENER = E1 + (IHIST - 1) * DELTA
              DO 170 IBAND = 1, NUO
                call LocalToGlobalOrb(IBAND,Node,Nodes,IBANDG)
                DIFF = (ENER - EO(IBANDG,ISPIN,IK))**2 / (SIGMA ** 2)
                IF (DIFF .GT. 15.0D0) THEN
                  GOTO 170
                ELSE
                  GAUSS = ( EXP(-DIFF) )
                  if (Node.eq.BNode) then
C Only add once to DTOT - not everytime loop over processors is executed
                    DTOT(IHIST,ISPIN) = DTOT(IHIST,ISPIN) + GAUSS*WK(IK)
                  endif
                  DO jo = 1, Bnuo
                    call LocalToGlobalOrb(jo,BNode,Nodes,juo)
                    DO IUO = 1, NUOTOT
C Solo para los Juo que satisfagan el criterio del record...
                      PIPJ1 = PSI(1,IUO,IBAND) * PSI(1,JUO,IBAND) +
     .                        PSI(2,IUO,IBAND) * PSI(2,JUO,IBAND)
                      PIPJ2 = PSI(1,IUO,IBAND) * PSI(2,JUO,IBAND) -
     .                        PSI(2,IUO,IBAND) * PSI(1,JUO,IBAND)
                      PIPJS1= PIPJ1*Sloc(1,IUO,JO)-PIPJ2*Sloc(2,IUO,JO)
                      PIPJS2= PIPJ1*Sloc(2,IUO,JO)+PIPJ2*Sloc(1,IUO,JO)
                      DPR(IHIST,JUO,ISPIN)= DPR(IHIST,JUO,ISPIN) + 
     .                                       PIPJS1*GAUSS*WK(IK)
                    ENDDO
                  ENDDO
                ENDIF
 170          ENDDO

            ENDDO

C End loop over broadcast nodes
          enddo

C Free workspace array for overlap
          call memory('D','D',size(Sloc),'pdosk')
          deallocate(Sloc)

#else
C Loop over all the energy range
          DO IHIST = 1, NHIST
            ENER = E1 + (IHIST - 1) * DELTA
            DO 170 IBAND = 1, NUO
              DIFF = (ENER - EO(IBAND,ISPIN,IK))**2 / (SIGMA ** 2)
              IF (DIFF .GT. 15.0D0) THEN
                GOTO 170
              ELSE
                GAUSS = ( EXP(-DIFF) )
                DTOT(IHIST,ISPIN) = DTOT(IHIST,ISPIN) + GAUSS*WK(IK)
                DO IUO = 1, NUOTOT
C Solo para los Juo que satisfagan el criterio del record...
                  DO JUO = 1, NUOTOT
                    PIPJ1 = PSI(1,IUO,IBAND) * PSI(1,JUO,IBAND) +
     .                      PSI(2,IUO,IBAND) * PSI(2,JUO,IBAND)
                    PIPJ2 = PSI(1,IUO,IBAND) * PSI(2,JUO,IBAND) -
     .                      PSI(2,IUO,IBAND) * PSI(1,JUO,IBAND)
                    PIPJS1= PIPJ1*SAUX(1,IUO,JUO)-PIPJ2*SAUX(2,IUO,JUO)
                    PIPJS2= PIPJ1*SAUX(2,IUO,JUO)+PIPJ2*SAUX(1,IUO,JUO)
                    DPR(IHIST,JUO,ISPIN)= DPR(IHIST,JUO,ISPIN) + 
     .                                     PIPJS1*GAUSS*WK(IK)
                  ENDDO
                ENDDO
              ENDIF
 170        ENDDO

          ENDDO
#endif

        enddo

      enddo

#ifdef MPI
C Allocate workspace array for global reduction
      allocate(Sloc(nhist,max(nuotot,nspin),nspin))
      call memory('A','D',nhist*nuotot*nspin,'pdosk')

C Global reduction of DPR matrix
      Sloc(1:nhist,1:nuotot,1:nspin) = 0.0d0
#ifdef NODAT
      call MPI_AllReduce(dpr(1,1,1),Sloc(1,1,1),nhist*nuotot*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(dpr(1,1,1),Sloc(1,1,1),nhist*nuotot*nspin,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif
      dpr(1:nhist,1:nuotot,1:nspin) = Sloc(1:nhist,1:nuotot,1:nspin)

C Global reduction of DTOT matrix
      Sloc(1:nhist,1:nspin,1) = 0.0d0
#ifdef NODAT
      call MPI_AllReduce(dtot(1,1),Sloc(1,1,1),nhist*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(dtot(1,1),Sloc(1,1,1),nhist*nspin,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif
      dtot(1:nhist,1:nspin) = Sloc(1:nhist,1:nspin,1)

C Free workspace array for global reduction
      call memory('D','D',size(Sloc),'pdosk')
      deallocate(Sloc)
#endif

      WKSUM = 0.0d0
      DO IK = 1,NK
        WKSUM = WKSUM + WK(IK)
      ENDDO

      NORM = SIGMA * SQRT(PI) * NK * WKSUM

      DO IHIST = 1, NHIST
        DO ISPIN = 1, NSPIN
          DTOT(IHIST,ISPIN) = DTOT(IHIST,ISPIN) / NORM
          DO IUO = 1, NUOTOT
            DPR(IHIST,IUO,ISPIN) = DPR(IHIST,IUO,ISPIN) /NORM
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
