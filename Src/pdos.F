      SUBROUTINE PDOS( NO, NSPIN, MAXSPN, MAXUO, MAXNH, 
     .                 MAXO, NUMH, LISTHPTR, LISTH, H, S, 
     .                 E1, E2, SIGMA, NHIST,
     .                 GAMMA, XIJ, INDXUO, NK, KPOINT, WK, EO, NUOTOT )
C **********************************************************************
C Subroutine to calculate the proyected density of states on the
C atomic orbitals for a given eigenvalue spectra
C Written by J. Junquera and E. Artacho, November 1999.
C ***********  INPUT  **************************************************
C INTEGER NO                  : Number of basis orbitals in the supercell
C INTEGER NSPIN               : Spin polarizations (1 or 2)
C INTEGER MAXSPN              : Second dimension of eo and qo
C                               (Max number of different spin polarizations)
C INTEGER MAXUO               : Maximum number of atomic orbitals in the unit 
C                               cell. First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
C INTEGER MAXNH               : Maximum number of orbitals interacting 
C                               with any orbital
C INTEGER MAXO                : First dimension of eo 
C INTEGER NUMH(NUO)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C INTEGER LISTH(MAXNH)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTHPTR(NUO)       : Pointer to each row (-1) of the
C                               density matrix
C REAL*8  H(MAXNH,NSPIN)      : Hamiltonian in sparse format
C REAL*8  S(MAXNH)            : Overlap in sparse format
C REAL*8  E1, E2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C REAL*8  SIGMA               : Width of the gaussian to expand the eigenvalues
C INTEGER NHIST               : Number of subdivisions of the histogram
C LOGICAL GAMMA               : Only gamma point?
C REAL*8  XIJ(3,MAXNH)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C INTEGER INDXUO(NO)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the nuber of orbitals in the unit cell
C INTEGER NK                  : Number of k points
C REAL*8  KPOINT(3,NK)        : k point vectors
C REAL*8  WK(NK)              : k point weights (must sum one)
C REAL*8  EO(MAXUO,MAXSPN,NK) : Eigenvalues
C INTEGER NUOTOT              : Total number of orbitals in unit cell
C **********************************************************************

      use precision,   only : dp
      use parallel,    only : Node, Nodes
      use fdf
      use atomlist,    only : xa, isa, iphorb, iaorb
      use atmfuncs,    only : zetafio, mofio, lofio, cnfigfio, labelfis
      use flib_wxml
      use xml,         only : xml_dump_attribute, xml_dump_element
      use ionew,       only : IOnode
      use densematrix, only : Haux, Saux, psi
      use alloc
#ifdef MPI
      use parallelsubs, only : GetNodeOrbs
#endif

      implicit none

      INTEGER
     .  NO, NSPIN, MAXSPN, MAXUO, MAXNH, NK, NHIST, MAXO, NUOTOT

      INTEGER 
     .  NUMH(*), LISTH(MAXNH), LISTHPTR(*), INDXUO(NO)

      real(dp)
     .  H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA,  
     .  XIJ(3,MAXNH), KPOINT(3,NK), WK(NK), EO(MAXO,MAXSPN,NK)

      LOGICAL
     .  GAMMA

C Dynamic arrays -------------------------------------------------------
      real(dp), DIMENSION(:,:)  , ALLOCATABLE :: DTOT
      real(dp), DIMENSION(:,:,:), ALLOCATABLE :: DPR

C Internal variables ---------------------------------------------------
      INTEGER
     .  NUO, NHS, NPSI, IUO, IHIST, ISPIN, 
     .  IUNIT1, IUNIT2, I

      integer iat, spec, ii, iorb

      real(dp), dimension(:), allocatable  :: tmp

      character*40 pos_string

      CHARACTER
     .  SNAME*30, FNAMETOT*37, FNAMEPRO*37, PASTE*53,  POINT

      real(dp)
     .  DELTA, ENER, EV

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE, 
     .  PDOSK, TIMER

      type(xmlf_t) :: xf            ! For new XML output

      CALL TIMER( 'pdos', 1)

C Find the intervals between the subdivisions in the energy scale ------
      DELTA = (E2 - E1) / (NHIST-1)
      EV = 13.6058D0
      if (Node.eq.0) then
        SNAME = FDF_STRING('SystemLabel','siesta')
      endif
      POINT = '.'

C Find number of orbitals per unit cell 
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

C Check internal dimensions --------------------------------------------
      IF ( NSPIN.LE.2 .AND. GAMMA) THEN
        NHS  = NUOTOT * NUO
        NPSI = NUOTOT * MAXUO * NSPIN
      ELSE IF ( NSPIN.LE.2 .AND. .NOT.GAMMA) THEN
        NHS  = 2 * NUOTOT * NUO
        NPSI = 2 * NUOTOT * NUO
      ENDIF

C Allocate local arrays ------------------------------------------------
      call re_alloc(Haux,1,nhs,name='Haux',routine='pdos')
      call re_alloc(Saux,1,nhs,name='Saux',routine='pdos')
      call re_alloc(psi,1,npsi,name='psi',routine='pdos')
      ALLOCATE(DTOT(NHIST,2))
      CALL MEMORY('A','D',2*NHIST,'pdos')
      ALLOCATE(DPR(NHIST,NUOTOT,2))
      CALL MEMORY('A','D',2*NHIST*NUOTOT,'pdos')

C Initialize the projected density of states ---------------------------
      DO ISPIN = 1, 2
        DO IHIST = 1, NHIST
          DTOT(IHIST,ISPIN) = 0.D0
          DO IUO = 1, NUOTOT
            DPR(IHIST,IUO,ISPIN) = 0.D0
          ENDDO
        ENDDO
      ENDDO

C Call appropiate routine ----------------------------------------------
      IF (NSPIN.LE.2 .AND. GAMMA) THEN
        CALL PDOSG( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH,
     .              MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .              E1, E2, NHIST, SIGMA, INDXUO, EO,
     .              HAUX, SAUX, PSI, DTOT, DPR, NUOTOT )
      ELSE IF ( NSPIN.LE.2 .AND. .NOT.GAMMA) THEN
        CALL PDOSK( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH,
     .              MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .              E1, E2, NHIST, SIGMA, 
     .              XIJ, INDXUO, NK, KPOINT, WK, EO,
     .              HAUX, SAUX, PSI, DTOT, DPR, NUOTOT )
      ENDIF

      if (IOnode) then
C Open file for write on I/O node
        FNAMETOT = PASTE(SNAME,'.DOS')
        CALL IO_ASSIGN(IUNIT1)
        OPEN(UNIT=IUNIT1, FILE=FNAMETOT, FORM='formatted', 
     .       STATUS='unknown') 

C Output histogram
        DO IHIST = 1, NHIST
          ENER = E1 + (IHIST-1) * DELTA
          WRITE(IUNIT1,'(3f20.5)') ENER*EV,DTOT(IHIST,1)/EV,
     .         DTOT(IHIST,2)/EV
        ENDDO

C Close file for write
        CALL IO_CLOSE(IUNIT1)
      endif

C New writing
      if (IOnode) then
         call xml_OpenFile("pdos.xml",xf)
         call xml_AddXMLDeclaration(xf,"UTF-8")

        FNAMEPRO = PASTE(SNAME,'.PDOS')
        CALL IO_ASSIGN(IUNIT2)
        open(iunit2,file=FNAMEPRO,form='formatted',status='unknown')
        call xml_NewElement(xf,"pdos")
        call xml_NewElement(xf,"nspin")
        call xml_AddPcdata(xf,str(nspin), line_feed=.false.)
        call xml_EndElement(xf,"nspin")

        write(iunit2,'(a)') '<pdos>'
        write(iunit2,'(a,i1,a)') '<nspin>', nspin, '</nspin>'
        write(iunit2,'(a,i4,a)') '<norbitals>', nuotot, '</norbitals>'
        call xml_NewElement(xf,"energy_values")
        call xml_AddAttribute(xf,"units","eV")
        write(iunit2,'(a)') '<energy_values units="eV">'

        allocate(tmp(2*nhist))
        do ihist = 1,nhist
           ENER = E1 + (IHIST-1) * DELTA
           tmp(ihist) = ener*eV
           write(iunit2,'(f20.5)') ener*eV
        enddo
        write(iunit2,'(a)') '</energy_values>'

        call xml_AddArray(xf,tmp(1:nhist))
        call xml_EndElement(xf,"energy_values")

        do i = 1, nuotot
          iat = iaorb(i)
          iorb = iphorb(i)
          spec = isa(iat)

          call xml_NewElement(xf,"orbital")
          write(iunit2,'(a)') '<orbital '
          call xml_AddAttribute(xf,"index",str(i))
          call xml_dump_attribute(iunit2,"index",str(i))
          call xml_AddAttribute(xf,"atom_index",str(iat))
          call xml_dump_attribute(iunit2,"atom_index",str(iat))
          call xml_AddAttribute(xf,"species",trim(labelfis(spec)))
          call xml_dump_attribute(iunit2,"species",
     .                           trim(labelfis(spec)))
          write(pos_string,'(3f11.6)') (xa(ii,iat),ii=1,3)
          call xml_AddAttribute(xf,"position",pos_string)
          call xml_dump_attribute(iunit2,"position",pos_string)
          call xml_AddAttribute(xf,"n",str(cnfigfio(spec,iorb)))
          call xml_dump_attribute(iunit2,"n",str(cnfigfio(spec,iorb)))
          call xml_AddAttribute(xf,"l",str(lofio(spec,iorb)))
          call xml_dump_attribute(iunit2,"l",str(lofio(spec,iorb)))
          call xml_AddAttribute(xf,"m",str(mofio(spec,iorb)))
          call xml_dump_attribute(iunit2,"m",str(mofio(spec,iorb)))
          call xml_AddAttribute(xf,"z",str(zetafio(spec,iorb)))
          call xml_dump_attribute(iunit2,"z",str(zetafio(spec,iorb)))
          write(iunit2,'(a)') '> '

!------------------------------------------------------------
          call xml_NewElement(xf,"data")
          if (nspin.eq.1) then
             call xml_AddArray(xf,dpr(1:nhist,i,1)/eV)
          elseif (nspin .eq. 2) then
             tmp(1:2*nhist-1:2) = dpr(1:nhist,i,1)/eV
             tmp(2:2*nhist:2) = dpr(1:nhist,i,2)/eV
             call xml_AddArray(xf,tmp(1:2*nhist),format="(2f20.6)")
          endif
          call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"orbital")
!------------------------------------------------------------

          write(iunit2,'(a)') '<data>'
          do ihist = 1,nhist
            if (nspin.eq.1) then
              write(iunit2,'(f10.5)') dpr(ihist,i,1)/eV
            elseif (nspin .eq. 2) then
              write(iunit2,'(2f10.5)') dpr(ihist,i,1)/eV,
     .                                 dpr(ihist,i,2)/eV
            endif
          enddo
          write(iunit2,'(a)') '</data>'
          write(iunit2,'(a)') '</orbital>'
        enddo

C Close file
        call xml_EndElement(xf,"pdos")
        call xml_Close(xf)
        write(iunit2,'(a)') '</pdos>'
        call io_close(iunit2)

C Free local workspace array
        deallocate(tmp)

      endif
            

C Free local arrays ----------------------------------------------------
      CALL MEMORY('D','D',SIZE(DPR),'pdos')
      DEALLOCATE(DPR)
      CALL MEMORY('D','D',SIZE(DTOT),'pdos')
      DEALLOCATE(DTOT)

      CALL TIMER( 'pdos', 2)

      RETURN
      END

