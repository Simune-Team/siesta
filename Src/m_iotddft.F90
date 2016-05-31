      MODULE m_iotddft

!     This module has different subroutines to write the total 
!     E_KS, the instantaneous so-called Eigen values and dipole moment
!     at every time step and time-dependent density (rho) after every 
!     given number of steps in case of TDDFT calculations depending 
!     on the user's choice.
!     Based on the modified version of Daniel Sanchez Portal's 
!     original  subroutines. 
!     Rafi Ullah November 2014.
!

      USE m_dipol,          ONLY: dipol
      USE m_steps,          ONLY: fincoor
      USE files,            ONLY: slabel, label_length
      USE siesta_options,   ONLY: eigen_time, dip_time, etot_time, tdsaverho, &
                                  tdsaverho, ntdsaverho
      USE wavefunctions,    ONLY: wavef_ms
      USE parallel,         ONLY: IONode
      USE files,            ONLY: filesOut_t
      
      IMPLICIT NONE
       
      CHARACTER(LEN=label_length+7), EXTERNAL  :: paste, npaste
      CHARACTER(LEN=15)  :: fform


      CONTAINS

      SUBROUTINE write_tddft(totime,istp,itd, ntd,rstart_time,         &
                             etot,eigen, maxo, nspin, nk)

      INTEGER, INTENT(IN)           :: istp, itd, ntd, maxo, nk, nspin
      DOUBLE PRECISION, INTENT(IN)  :: totime, rstart_time, etot
      DOUBLE PRECISION, INTENT(IN)  :: eigen(maxo,nspin,nk) 
      TYPE (filesOut_t)             :: filesOut 
      LOGICAL, SAVE                 :: laststp = .false.
      
      IF (istp .gt. fincoor .AND. itd .gt. ntd) THEN
         laststp = .true.
      END IF
      
      IF (dip_time) THEN
      CALL iodipole (totime, dipol, laststp, rstart_time)
      END IF 
      
      IF (etot_time) THEN
      CALL ioetot   (totime, etot, laststp, rstart_time)
      END IF
      
      IF (eigen_time) THEN
      CALL ioeigenvalues (totime, eigen, laststp, rstart_time, &
                                 maxo, nspin, nk)
      END IF

      IF (tdsaverho) THEN
        IF (mod(istp,ntdsaverho) .eq. 0) THEN
          filesOut%tdrho = npaste (istp, '.TDRho')
        ELSE
          filesOut%tdrho = ' '
        END IF
      END IF

      END SUBROUTINE write_tddft
!----------------------------------------------------------------
       
      SUBROUTINE  iodipole (totime, dipole,lastistp,rstart_time)
       
       
       CHARACTER(LEN=label_length+3) :: dipolefile
       DOUBLE PRECISION   :: dipole(3), extfield(3), totime, rstart_time
       INTEGER, SAVE      :: iu
!       CHARACTER(LEN=15)  :: fform
       LOGICAL,    SAVE      :: frstme  = .true.
       LOGICAL, INTENT(IN)  :: lastistp

!      Only first node writes
       IF(IONode) THEN

       IF (frstme) THEN
         dipolefile = paste (slabel, '.dipol_vs_time')  
         call io_assign( iu )
         fform='formatted'
         OPEN( iu, FILE=dipolefile, FORM=fform, STATUS='unknown' )
!        write(iu,'(a,3f15.6)') '#',extfield(1), extfield(2), extfield(3)
         frstme = .false.
       END IF
       IF (totime .gt. rstart_time) THEN
          WRITE(iu,'(4f15.6)')                                         &
          totime,                                                      &
          dipole(1),dipole(2), dipole(3)
       END IF
       IF (lastistp) call io_close(iu)
        
       END IF ! IONode
      END SUBROUTINE iodipole
!----------------------------------------------------------------------
      SUBROUTINE ioetot (totime, etot, lastistp, rstart_time)
        
       DOUBLE PRECISION         ::totime, etot, eV, rstart_time
       INTEGER                  :: iu, istp, itd, ntd
       LOGICAL                  :: lastistp
       LOGICAL, SAVE            :: frstme = .true. 
       SAVE                     :: iu, eV
       CHARACTER(LEN=70)        :: etotfile
        
        IF(IONode) THEN

        IF (frstme) THEN
          etotfile = paste (slabel, '.etot_vs_time')
          CALL io_assign ( iu )
          fform = 'formatted'
          OPEN (iu, FILE=etotfile, FORM=fform, POSITION='APPEND',      &
                STATUS='UNKNOWN')
          frstme = .false.
          eV     = 1.d0/13.60580d0
        END IF
        IF (totime .gt. rstart_time) THEN ! To avoid rewriting the already written data in case of restart.
           WRITE (iu, '(2f15.6)') totime, etot/eV
        END IF
        IF (lastistp) CALL io_close(iu)
        END IF ! IONode
      END SUBROUTINE ioetot
!------------------------------------------------------------------------

      SUBROUTINE ioeigenvalues (totime, eigen, lastistp, rstart_time, &
                                 maxo, nspin, nk)
    
       INTEGER            :: maxo, nspin, nk, ik, ispin, ie, iu
       INTEGER            :: nocc(nk,nspin)
       DOUBLE PRECISION   :: eV, totime, rstart_time
       DOUBLE PRECISION   :: eigen(maxo,nspin,nk)
       LOGICAL            :: lastistp
       LOGICAL, SAVE      :: frstme = .true.
       CHARACTER(LEN=70)  :: eigenfile

        IF (IONode) THEN 

        IF (frstme) THEN
          eV =1.0d0/13.60580d0
          frstme = .false.
          eigenfile = paste (slabel, '.eigen_vs_time')
          CALL io_assign (iu) 
          fform = 'formatted'
          OPEN(iu, FILE = eigenfile, FORM=fform, STATUS='unknown')
          WRITE(iu,*) '#  ', nspin, nk
          DO ispin=1,nspin          
             WRITE(iu,*) '#  ',((wavef_ms(ik,ispin)%dim2), ik=1,nk)
          END DO
        END IF

        IF (totime .gt. rstart_time) THEN
          WRITE(iu,'(11f12.5,/,(5x,10f12.5))') totime,                 &
                (((eigen(ie,ispin,ik)/eV, ie=1,(wavef_ms(ik,ispin)%dim2)),        &
                ispin=1,nspin), ik=1,nk)
        END IF
        
        IF (lastistp) CALL io_close(iu)
        END IF ! IONode
       END SUBROUTINE ioeigenvalues 




       END MODULE m_iotddft
