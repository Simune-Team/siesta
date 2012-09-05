      MODULE m_ts_in_siesta
!
! Routines that are called from SIESTA routines
!
!=============================================================================
! CONTAINS:
!          1) transiesta
!          2) TSVHFix

  
      implicit none

      public :: transiesta, TSVHFix

      private

      CONTAINS


! ##################################################################
! ##                                                              ##       
! ##                       "TRANSIESTA"                           ##
! ##                                                              ##       
! ##          Non-equilibrium Density Matrix Subroutine           ##
! ##                   to be called from SIESTA                   ##
! ##                                                              ## 
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##               Kurt Stokbro, ks@mic.dtu.dk                    ## 
! ##               Mikroelektronik Centret (MIC)                  ##
! ##           Technical University of Denmark (DTU)              ##
! ##                                                              ##
! ##   Modified by Nick Papior Andersen (electrode repetitions)   ##
! ##   Modified by Nick Papior Andersen, major clean up           ##
! ##                                                              ##
! ##################################################################
!
!
! Tight-binding density-matrix/transport program for the SIESTA
! package.
! Copyright by Mads Brandbyge, 1999, 2000, 2001, 2002.
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the authors.
!
!=================================================================
! Modified by F. D. Novaes to allow calculations with k-point 
! samplings other than gamma point.
!=================================================================
!
!
! Assumes listh and listd are the same
! Assumes that indxuo is consecutive in numbers:
!     do i = 1 , no_u
!       indxuo(i) = i
!     end do
!   both in Gamma and (.not. Gamma) point.
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! UNITS:
!     Ry, and Bohr is used throughout the program unless explicitly 
!     stated.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine transiesta(ucell,gamma,lasto,iaorb,
     &     no_l, no_s, no_u, nspin, xa, na_u, maxnh, numh, listhptr,
     &     listh,  xij, indxuo, Hs, Ss, Dnew, Enew, Ef,
     &     TSiscf,Qtot)

      use parallel  , only: Node, Nodes, IONode
#ifdef MPI
      use mpi_siesta, only: MPI_Comm_World,MPI_SUM
      use mpi_siesta, only: MPI_Bcast
      use mpi_siesta, only: DAT_dcomplex, DAT_double,MPI_Logical
      use mpi_siesta, only: MPI_Integer

      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use parallelsubs, only : LocalToGlobalOrb
#endif
      use m_ts_options, only : ReUseGF
      use m_ts_options, only : HSFileL,HSFileR
      use m_ts_options, only : NUsedAtomsL,NUsedAtomsR
      use m_ts_options, only : NRepA1L,NRepA2L,NRepA1R,NRepA2R
      use m_ts_options, only : NBufAtL, NBufAtR
      use m_ts_options, only : GFFileL, GFFileR
      use m_ts_options, only : IsVolt, VoltFDF, VoltL, VoltR
      use m_ts_options, only : UseBulk, TriDiag, updatedmcr
      use m_ts_options, only : nvolt, ElecValenceBandBot, GFTitle
      use m_ts_options, only : ChargeCorr, ChargeCorr_factor
      use m_ts_options, only : CCEmin, GFEta, kT, smethod
      use m_ts_options, only : nline,ncircle,npol,nvolt

      use m_ts_contour, only : setupcontour, distrcontour
      use m_ts_GF,      only : read_Green
      use m_ts_scattering, only : getsfe, getGF_Tri
      use m_ts_scattering, only : GFGammaGF, weightDM, weightDMC
      use m_ts_scattering, only : get_scat_region

      use m_ts_kpoints, only : ts_nkpnt,ts_kpoint,ts_kweight
      use m_ts_kpoints, only : ts_gamma_scf

#ifdef MPI
      use m_glob_sparse, only: glob_sparse_arrays, glob_sparse_matrix
#endif
      use m_hs_matrix, only : set_HS_matrix, matrix_rem_left_right
      use m_hs_matrix, only : matrix_symmetrize

      use units, only : eV,pi
      use sys, only : die
      use precision, only : dp


!======================================================================
!   FROM SIESTA:
!======================================================================
      real(dp) ucell(3,3)       ! Unit cell parameters
      logical gamma             ! True if gamma point calculation
      integer na_u               ! No. atoms in unit cell
      
      integer lasto(0:na_u)      ! last-orbital index on atoms
      integer nspin             ! Number of spin components (1 or 2)
!                               ! Can only be 1 for the time being!

      integer no_s              ! no. orbs. in supercell
      integer no_l              ! no. orbs. in unit cell (local)
      integer no_u              ! no. orbs. in unit cell (global)
      integer maxnh             ! Maximum number of nonzero elements of 
!                               ! each row of hamiltonian matrix
      integer  iaorb(no_u)      ! The atomic index of the orbitals
      real(dp) xa(3,na_u)        ! Atomic coordinates in unit cell
      integer numh(no_l)        ! Number of nonzero elements of each row
      integer listh(maxnh)      ! Nonzero hamiltonian-matrix element
!                               ! column indexes for each matrix row
!                               ! of hamiltonian matrix
      integer listhptr(no_l)    ! Pointer to each row (-1) of the
!                               ! hamiltonian matrix
!                               ! column indexes for each matrix row
      real(dp)  Hs(maxnh,nspin) !Hamiltonian in sparse form
      real(dp)  Ss(maxnh)       ! Overlap in sparse form
      real(dp)  xij(3,maxnh)    ! Vectors between orbital centers (sparse)
      integer indxuo(no_s)      ! Index of equivalent orbital in super cell
!                               ! Unit cell orbitals must be the first in
!                               ! orbital lists, i.e. indxuo(i).le.no_l,
!                               ! for all i with no_l the number of 
!                               ! orbitals in unit cell
      
      integer TSiscf ! Current iteration

! The fermi-energy of the CONTACT (without electrodes) is initially determined
! by SIESTA. This is taken as the starting energy shift of the contact region:
! The electrodes in equilibrium has per definition Fermi energy at ZERO.
      real(dp) Ef                 ! Fermi energy (used in CONTACT) 
      real(dp) Qtot
!======================================================================
!   BACK TO SIESTA:
!======================================================================

      real(dp) Dnew(maxnh,nspin) ! Output Density Matrix
      real(dp) Enew(maxnh,nspin) ! Output Energy-Density Matrix

!======================================================================
!     Greens functions:
!======================================================================

! Full Green's function, H and S of contact:
      complex(dp), dimension(:,:), pointer :: H, S
      complex(dp), dimension(:),   pointer :: GF

!     Energy Contour:
      complex(dp), dimension (:), pointer, save :: contour,wgf
!     NCONTOUR0 is the number of contour point if it is an equilibrium
!     situation
!     NCONTOUR is the total number of contour points
! No. contour points (NCONTOUR if VOLT=0 else NCONTOUR + NVOLT)
      integer, save :: NCONTOUR0,NCONTOUR
! The number of energy points on the contour
      integer, save :: NEn

      integer, save :: NPARACONTOUR   ! No. points on parallel contour

!     Energy points and their weight distributed on nodes
      complex(dp), dimension(:,:), allocatable, save :: paracontour
      complex(dp), dimension(:,:), allocatable, save :: paraWGF

!     (point,node) --> part of contour :
!              ('L','R','V','N') Left,Right,Voltage, No voltage
      character, dimension(:,:), allocatable, save ::  contourpart
!     When looping a single char is needed
      character                                    :: Part

!     No. actual energy points to read-in for given parallel energy point

      integer, dimension(:), allocatable, save ::  points2read

!==================================================================
!     Left/Right electrode variables:
!==================================================================

!     All variables are denoted by L/R for Left/Right, respectively.

!     File units for Left/Right GF
      integer :: uGFL,uGFR

!     Here we write:
!     Hamiltonian, overlap, Green's function, scattering matrix,
!     k-points, q-points, k/q-weights, fermi level shift due to V,
!     etc.

!     Fermi levels
!       Left   = 0.5 * Volt
!       Device = 0.0
!       Right = -0.5 * Volt
      real(dp), save :: EFermiL,EFermiR,EFermi0

!     Number of atoms in the electrode (before expansion)
      integer, save :: nuaL_GF, nuaR_GF
!     Number of atoms in the electrode (after expansion)
      integer, save :: nuaL, nuaR
!     Number of orbitals in the electrodes (before expansion)
      integer, save :: noL_GF, noR_GF
!     Number of orbitals in the electrodes (after expansion)
      integer, save :: noL,noR

!     The lasto shortened to only the electrodes
      integer, dimension(:), allocatable, save :: lastoL,lastoR

!     Read-in bulk H, S, for a given q-points
      complex(dp), dimension(:,:,:), allocatable :: HAAL, SAAL
      complex(dp), dimension(:,:,:), allocatable :: HAAR, SAAR
!     Self energies and scattering matrices
!     Gamma=i*0.5*(Sigma - Sigma^dagger) 
      complex(dp), dimension(:,:),   allocatable :: SFEL, GAMMAL
      complex(dp), dimension(:,:),   allocatable :: SFER, GAMMAR

!     Electrode k-points kpar-points and their weights:
!     Notice that these are at the moment completely identical
!     to the TranSIESTA k-points at time of creation.
!     Therefore they should be the same as ts_kpoints until 
!     implementation is completed to handle different electrode k-point
!     sampling
      integer                                :: nkparL,nkparR
      real(dp), dimension (:,:), allocatable :: kparL ,kparR
      real(dp), dimension (:)  , allocatable :: wkparL,wkparR

!     q-points and their weights:
!     The notation of b tells the programmer that
!     q is in units of b_i**-1 and NOT in Bohr**-1
      integer , save                         :: nqL,nqR
      real(dp), dimension (:,:), allocatable :: qLb,qRb
      real(dp), dimension (:),   allocatable :: wqL,wqR

!==================================================================

!==================================================================
!     Density matrix variables:
!==================================================================

!     To save space we will save the temporary DML,EDML in DM,EDM. 
!     contour part of DM integration 
!     real axis part of DM integration
!     L,R estimates of Edm

      real(dp), dimension(:), allocatable :: DM,DMR
      real(dp), dimension(:), allocatable :: DMneqL,DMneqR
      real(dp), dimension(:), allocatable :: EDM, EDMR

      complex(dp), dimension(:), allocatable :: DMCplx,DMRCplx
      complex(dp), dimension(:), allocatable :: DMneqLCplx,DMneqRCplx
      complex(dp), dimension(:), allocatable :: EDMCplx,EDMRCplx
      complex(dp), dimension(:), allocatable :: work1Cplx
!==================================================================

!==================================================================
!     Size of problem:
!==================================================================
      
!     Number of orbitals in the buffer regions
      integer, save :: noBufL,noBufR

      integer, save :: nou ! no. of orbs. from SIESTA minus buffer-orbitals
!                          ! i.e. no. used orbitals


!==================================================================
!   Bookkeeping of charges in different regions:
!==================================================================

      real(dp), save :: QC0  ! Init(SIESTA) no. electrons in update region
      real(dp), save :: QCn0 ! Init(SIESTA) no. electrons outside update region
      real(dp)       :: QC   ! no. electrons in update region

      ! Collecting charges in updated region on k-point resolution
      real(dp), dimension(:,:,:), allocatable :: Q_kpt
#ifdef MPI
      real(dp), dimension(:,:,:), allocatable :: Q_kpt_loc
      real(dp)                                :: Qcloc ! Used for reduction
#endif
      ! Integer array holding information on the region content.
      ! It is built together with listud[g] and has the following
      ! values:
      !   0  :  Unknown (actually, all other regions which is NOT unknown)
      !   1  :  Left buffer region
      !   2  :  Left buffer-electrode coupling
      !   3  :  Left electrode
      !   4  :  Left-device coupling
      !   5  :  Device
      !   6  :  Right-device coupling
      !   7  :  Right electrode
      !   8  :  Right buffer-electrode coupling
      !   9  :  Right buffer region
      integer, dimension(:), allocatable, save :: region
#ifdef MPI
      integer, dimension(:), allocatable, save :: regiong
#endif
      ! Variable for doing a renormalization of the charge
      real(dp), dimension(:), allocatable :: reD
!==================================================================

!==================================================================
! Helpers, workspace, tempos, auxilaries etc...
!==================================================================
      ! Temporary k-points
      real(dp) :: k(3) ! Temp. k-point
      real(dp) :: kxij ! exp(k*xij)
      complex(dp) :: cphase, ccphase ! phases
      
      ! Work arrays
      real(dp),    dimension(:), allocatable :: work1
      complex(dp), dimension(:), allocatable :: work2
      complex(dp), dimension(:), allocatable :: dosdum

! Energies on the contour
      complex(dp) :: Z,W,ZW
!==================================================================

! Different logicals
      logical :: ioinC,joinC,updatedm

! Loop counters:
      integer :: iPE,ispin
      integer :: ts_ikpt
      integer :: i,j,ic,jc,ind
      integer :: ia,ia2
      integer :: io,jo,iuo,juo,iio,is
! error logicals or integers
      logical :: errorGS
      integer :: ierr,itmp

! array for determining whether an object should be updated
      logical, dimension (:), allocatable, save :: listud
#ifdef MPI
      logical, dimension (:), allocatable, save :: listudg
#endif

#ifdef MPI
      ! Used for obtaining location of information on target node
      integer :: BNode
! All these variables are the globalised versions of their equivalents
! without the notation of "g"
      integer                               :: maxnhg
      integer, dimension(:),    allocatable :: numhg
      integer, dimension(:),    allocatable :: listhptrg, listhg
      real(dp), dimension(:,:), allocatable :: xijg
      real(dp), dimension(:,:), allocatable :: Dnewg, Enewg
      real(dp), dimension(:),   allocatable :: Hsg, Ssg
!     The obvious error for MPI
      integer :: MPIerror
#endif

      logical, save :: firsttime = .true.

!==================================================================
!     External routines and interfaces:
!==================================================================
      external io_assign, io_close
      
      interface
         function paste(s1,s2)
         character(LEN=*), intent(in) :: s1,s2
         character(LEN=33) :: paste
         end function paste
      end interface

!  I need to create a general case of electrodes with potentially 
!  different numbers of k points, such as an electrode and a bulk nanotube.

#ifdef DEBUG
      call write_debug( 'PRE transieasta' )
#endif

#ifdef MPI
      call MPI_Bcast(Ef,1,DAT_double,0, MPI_Comm_World,MPIerror)
#endif

! ============================================================
!       ----------- BEGIN: first time only ------------
      if (firsttime) then
! ============================================================
!
!     We take the common Fermi energy of the electrodes to be zero!
!     
         EFermiL = VoltL !Ry
         EFermiR = VoltR !Ry
         EFermi0 = 0.0_dp ! (EfermiL + EfermiR)*.5_dp !Ry

!
!     Setup contour
!
         call setupcontour(NEn,EFermi0,EFermiL,EFermiR,
     &        WGF,contour,NCONTOUR,NCONTOUR0,
     &        nline,ncircle,npol,nvolt,
     &        CCEmin,GFEta,kT,smethod)

!
!     Distribute contour points on nodes
!
         allocate( paracontour(1:NCONTOUR,0:Nodes-1)  )
         allocate( paraWGF(1:NCONTOUR,0:Nodes-1)      )
         allocate( contourpart(1:NCONTOUR,0:Nodes-1)  ) ! is CHAR
         allocate( points2read(1:NCONTOUR)            )
         call memory('A','Z',2*NCONTOUR*Nodes,'transiesta')
         call memory('A','I',NCONTOUR,'transiesta')


! (note we're allocating more than we need here if there're more than 1 node)

         call distrcontour(
     &        WGF,contour,NCONTOUR,NCONTOUR0,NVOLT,IsVolt,Nodes,
     &        NPARACONTOUR,paracontour,paraWGF,
     &        contourpart,points2read)


!     the first NBufAtL atoms will be removed
         noBufL = 0
         do ia = 1 , NBufAtL
            noBufL = noBufL+(lasto(ia)-lasto(ia-1))
         end do                 !ia
!     the last NBufAtR atoms will be removed
         noBufR = 0
         do ia = na_u - NBufAtR+1 , na_u
            noBufR = noBufR+(lasto(ia)-lasto(ia-1))
         end do                 !ia

!
! No. states minus buffers
!     
#ifdef MPI
         nou = no_u - (noBufL+noBufR)
#else
         nou = no_l - (noBufL+noBufR)
#endif

! -------------------------------------------------------
!
!     GF generation:
!
         allocate(dosdum(NEn))
         call memory('A','Z',NEn,'transiesta')

         call ts_do_Green('L',HSFileL, GFFileL, GFTitle, 
     &        ElecValenceBandBot, ReUseGF,
     &        ts_nkpnt,ts_kpoint,ts_kweight,
     &        NBufAtL,NUsedAtomsL,NRepA1L,NRepA2L,
     &        ucell,xa,na_u,NEn,contour,wGF,EFermiL,dosdum,nspin,
     &        nuaL_GF,noL_GF)

!     Create the Right GF file
         call ts_do_Green('R',HSFileR,GFFileR, GFTitle,
     &        ElecValenceBandBot, ReUseGF,
     &        ts_nkpnt,ts_kpoint,ts_kweight,
     &        NBufAtR,NUsedAtomsR,NRepA1R,NRepA2R,
     &        ucell,xa,na_u,NEn,contour,wGF,EFermiR,dosdum,nspin,
     &        nuaR_GF,noR_GF)

         call memory('D','Z',NEn,'transiesta')
         deallocate(dosdum)

! Make the expansion of the electrodes in atoms and orbitals
         nuaL = nuaL_GF*NRepA1L*NRepA2L
         nuaR = nuaR_GF*NRepA1R*NRepA2R
         noL  = noL_GF *NRepA1L*NRepA2L
         noR  = noR_GF *NRepA1R*NRepA2R

         ! Do a crude check of the sizes
         if ( nou < noL .or. nou < noR ) then
            call die("The contact region size is &
     &smaller than the electrode size. Please correct.")
         end if

! -------------------------------------------------------
         
!     Print out information in Green's function files
!     Show the number of used atoms and orbitals
         if ( IONode ) then
            write(*,'(/,a,i6,'' / '',i6)')
     &           'Left : GF atoms    / Expanded atoms    : ',nuaL_GF
     &           ,nuaL
            write(*,'(a,i6,'' / '',i6)')
     &           'Left : GF orbitals / Expanded orbitals : ',noL_GF,noL
            write(*,'(a,i6,'' / '',i6)')
     &           'Right: GF atoms    / Expanded atoms    : ',nuaR_GF
     &           ,nuaR
            write(*,'(a,i6,'' / '',i6)')
     &           'Right: GF orbitals / Expanded orbitals : ',noR_GF,noR
         end if
         
! expected no. states on Electrode atoms within the Green's function file: 
! Left
         allocate(lastoL(0:nuaL_GF)) 
         call memory('A','I',nuaL_GF+1,'transiesta')

         lastoL(0)=0
         ia2=0
         do ia=NBufAtL+1,NBufAtL + nuaL, NRepA1L*NRepA2L
            ia2=ia2+1
            lastoL(ia2)=lastoL(ia2-1) + (lasto(ia) - lasto(ia-1))
         end do                 !ia

         if(lastoL(nuaL_GF) .ne. noL_GF) then
            if(IONode)
     &           write(*,*) 'ERROR: lastoL,noL_GF',lastoL,noL_GF
            call die('ERROR: Unexpected no. orbs. in L elec.')
         end if

! Right
         allocate(lastoR(0:nuaR_GF))    
         call memory('A','I',nuaR_GF+1,'transiesta')

         lastoR(0)=0
         ia2=0
         do ia=na_u-(nuaR + NBufAtR)+1,na_u-NBufAtR, NRepA1R*NRepA2R
            ia2=ia2+1
            lastoR(ia2)=lastoR(ia2-1) + (lasto(ia) - lasto(ia-1))
         end do                 !ia

         if(lastoR(nuaR_GF) .ne. noR_GF) then
            if(IONode)
     &           write(*,*) 'ERROR: lastoR,noR_GF',lastoR,noR_GF

            call die('ERROR: Unexpected no. orbs. in R elec.')
         end if

      end if                    ! firsttime
! -------------------------------------------------------
      call timer("TS_calc",1)

#ifdef MPI
      call glob_sparse_arrays(no_l,no_u,no_s,maxnh, 
     &     numh ,listh ,listhptr ,xij , Gamma,
     &     maxnhg,
     &     numhg,listhg,listhptrg,xijg)
      
      allocate(Dnewg(maxnhg,nspin))
      allocate(Enewg(maxnhg,nspin))
      call memory('A','D',2*maxnhg*nspin,'transiesta')
      
      do is = 1, nspin
         call glob_sparse_matrix(no_l,no_u,no_s,  
     &        maxnh,  numh , listhptr ,Dnew(:,is) ,Enew(:,is),
     &        maxnhg, numhg, listhptrg,Dnewg(:,is),Enewg(:,is))
      end do
#endif
      call timer("TS_calc",2)
! FIRSTTIME

      if(TSiscf.eq.1) then

! -------------------------------------------------------
!
!     IMPORTANT: Here we set the region which we do not change in DM
!                and initialize.
!   Furthermore we determine the regions of charges.
!   This helps in the discovery of charge redistribution during the
!   self-consistent runs of TranSIESTA.
         if(IONode) then
            write(*,'(/,a)') 'Creating update lists and' //
     &           ' calculating SIESTA charge distribution...'
         endif         


         if(allocated(listud)) then
            call memory('D','L',size(listud),'transiesta')
            deallocate(listud)
         endif
         allocate(listud(maxnh))
         call memory('A','L',maxnh,'transiesta')
         listud = .false.

         if(allocated(region)) then
            call memory('D','I',size(region),'transiesta')
            deallocate(region)
         endif
         allocate(region(maxnh))
         call memory('A','I',maxnh,'transiesta')

         do io = 1,no_l
#ifdef MPI
            call  LocalToGlobalOrb(io, Node, Nodes, iio)
            iuo = indxuo(iio)
#else
            iuo = indxuo(io)
#endif
            ic = iuo - noBufL
            do j = 1,numh(io)
               ind = listhptr(io) + j
               jo = listh(ind)
               juo = indxuo(jo)

               jc = juo - noBufL

               ! Retrieve the scattering position
               region(ind) = get_scat_region(ic,noL,jc,noR,nou)

               ioinC=.false.
               joinC=.false.
               if(UseBulk) then !update only where Ham is updated:
                  if(ic.gt.noL .and. ic.lt.nou-noR+1) ioinC=.true.
                  if(jc.gt.noL .and. jc.lt.nou-noR+1) joinC=.true.

! FDN From fdf Now
                if( .not. updatedmcr ) then                
!
!              update the cross DM terms between R,L and  C
!
                  updatedm = ioinC .or. joinC
                else
!
!              update DM only in C
!
                 updatedm = ioinC .and. joinC
! FDN
                end if
               else
                  if(ic.ge.1 .and. ic.le.nou) ioinC=.true.
                  if(jc.ge.1 .and. jc.le.nou) joinC=.true.
                  updatedm = ioinC .and. joinC
               end if           !UseBulk

!     Always fix bonds between Left/Right electrodes:
               if(jc.le.noL .and. ic.ge.(nou-noR+1))
     &              updatedm=.false.
               if(ic.le.noL .and. jc.ge.(nou-noR+1))
     &              updatedm=.false.
               listud(ind)=updatedm
            end do              !j
         end do                 !io

#ifdef MPI
         if(allocated(listudg)) then
            call memory('D','L',size(listudg),'transiesta')
            deallocate(listudg)
         endif
         allocate(listudg(maxnhg))
         call memory('A','L',maxnhg,'transiesta')
         listudg = .false.

         if(allocated(regiong)) then
            call memory('D','I',size(regiong),'transiesta')
            deallocate(regiong)
         endif
         allocate(regiong(maxnhg))
         call memory('A','I',maxnhg,'transiesta')

         do io = 1,no_u
            iuo = indxuo(io)
            ic = iuo - noBufL

            do j = 1,numhg(io)
               ind = listhptrg(io) + j
               jo = listhg(ind)
               juo = indxuo(jo)
               jc = juo - noBufL
               
               ! Retrieve the scattering position
               regiong(ind) = get_scat_region(ic,noL,jc,noR,nou)

               ioinC=.false.
               joinC=.false.
               if(UseBulk) then !update only where Ham is updated:
                  if(ic.gt.noL .and. ic.lt.nou-noR+1) ioinC=.true.
                  if(jc.gt.noL .and. jc.lt.nou-noR+1) joinC=.true.
! FDN From fdf Now
                  if( .not. updatedmcr ) then                
!
!              update the cross DM terms between R,L and  C
!
                     updatedm = ioinC .or. joinC
                  else
!
!              update DM only in C
!
                     updatedm = ioinC .and. joinC
! FDN
                  end if
               else
                  if(ic.ge.1 .and. ic.le.nou) ioinC=.true.
                  if(jc.ge.1 .and. jc.le.nou) joinC=.true.
                  updatedm = ioinC .and. joinC
               end if           !UseBulk
               if(jc.le.noL .and. ic.ge.(nou-noR+1))
     &              updatedm=.false.
               if(ic.le.noL .and. jc.ge.(nou-noR+1))
     &              updatedm=.false.
               listudg(ind)=updatedm
            end do              !j
         end do                 !io
#endif

! -------------------------------------------------------
!
!     Calculate initial region charges
!
         allocate(Q_kpt(0:9,nspin,1))
         Q_kpt(:,:,:) = 0.0_dp
         call memory('A','D',nspin*10,'transiesta')

         Qcn0 = 0.d0
         do ispin = 1,nspin
            do ind = 1,maxnh
               if(.not.listud(ind)) then
                  Qcn0 = Qcn0 + Dnew(ind,ispin)*Ss(ind)
               endif
               Q_kpt(region(ind),ispin,1) =
     &              Q_kpt(region(ind),ispin,1)+
     &              Dnew(ind,ispin)*Ss(ind)
            end do              !io
         end do                 !ispin

#ifdef MPI
         call MPI_AllReduce(Qcn0,qcloc,1,
     &        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         Qcn0 = Qcloc

         allocate(Q_kpt_loc(0:9,nspin,1))
         call memory('A','D',nspin*10,'transiesta')
         call MPI_AllReduce(Q_kpt,Q_kpt_loc,nspin*10,
     &        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         Q_kpt = Q_kpt_loc
         call memory('D','D',nspin*10,'transiesta')
         deallocate(Q_kpt_loc)
#endif

!     Write out charges in each region
!     If a buffer region exists it will print out information about
!     those as well
         Qc0 = Qtot - Qcn0

         if(IONode) then 
            write(*,'(4(/,a,f12.5))')
     &           'Efermi from SIESTA                   :',ef/eV,
     &           'Total charge                  [Qt0]  :',Qtot,
     &           'Charge in update region       [Qc0]  :',Qc0,
     &           'Charge outside update region  [Qcn0] :',Qcn0
            if ( NBufAtL > 0 ) write(*,'(a,f12.5,/,a,f12.5)')
     &           'Left buffer                          :',Q_kpt(1,1,1),
     &           'Left buffer/left electrode           :',Q_kpt(2,1,1)
            write(*,'(a,f12.5,4(/,a,f12.5))')
     &           'Left electrode                       :',Q_kpt(3,1,1),
     &           'Left electrode/device                :',Q_kpt(4,1,1),
     &           'Device                               :',Q_kpt(5,1,1),
     &           'Device/right electrode               :',Q_kpt(6,1,1),
     &           'Right electrode                      :',Q_kpt(7,1,1)
            if ( NBufAtR > 0 ) write(*,'(a,f12.5,/,a,f12.5)')
     &           'Right electrode/right buffer         :',Q_kpt(8,1,1),
     &           'Right buffer                         :',Q_kpt(9,1,1)
            write(*,'(a,f12.5)')
     &           'Other                                :',Q_kpt(0,1,1)
         end if    
         ! Deallocate Q_kpt as we will reuse it in the update of the region
         call memory('D','D',nspin*10,'transiesta')
         deallocate(Q_kpt)
      end if                    !TSiscf ==1

      call timer("TS_calc",1)


! -------------------------------------------------------
!     Read-in header of Green's functions
!     Prepare for the calculation
!     We read in the k-points that the electrode was generated with.
!     Furthermore we read in the expansion q-points
!     They are communicated in the routine


      if(IONode) then
         call io_assign(uGFL)
         open(file=GFFileL,unit=uGFL,form='unformatted')
         call io_assign(uGFR)
         open(file=GFFileR,unit=uGFR,form='unformatted')
      endif

! Left
      call read_Green(uGFL,EfermiL,ts_nkpnt,NEn,nuaL_GF, 
     &     NRepA1L,NRepA2L,noL_GF,nspin,
     &     nkparL,kparL,wkparL,
     &     nqL,wqL,qLb)

! Right
      call read_Green(uGFR,EfermiR,ts_nkpnt,NEn,nuaR_GF,
     &     NRepA1R,NRepA2R,noR_GF,nspin, 
     &     nkparR,kparR,wkparR,
     &     nqR,wqR,qRb)

!----------------------------------------------
! Allocation
#ifdef MPI
      ! Allocate for the globalized sparse matrices
      allocate(Hsg(maxnhg),Ssg(maxnhg))
      call memory('A','D',maxnhg*2,'transiesta')
#endif

      nullify(H,S)
      allocate(H(nou,nou))
      allocate(S(nou,nou))
      call memory('A','Z',2*nou*nou,'transiesta')
! FDN Complex for k point sampling
      if ( ts_gamma_scf ) then
         allocate(DM(nou*nou))
         allocate(EDM(nou*nou))
         call memory('A','D',2*nou*nou,'transiesta')
      else
         allocate(DMCplx(nou*nou))
         allocate(EDMCplx(nou*nou))
         call memory('A','Z',2*nou*nou,'transiesta')
      end if
      
      allocate(SFEL(noL,noL))
      allocate(SFER(noR,noR))
      allocate(GAMMAL(noL,noL))
      allocate(GAMMAR(noR,noR))
      call memory('A','Z',2*noL**2+2*noR**2,'transiesta')
      allocate(HAAL(noL_GF,noL_GF,nqL))
      allocate(SAAL(noL_GF,noL_GF,nqL))
      allocate(HAAR(noR_GF,noR_GF,nqR))
      allocate(SAAR(noR_GF,noR_GF,nqR))
      call memory('A','Z',2*noL_GF**2*nqL+2*noR_GF**2*nqR,'transiesta')

      if(isvolt)then
! FDN Complex for k point sampling
         allocate(work2(nou*nou))
         call memory('A','D',nou*nou,'transiesta')
         if ( ts_gamma_scf ) then      
            allocate(DMR(nou*nou))
            allocate(DMneqL(nou*nou))
            allocate(DMneqR(nou*nou))
            allocate(EDMR(nou*nou))
            call memory('A','D',4*nou*nou,'transiesta')
         else
            allocate(DMRCplx(nou*nou))
            allocate(DMneqLCplx(nou*nou))
            allocate(DMneqRCplx(nou*nou))
            allocate(EDMRCplx(nou*nou))
            call memory('A','Z',4*nou*nou,'transiesta')
         end if
      endif

! Prepare for the gathering of the charges in the different regions
      allocate(Q_kpt(0:ts_nkpnt,0:9,nspin))
      Q_kpt(:,:,:) = 0.0_dp
      call memory('A','D',(1+ts_nkpnt)*nspin*10,'transiesta')

      nullify(GF)
      allocate(GF(nou*nou))
      call memory('A','Z',nou*nou,'transiesta')

! FDN Spin loop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ispin = 1,nspin        !not really implemented yet!!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
         ! Communicate out the globalized sparse matrices
         ! There is no need to do this within the k-point loop
         ! However, it will require maxnhg*2 more disk space...
         call glob_sparse_matrix(no_l,no_u,no_s,
     &        maxnh,  numh , listhptr ,Hs(:,ispin) ,Ss,
     &        maxnhg, numhg, listhptrg,Hsg,Ssg)
         
#endif

!     Initialize density here instead of in the k-point loop
#ifdef MPI
       do ind = 1,maxnhg
          if(listudg(ind)) then
             Dnewg(ind,ispin) = 0d0
             Enewg(ind,ispin) = 0d0
          end if
       end do
#else
       do ind = 1,maxnh
          if(listud(ind)) then
             Dnew(ind,ispin)  = 0d0
             Enew(ind,ispin)  = 0d0
          end if
       end do
#endif


! we wish to loop over the large k-points... 
!     other sub calls that ts_kpoint is the correct array...
         DO ts_ikpt = 1 , ts_nkpnt
          k(:) = ts_kpoint(:,ts_ikpt)

!     
!     Init DM, EDM
!
          if (ts_gamma_scf) then 
            DM    =  0.d0
            EDM   =  0.d0
            if (isvolt) then
               DMR    = 0.d0
               DMneqL = 0.d0
               DMneqR = 0.d0
               EDMR   = 0.d0
            end if
          else
!     FDN Complex
            DMCplx        = dcmplx(0.d0,0.d0)
            EDMCplx       = dcmplx(0.d0,0.d0)
            if (isvolt) then
               DMRCplx    = dcmplx(0.d0,0.d0)
               DMneqLCplx = dcmplx(0.d0,0.d0)
               DMneqRCplx = dcmplx(0.d0,0.d0)
               EDMRCplx   = dcmplx(0.d0,0.d0)
            end if
          end if
! FDN

! Here we could employ the same scheme as in the TBTrans utility
! This should yield a faster communication as each k-point
! is calculated on separate nodes...
! Then it is only the communication overhead which "slows" the
! computation down.
#ifdef MPI
          call set_HS_matrix(Gamma,ucell,na_u,no_u,no_s,
     &         maxnhg,xijg,numhg,listhptrg,listhg,indxuo,Hsg,Ssg, 
     &         k,H,S,
     &         RemNFirstOrbitals=noBufL,RemNLastOrbitals=noBufR)
#else
          call set_HS_matrix(Gamma,ucell,na_u,no_u,no_s,
     &         maxnh,xij,numh,listhptr,listh,indxuo,Hs(:,ispin),Ss,
     &         k,H,S,
     &         RemNFirstOrbitals=noBufL,RemNLastOrbitals=noBufR)
#endif
          call matrix_rem_left_right(nou,H,S,noL,noR)
          ! The symmetrization has been moved to this routine
          call matrix_symmetrize(nou,H,S,Ef)
          
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
         DO iPE = 1, NPARACONTOUR ! loop over parallel contour points
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=

         Z = paracontour(iPE,Node)    ! ZEnergy point for this node
         Part = contourpart(iPE,Node) ! and place on the contour

!     Weight of energy point:
!     This is where we include the factor-of-two for spin and
!     and the (1/Pi) from DM = Im[G]/Pi
         if(nspin.eq.1) then
            W = (2d0/Pi)*paraWGF(iPE,Node)
         else
            W = (1d0/Pi)*paraWGF(iPE,Node)
         end if

         ! Save Z*W, used for E-arrays
         ZW = Z*W
!----------------------------------------------------------------------
! At this point we need the transformation from
! G(q) in (1x1) to G(0) in (NA1xNA2)
!
! SFEL/SFER  <------ GAAL/GAAR

!         if(IONode) write(6,*) "hola getsfe L", ipe
! LEFT:

         call getSFE(UseBulk,uGFL,HAAL,SAAL,Z,ts_ikpt,
     &        nqL,qLb,wqL,
     &        noL_GF,
     &        nuaL_GF,lastoL,NRepA1L,NRepA2L,
     &        noL,SFEL,GAMMAL,
     &        points2read(iPE),errorgs)

         if(errorgs) call die('ERROR in getSFE Left') 

! RIGHT:
         call getSFE(UseBulk,uGFR,HAAR,SAAR,Z,ts_ikpt,
     &      nqR,qRb,wqR,
     &      noR_GF,
     &      nuaR_GF,lastoR,NRepA1R,NRepA2R,
     &      noR,SFER,GAMMAR,
     &      points2read(iPE),errorgs)

         if(errorgs) call die('ERROR in getSFE Right') 

       
! Calculate Full GF at energy point
         ierr=0
         call getGF_Tri(UseBulk,updatedmcr,VoltFDF,nou,noL,noR,H,S,SFEL
     &        ,SFER,Z,GF,ierr,TriDiag)
      
#ifdef MPI
         itmp = 0
         call MPI_AllReduce(ierr,itmp,1,MPI_integer,MPI_sum,
     &        MPI_Comm_World,MPIerror)      
         ierr=itmp 
#endif

         if(ierr.ne.0) then
            if(IONode) then
               write(*,*) 'ERROR: getGF MATRIX INVERSION FAILED'
            endif 
         end if

! -------------------------------------------------------------
!     No voltage/Left Equilibrium Contour:
         if(Part .eq. 'N' .or. Part .eq. 'L') then
            if (ts_gamma_scf) then
              do i=1,nou*nou
               DM(i)  =  DM(i) - DIMAG(GF(i)*W)
               EDM(i) = EDM(i) - DIMAG(GF(i)*ZW)
              end do ! i
           else
              do i=1,nou*nou
               DMCplx(i)  =  DMCplx(i) - GF(i)*W
               EDMCplx(i) = EDMCplx(i) - GF(i)*ZW
              end do ! i
           end if
! -------------------------------------------------------------
!     Right Equilibrium Contour:
         else if(Part .eq. 'R') then
            if (ts_gamma_scf) then
               do i=1,nou*nou
                DMR(i)  =  DMR(i) - DIMAG(GF(i)*W)
                EDMR(i) = EDMR(i) - DIMAG(GF(i)*ZW)
               end do ! i
            else 
               do i=1,nou*nou
                DMRCplx(i)  =  DMRCplx(i) - GF(i)*W
                EDMRCplx(i) = EDMRCplx(i) - GF(i)*ZW
               end do ! i
            end if
! -------------------------------------------------------------
!     Non-equilibrium part:
         else if(Part .eq. 'V') then

!   Calculate rho^R(e) -> work2 ! NOTE R!!
            call GFGammaGF(.false.,nou,noR,GF,GAMMAR,work2)
! note DM = DML+DMneqR or DM=DMR+DMneqL

! note we use "+" here!
            if (ts_gamma_scf) then
               do i=1,nou*nou
                DMneqR(i) = DMneqR(i) +DIMAG(work2(i)*W)
                EDM(i)    = EDM(i) +   DIMAG(work2(i)*ZW)
               end do ! i
            else
               do i=1,nou*nou
                DMneqRCplx(i) = DMneqRCplx(i) +work2(i)*W
                EDMCplx(i)    = EDMCplx(i) +   work2(i)*ZW
               end do ! i
            end if

!  Calculate rho^L(e) -> work2 ! NOTE L!!
            call GFGammaGF(.true.,nou,noL,GF,GAMMAL,work2)
        
!     note we use "-" here!
            if (ts_gamma_scf) then
               do i=1,nou*nou
                DMneqL(i) = DMneqL(i) - DIMAG(work2(i)*W)
                EDMR(i)   = EDMR(i) +   DIMAG(work2(i)*ZW)
               end do ! i
            else
               do i=1,nou*nou
                DMneqLCplx(i) = DMneqLCplx(i) - work2(i)*W
                EDMRCplx(i)   = EDMRCplx(i) +   work2(i)*ZW
               end do
            end if
! -------------------------------------------------------------
         else
            if(IONode) 
     &           write(6,*) 'TRANSIESTA: ERROR in contour setup!'
           call die ('TRANSIESTA: ERROR in contour setup!')
         end if                 !which part of the contour ??
         

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
      END DO                    ! the loop over ALL contour points
!                               ! on each processor
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=

!     on return the result is in DM and EDM
! We need a reduction
! There is no need to initialize the work1 er work1cmplx arrays.
! The interface for the MPI routines are INTENT(OUT) which ensures
! correct handling, however, we have them still.
#ifdef MPI

! Global reduction of density matrices
      call timer("TS_comm",1)
! FDN ts_gamma_scf
      if (ts_gamma_scf) then
         allocate(work1(nou*nou))
         call memory('A','D',nou*nou,'transiesta')
         
         work1 = 0.0d0
         call MPI_AllReduce(DM,work1,nou*nou,
     &        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         DM = work1

         work1 = 0.0d0
         call MPI_AllReduce(EDM,work1,nou*nou,
     &        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         EDM = work1
      
         if (IsVolt) then

            work1 = 0.0d0
            call MPI_AllReduce(DMR,work1,nou*nou,
     &           DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
            DMR = work1

            work1 = 0.0d0
            call MPI_AllReduce(DMneqL,work1,nou*nou,
     &           DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
            DMneqL = work1

            work1 = 0.0d0
            call MPI_AllReduce(DMneqR,work1,nou*nou,
     &           DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
            DMneqR = work1

            work1 = 0.0d0
            call MPI_AllReduce(EDMR,work1,nou*nou,
     &           DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
            EDMR = work1
            
         endif
       
         call memory('D','D',size(work1),'transiesta')
         deallocate(work1)
      else

         allocate(work1Cplx(nou*nou))
         call memory('A','Z',nou*nou,'transiesta')

         work1Cplx = 0.0d0
         call MPI_AllReduce(DMCplx,work1Cplx,nou*nou,
     &        DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
         DMCplx = work1Cplx

         work1Cplx = 0.0d0
         call MPI_AllReduce(EDMCplx,work1Cplx,nou*nou,
     &        DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
         EDMCplx = work1Cplx

         if (IsVolt) then

            work1Cplx = 0.0d0
            call MPI_AllReduce(DMRCplx,work1Cplx,nou*nou,
     &           DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
            DMRCplx = work1Cplx

            work1Cplx = 0.0d0
            call MPI_AllReduce(DMneqLCplx,work1Cplx,nou*nou,
     &           DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
            DMneqLCplx = work1Cplx
            
            work1Cplx = 0.0d0
            call MPI_AllReduce(DMneqRCplx,work1Cplx,nou*nou,
     &           DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
            DMneqRCplx = work1Cplx

            work1Cplx = 0.0d0
            call MPI_AllReduce(EDMRCplx,work1Cplx,nou*nou,
     &           DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
            EDMRCplx = work1Cplx
            
         endif

         call memory('D','Z',size(work1Cplx),'transiesta')
         deallocate(work1Cplx)

      end if
      call timer("TS_comm",2)
#endif MPI

      IF ( IsVolt ) THEN 
         if (ts_gamma_scf) then
            call weightDM(nou,noL,noR,DM,DMR,DMneqL,
     &           DMneqR,EDM,EDMR)
         else
            call weightDMC(nou,noL,noR,DMCplx,DMRCplx,
     &           DMneqLCplx,DMneqRCplx,EDMCplx,EDMRCplx)
         end if
      endif
      
!--------------------------------------------------------------------
!
!     The original Hamiltonian from SIESTA was shifted ef: 
!  -- thus we have a shift in EDM:

      if (ts_gamma_scf) then 
         call daxpy(nou*nou,Ef,DM,1,EDM,1)
      else 
         call zaxpy(nou*nou,dcmplx(Ef,0.d0),DMCplx,1,EDMCplx,1)
      end if

!
!     Update Dnew, Enew written in sparse SIESTA-format:
!     
#ifdef MPI
      if(ts_gamma_scf) then
         do io = 1,no_u
            iuo = indxuo(io)
            ic = iuo - noBufL
            do j = 1,numhg(io)
               ind = listhptrg(io) + j
               if(listudg(ind)) then ! update
                  jo = listhg(ind)
                  juo = indxuo(jo)
                  jc = juo - noBufL
                  
                  Dnewg(ind,ispin) = Dnewg(ind,ispin) +
     &                 DM(ic+nou*(jc-1))
                  Enewg(ind,ispin) = Enewg(ind,ispin) +
     &                 EDM(ic+nou*(jc-1))
                  
               end if ! update

               ! Collect charges in the Q_kpt array
               call GlobalToLocalOrb(io,Node,Nodes,iio)
               if ( iio /= 0 ) then
                  Q_kpt(ts_ikpt,regiong(ind),ispin) =
     &                 Q_kpt(ts_ikpt,regiong(ind),ispin)+
     &                 Dnewg(ind,ispin)*Ss(listhptr(iio)+j)
               endif
            end do ! j
         end do    ! io
      else
         do io = 1,no_u
            iuo = indxuo(io)
            ic = iuo - noBufL
            do j = 1,numhg(io)
               ind = listhptrg(io) + j
               if(listudg(ind)) then ! update
                  jo = listhg(ind)
                  juo = indxuo(jo)
                  jc = juo - noBufL
                  
                  kxij = k(1)*xijg(1,ind) +
     &                 k(2)*xijg(2,ind) +
     &                 k(3)*xijg(3,ind)
                  cphase =  cdexp( dcmplx(0d0,-1.0d0)*kxij )
                  ccphase = cdexp( dcmplx(0d0,+1.0d0)*kxij )
                  
                  Dnewg(ind,ispin) = Dnewg(ind,ispin) + DIMAG(
     &                 dcmplx(0.5d0*ts_kweight(ts_ikpt),0d0)*
     &                 (cphase  *DMCplx(ic+nou*(jc-1)) +
     &                  ccphase *DMCplx(jc+nou*(ic-1))))
                  Enewg(ind,ispin) = Enewg(ind,ispin) + DIMAG(
     &                 dcmplx(0.5d0*ts_kweight(ts_ikpt),0d0)*
     &                 (cphase  *EDMCplx(ic+nou*(jc-1))+
     &                  ccphase *EDMCplx(jc+nou*(ic-1))))

               end if ! update

               ! Collect charges in the Q_kpt array
               call GlobalToLocalOrb(io,Node,Nodes,iio)
               if ( iio /= 0 ) then
                  Q_kpt(ts_ikpt,regiong(ind),ispin) =
     &                 Q_kpt(ts_ikpt,regiong(ind),ispin)+
     &                 Dnewg(ind,ispin)*Ss(listhptr(iio)+j)
               endif
            end do ! j
         end do    ! io
      end if
!     NPA The MPI part and the sequential part has been seperated to limit
!     memory usage!
#else // MPI >>>>>>>>>>>>>>>>>>> MPI ELSE <<<<<<<<<<<<<<<<<<<

        do io = 1,no_l
           iuo = indxuo(io)
           ic = iuo - noBufL 
           do j = 1,numh(io)
              ind = listhptr(io) + j
              jo = listh(ind)
              juo = indxuo(jo)
              jc = juo - noBufL 
              if(listud(ind)) then ! update

                 if ( ts_gamma_scf ) then

                    Dnew(ind,ispin) = Dnew(ind,ispin) +
     &                   DM(ic+nou*(jc-1))
                    Enew(ind,ispin) = Enew(ind,ispin) +
     &                   EDM(ic+nou*(jc-1))
                    

                 else
                    kxij = (k(1)*xij(1,ind) +
     &                   k(2)*xij(2,ind) +
     &                   k(3)*xij(3,ind) )
                    cphase =  exp( dcmplx(0d0,-1.0d0)*kxij )
                    ccphase = exp( dcmplx(0d0,+1.0d0)*kxij )
                    Dnew(ind,ispin) = Dnew(ind,ispin) + Dimag(
     &                   dcmplx(0.5d0*ts_kweight(ts_ikpt),0d0)*(
     &                   cphase*DMCplx(ic+nou*(jc-1))+ccphase*DMCplx(jc
     &                   +nou*(ic-1)))          )

                    Enew(ind,ispin) = Enew(ind,ispin) + Dimag(
     &                   dcmplx(0.5d0*ts_kweight(ts_ikpt),0d0)*(
     &                   cphase*EDMCplx(ic+nou*(jc-1))+ccphase
     &                   *EDMCplx(jc+nou*(ic-1)))          )
                 end if 
              end if            !update DM,EDM
              
              ! Collect charges in the Q_kpt array
              Q_kpt(ts_ikpt,region(ind),ispin) = 
     &             Q_kpt(ts_ikpt,region(ind),ispin) +
     &             Dnew(ind,ispin)*Ss(ind)
           end do               !j
       end do                   !io
#endif // MPI


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END DO                    ! ts_kpoint
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END DO                    !ispin
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      call memory('D','Z',nou*nou,'transiesta')
      deallocate(GF)
      call memory('D','Z',size(work2),'transiesta')
      deallocate(work2)

      
#ifdef MPI
      allocate(Q_kpt_loc(0:ts_nkpnt,0:9,nspin))
      call memory('A','D',(1+ts_nkpnt)*nspin*10,'transiesta')
      call MPI_AllReduce(Q_kpt,Q_kpt_loc,(1+ts_nkpnt)*nspin*10,
     &     DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
      Q_kpt = Q_kpt_loc
      call memory('D','D',(1+ts_nkpnt)*nspin*10,'transiesta')
      deallocate(Q_kpt_loc)

! Distribute Dnew/Enew 
      do io = 1,no_l
        call LocalToGlobalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(io)
               Dnew(listhptr(io)+jo,is)= Dnewg(listhptrg(iio)+jo,is)
               Enew(listhptr(io)+jo,is)= Enewg(listhptrg(iio)+jo,is)
            enddo
          enddo
      enddo
      call memory('D','D',size(Dnewg)+size(Enewg),'transiesta')
      deallocate(Dnewg)
      deallocate(Enewg)

#endif

      ! Charge in device
      Qc = sum(Q_kpt(ts_nkpnt,5,:))

! 
!     Charge Renormalization 
!
      if ( ChargeCorr == 1 ) then ! Correct charges in the buffer region
         allocate(reD(2))
         call memory('A','D',2,'transiesta')
         reD(1) = REAL(noBufL+noBufR,dp)/REAL(noBufL,dp)
         reD(1) = (Qtot-sum(Q_kpt(ts_nkpnt,:,:)))/(reD(1)
     &        *sum(Q_kpt(ts_nkpnt,1,:)))
         reD(2) = REAL(noBufL+noBufR,dp)/REAL(noBufR,dp)
         reD(2) = (Qtot-sum(Q_kpt(ts_nkpnt,:,:)))/(reD(2)
     &        *sum(Q_kpt(ts_nkpnt,9,:)))
         ! Apply charge-correction factor 
         ! This will reduce "heavy" charge fluctuations and
         ! should guard against this.
         reD(:) = reD(:) * ChargeCorr_factor
         do ispin = 1 , nspin
            do io = 1,no_l
               do j = 1,numh(io)
                  ind = listhptr(io) + j
                  if ( region(ind) == 1 ) then
                     Dnew(ind,ispin) = Dnew(ind,ispin)*reD(1)+Dnew(ind
     &                    ,ispin)
                     Enew(ind,ispin) = Enew(ind,ispin)*reD(1)+Enew(ind
     &                    ,ispin)
                  else if ( region(ind) == 9 ) then
                     Dnew(ind,ispin) = Dnew(ind,ispin)*reD(2)+Dnew(ind
     &                    ,ispin)
                     Enew(ind,ispin) = Enew(ind,ispin)*reD(2)+Enew(ind
     &                    ,ispin)
                  end if
               end do           ! numh
            end do              ! no_l
         end do                 !ispin
      end if

      if(IONode) then
         write(*,'(a,I6)') 'TranSIESTA charge on k-points:'
     &        ,ts_nkpnt
         if ( nspin == 2 ) then
            if ( mod(ts_nkpnt,2) == 0 ) then
               io = ts_nkpnt / 2 
            else
               io = (ts_nkpnt+2-mod(ts_nkpnt,2)) / 2
            end if
            do i = 1 , io
               jo = 1
               if ( i + jo * io > ts_nkpnt ) then
                  jo = 0
               end if
               write(*,'(tr2,2(i3,2(tr1,f8.3),tr3))') (i+io*j,
     &              sum(Q_kpt(i+io*j,:,:),1)-sum(Q_kpt(i+io*j-1,:,:),1),
     &              sum(Q_kpt(i+io*j,:,:),2)-sum(Q_kpt(i+io*j-1,:,:),2),
     &              j=0,jo)
            end do
         else
            if ( mod(ts_nkpnt,4) == 0 ) then
               io = ts_nkpnt / 4 
            else
               io = (ts_nkpnt+4-mod(ts_nkpnt,4)) / 4 
            end if
            do i = 1 , io
               jo = 3
               do j = 0 , 3 
                  if ( i + j * io > ts_nkpnt ) then
                     jo = j - 1
                     exit
                  end if
               end do
               write(*,'(tr2,4(i3,tr1,f8.3,tr3))') (i+io*j,
     &              sum(Q_kpt(i+io*j,:,:),1)-sum(Q_kpt(i+io*j-1,:,:),1),
     &              j=0,jo)
            end do
         end if
         write(*,'(a)') "TranSIESTA region charges:"
         write(*,'(6(a9,1X))') 'U','L','L-C','C','C-R','R'
         if ( ChargeCorr == 1 ) then
            write(*,'(2(a9,1X))') 'Left Buf','Right Buf'
         end if
         write(*,'(6(f9.3,1X))') sum(Q_kpt(ts_nkpnt,0,:))
     &        ,sum(Q_kpt(ts_nkpnt,3,:)),sum(Q_kpt(ts_nkpnt,4,:))
     &        ,sum(Q_kpt(ts_nkpnt,5,:)),sum(Q_kpt(ts_nkpnt,6,:))
     &        ,sum(Q_kpt(ts_nkpnt,7,:))
         if ( ChargeCorr == 1 ) then
            write(*,'(2(f9.3,1X))') sum(Q_kpt(ts_nkpnt,1,:))
     &           +reD(1)*noBufL,sum(Q_kpt(ts_nkpnt,9,:))+reD(2)*noBufR
            call memory('D','D',2,'transiesta')
            deallocate(reD)
         end if

         write(*,'(a,1X,5(f9.3,1X))') 
     &        'TranSIESTA: Qc, Qt, Qt0 : ',Qc,Qc+Qcn0,Qtot
      endif 
!--------------------------------------------------------------------

!     Close Files

      if(IONode) then
         call io_close(uGFL)
         call io_close(uGFR)
      endif
      
!--------------------------------------------------------------------
!     Deallocate:
      call memory('D','D',(1+ts_nkpnt)*nspin*10,'transiesta')
      deallocate(Q_kpt)
      
      call memory('D','D',4*nqL,'transiesta')
      deallocate(qLb,wqL)
      call memory('D','D',4*nqR,'transiesta')
      deallocate(qRb,wqR)
      call memory('D','D',4*nkparL,'transiesta')
      deallocate(kparL,wkparL)
      call memory('D','D',4*nkparR,'transiesta')
      deallocate(kparR,wkparR)

!     Deallocate:
      call memory('D','Z',size(H)+size(S),'transiesta')

      deallocate(H)
      deallocate(S)

      if ( ts_gamma_scf ) then
         call memory('D','D',size(DM)+size(EDM),'transiesta')
         deallocate(DM)
         deallocate(EDM)
      else
         call memory('D','Z',size(DMCplx)+size(EDMCplx),'transiesta')
         deallocate(DMCplx)
         deallocate(EDMCplx)
      end if

      call memory('D','Z',size(SFEL)+size(SFER),'transiesta')
      deallocate(SFEL)
      deallocate(SFER)
      call memory('D','Z',size(GAMMAL)+size(GAMMAR),'transiesta')
      deallocate(GAMMAL)
      deallocate(GAMMAR)

      call memory('D','Z',size(HAAL)+size(SAAL),'transiesta')
      call memory('D','Z',size(HAAR)+size(SAAR),'transiesta')
      deallocate(HAAL)
      deallocate(SAAL)
      deallocate(HAAR)
      deallocate(SAAR)

      
      if(IsVolt) then
         if ( ts_gamma_scf ) then
            call memory('D','D',size(DMR)+size(DMneqL),'transiesta')
            call memory('D','D',size(EDMR)+size(DMneqR),'transiesta')
            deallocate(DMR)
            deallocate(DMneqL)
            deallocate(DMneqR)
            deallocate(EDMR)
         else
            call memory('D','Z',size(DMRCplx)+size(DMneqLCplx)
     &           ,'transiesta')
            call memory('D','Z',size(DMneqRCplx)+size(EDMRCplx)
     &           ,'transiesta')
            deallocate(DMRCplx)
            deallocate(DMneqLCplx)
            deallocate(DMneqRCplx)
            deallocate(EDMRCplx)
         end if
      endif


#ifdef MPI
      call memory('D','I',size(numhg)+size(listhptrg)+size(listhg)
     &     ,'transiesta')
      deallocate(numhg)
      deallocate(listhptrg)
      deallocate(listhg)
      if(.not.gamma) then
         call memory('D','D',size(xijg),'transiesta')
         deallocate(xijg)
      endif
      call memory('D','D',maxnhg*2,'transiesta')
      deallocate(Hsg,Ssg)
#endif

      call timer("TS_calc",2)
! ===================================================================
      firsttime = .false.       !first time -- not any more!
#ifdef DEBUG
      call write_debug( 'POS transiesta' )
#endif
      end subroutine transiesta


!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      subroutine TSVHFix(  mesh, nsm, v )
!
!  Modules
!
      use precision, only : dp, grid_p
      use sys, only : die
! FDN put the only attribute
      use parallel, only : ProcessorY, Node, Nodes
! FDN
#ifdef MPI
      use mpi_siesta, only : MPI_AllReduce, MPI_Comm_Rank, MPI_Sum
      use mpi_siesta, only : DAT_double, MPI_Integer, MPI_Comm_World
! FDN parallelsubs
      use parallelsubs, only : HowManyMeshPerNode
! FDN
#endif

      implicit          none
      real(grid_p)      v(*)
      integer           mesh(3), nsm


! Internal variables
      integer           i1, i2, i3, imesh, ntemp
      integer           nlp, meshl(3)
      integer           ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset
      integer           Zoffset, Py, Pz, i30, i20, meshnsm(3)
      integer           NRemY, NRemZ
#ifdef MPI
      integer           MPIerror, npl
#endif
      real(dp)          vav, vtot, temp

#ifdef DEBUG
      call write_debug( 'PRE TSVHfix' )
#endif

      call timer('TSVHFix',1)

! Find local number of mesh points
        meshnsm(1) = mesh(1)/nsm
        meshnsm(2) = mesh(2)/nsm
        meshnsm(3) = mesh(3)/nsm
#ifdef MPI
        !! print *, "N:",Node, "TSVHFix-- mesh:", mesh
        call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
        meshl(1) = meshl(1)*nsm
        meshl(2) = meshl(2)*nsm
        meshl(3) = meshl(3)*nsm
#else
        meshl(1) = mesh(1)
        meshl(2) = mesh(2)
        meshl(3) = mesh(3)
#endif

! Check that ProcessorY is a factor of the number of processors
        if (mod(Nodes,ProcessorY).gt.0) then
           call die('ERROR: ProcessorY must be a factor of the'//
     &          ' number of processors!')
        endif
        ProcessorZ = Nodes/ProcessorY

! Calculate blocking sizes
        BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
        NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
        BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
        NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

! Calculate coordinates of current node in processor grid
        Py = (Node/ProcessorZ)+1
        Pz = Node - (Py - 1)*ProcessorZ + 1

! Calculate starting point for grid
        Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
        Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

         vtot =0.d0
         nlp =0

        imesh = 0
        i30 = Zoffset - 1
        do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            do i1 = 0,meshl(1)-1
              imesh = imesh + 1
              if (i30.eq.0) then
                 nlp = nlp +1 
                 vtot = vtot + v(imesh)
              endif
            enddo
          enddo
        enddo

#ifdef MPI
      temp=0.d0
        call MPI_AllReduce(vtot,temp,1,DAT_double,MPI_Sum,
     &  MPI_Comm_World,MPIerror)
       vtot = temp
       ntemp=0
        call MPI_AllReduce(nlp,ntemp,1,MPI_integer,MPI_Sum,
     &    MPI_Comm_World,MPIerror)
       nlp = ntemp
 
#endif
        vav=vtot/float(nlp)

        imesh = 0
        
        do i3 = 0,meshl(3)-1
           do i2 = 0,meshl(2)-1
              do i1 = 0,meshl(1)-1
                 imesh = imesh + 1
                 v(imesh) = v(imesh) - vav
              enddo
           enddo
        enddo

        call timer('TSVHFix',2)

#ifdef DEBUG
      call write_debug( 'POS TSVHfix' )
#endif

      end subroutine TSVHFix

!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      END MODULE m_ts_in_siesta
