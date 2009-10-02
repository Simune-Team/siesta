MODULE m_ts_io
!
! Routines that are used for Input and Output of files 
!
!=============================================================================
! CONTAINS:
!          1) read_green
!          2) ts_iohs
!          3) TSiodm


  implicit none

  public :: read_green, ts_iohs, TSiodm

  private

  CONTAINS


! ##################################################################
! ##            read-in header of Greens function file            ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ## Changed to F90 by Jose-Luis Mozos, jlm@icmab.es              ##
! ##################################################################

      subroutine read_green(jgfu,EFermi,NEn,contour,wgf, &
                              nua,NA1,NA2,nq,ng,wq,q,errorgf)


      implicit none

      real*8 EPS
      parameter(EPS=1d-7)

      logical PRINTALOT
!      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)

! INPUT
      integer jgfu               !unit of gf-file


! these are the values expected:
! they will be compared to the read-in and ERROR messages will appear
      real*8 efermi             !The required Fermi energy
      integer NEn,ng
      complex*16 contour(NEn),wgf(NEn)
      integer nua,NA1,NA2       !no. atoms in uc and no. repetitions in A1,A2
      integer nq                !no. q-points
      logical tleft

      logical errorgf


! READ-IN values
      real*8 efermii            !The required Fermi energy
      integer NEni,ngi
      complex*16, dimension(:), allocatable :: contouri,wgfi

      integer nqi

!     q-point=k_|| point and weigth

      real*8, dimension (:,:), pointer:: q
      real*8, dimension (:), pointer:: wq

      character*70 gftitle

! Helpers..

      complex*16 ctmp
      integer iEn


!=======================================================================
! BEGIN:
!=======================================================================

       errorgf = .false.       


      read(jgfu) gftitle
      read(jgfu) EFermii,NEni

      read(jgfu) nua,NA1,NA2,nqi

!      write(6,*) 'read GF: ',gftitle

      if(NEni .ne. NEn)  then
         write(6,*) 'read GF: ERROR: NEn=',NEni,' expected:', NEn
            errorgf = .true. 
            return
      end if


      allocate(contouri(nen))
      allocate(wgfi(nen))

      nq = nqi

! FDN in general q// may have z component
!     dimension changed from 2 to 3
      allocate(q(3,nqi))
! FDN
      allocate(wq(nqi))

      read(jgfu) contouri,wgfi,q,wq
      read(jgfu) ng


!     check on contours.
      do iEn=1,NEn

         ctmp=contouri(iEn)-contour(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(6,*) ' Warning: contours differ by >', EPS
         end if
         if(cdabs(ctmp).GT.10d0*EPS) then 
            write(6,*) &
                 ' ERROR: contours differ by >', 10.d0*EPS
            errorgf = .true. 
            return
         end if
      end do

      do iEn=1,NEn
         ctmp=wgfi(iEn)-wgf(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(6,*)  &
                ' ERROR: contour weights differ by >',EPS
            errorgf = .true. 
            return
         end if
      end do

      deallocate(contouri)
      deallocate(wgfi)
      




      return
      end subroutine read_green





!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine ts_iohs(task, gamma, onlyS, no_u, no_s, Enspin,   &
      		 indxuo, maxnh, numh, listhptr, listh, H, S, qtot, &
		 temp, xij, fnlength, fname, na_u, lasto, isa, ef, &
		 ucell, ts_kscell_file, ts_kdispl_file,        	   &
		 ts_gamma_scf_file, xa, istep, ia1)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
! logical       gamma         : Is only gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,Enspin)     : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

!
!  Modules
!
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb, &
                               GlobalToLocalOrb, GetNodeOrbs
      use fdf
      use files,        only : slabel, label_length
      use sys,          only : die
! TSS Begin
! Added variables
!      use atomlist, only : lasto
!      use m_energies, only : ef
!      use m_spin, only : efs
!      use siesta_geom, only : xa, isa, ucell
      use m_ts_kpoints, only: ts_gamma_scf, ts_kscell, ts_kdispl
! TSS End
#ifdef MPI
      use mpi_siesta
#endif
! TSS Begin
! SIESTA Modules
      use m_spin, only: nspin       
! TSS End

      implicit          none

      character(8) :: task
      character(len=label_length+3) :: paste
      logical           gamma, onlyS, onlySfile
      integer           maxnh, no_u, no_s, Enspin
! TSS
!      integer           indxuo(no_s), listh(maxnh), numh(*), listhptr(*)
      integer, dimension(:), pointer :: indxuo, listh, numh, listhptr 
! TSS
!      real(dp)          H(maxnh,Enspin), S(maxnh), &
!                        qtot, temp, xij(3,maxnh)
      real(dp), dimension(:,:), pointer :: H, xij
      real(dp), dimension(:), pointer :: S
      real(dp) :: qtot, temp 
      integer :: istep, ia1
      external          io_assign, io_close, paste
! TSS Begin
! Added arguments for TRANSIESTA
      integer :: fnlength, na_u
      character(fnlength) :: fname
      integer, dimension(:), pointer :: lasto, isa
      integer, dimension(3,3) :: ts_kscell_file
      real(dp) :: ef
      real(dp), dimension(3,3) :: ucell
      real(dp), dimension(:,:), pointer :: xa
      real(dp), dimension(3) :: ts_kdispl_file 
      logical :: ts_gamma_scf_file
! TSS End

! Internal variables and arrays
      integer    im, is, iu, ju, k, mnh, ns
      integer    ih,hl,nuo,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
      real(dp), dimension(:,:), allocatable :: buffer2
#endif
      logical    baddim, found, gammaonfile
! TSS Begin
! Aux Variables
      integer :: ispin, ios, i, j, i2, l
! TSS End



! #################### READ ######################################
! Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then
    

! Check if input file exists
          inquire( file=fname, exist=found)

        if (found) then

! Open file
            call io_assign( iu )
            open( iu, file=fname, form='unformatted', status='old' )      

! Read dimensions
            read(iu) na_u, no_u, no_s, Enspin, maxnh

! Allocations



! Check dimensions
          baddim = .false.
          if (Enspin  .ne. nspin) baddim = .true.
           
          if (baddim) then
            if (Node.eq.0) then
              call io_assign( ju )
              open( ju, file='ts_iohs.h', status='unknown' )
              write(ju,'(a)') 'C Dimensions for input to ts_iohs'
              write(ju,'(6x,a,i8,a)') 'parameter ( nspin =', ns,  ' )'
              call io_close( ju )
              call die('ts_iohs: BAD DIMENSIONS')
            else
              call die()
            endif
          endif

! TSS Begin
! Allocate arrays that are going to be read now
          nullify(xa,isa)
          allocate(xa(3,na_u))
          allocate(isa(na_u)) 

! Read Geometry information
           read(iu) xa
           read(iu) isa   
       	   read(iu) ucell  

! Read k-point sampling information
           read(iu) gammaonfile   
	   read(iu) onlySfile
           read(iu) ts_gamma_scf_file       
           read(iu) ts_kscell_file
           read(iu) ts_kdispl_file  
	   read(iu) istep, ia1

! Check whether file is compatible from a gamma point of view
           if (onlySfile) then
             write(*,*) 'TSHS file does not contain Hamiltonian'
             call die('ts_iohs: onlyS flag, no H in file')
           endif
           if (.not.gamma.and.gammaonfile) then
             write(*,*) 'System Gamma:',gamma
             write(*,*) 'Electrode Gamma:',gammaonfile
             call die('ts_iohs: Non-TSgamma information not present')
           endif
           if (.not.ts_gamma_scf .and. ts_gamma_scf_file) then
             write(*,*) 'TS_Electrode Gamma:',ts_gamma_scf_file
             write(*,*) 'TS_System Gamma:',ts_gamma_scf
             call die('ts_iohs: Incompatible k samplings!') 
           end if
           do i = 1,3
             if ( ts_kdispl_file(i) /= ts_kdispl(i) ) then
               write(*,'(a,2F8.4)') 'TS_Electrode k displacements:',(ts_kdispl_file(j),j=1,2)
               write(*,'(a,2F8.4)') 'TS_System k displacements:', (ts_kdispl(j),j=1,2)
               call die('ts_iohs: Incompatibles k displacements') 
             end if
           enddo
           do i = 1,3
             do j = 1,3
               if ( ts_kscell(i,j)  /= ts_kscell_file(i,j) ) then
                 do i2=1,2
                   write(*,*) 'TS_Electrode k points:', (ts_kscell_file(i2,l),l=1,3)
                 enddo 
                 do i2=1,2
                   write(*,*) 'TS_System k points:', (ts_kscell(i2,l),l=1,3)
                 enddo
                 call die('ts_iohs: Incompatible k-point samplings')
               end if
             enddo
           enddo


! Read sparse listings
          nullify(lasto)
          allocate(lasto(0:na_u))
            read(iu) lasto


          if (.not.gamma) then
! Allocate arrays that are going to be read now
            nullify(indxuo)
            allocate(indxuo(1:no_s))
              read(iu) (indxuo(ih),ih=1,no_s)

          endif

! Allocate local array for global numh
          nullify(numh)
          allocate(numh(no_u))
          call memory('A','I',no_u,'iohs')

! Read numh and send to appropriate Node
	    read(iu) numh(1:no_u)



! Read Electronic Structure Information
            read(iu) qtot,temp
            read(iu) ef


! Create listhptr
          nullify(listhptr)
          allocate(listhptr(no_u))
          listhptr(1) = 0
! TSS nuo->no_u
! TSS hl->hi
          do ih = 2,no_u
            listhptr(ih) = listhptr(ih-1) + numh(ih-1)
          enddo

! Read listh
! Allocate lish
          nullify(listh)
          allocate(listh(maxnh))

           do ih = 1,no_u
             read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih))
           enddo

! Read Overlap matrix
! Allocate S
          nullify(S)
          allocate(S(maxnh))
           do ih = 1,no_u
             read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
           enddo

! Read Hamiltonian
! Allocate H
          nullify(H)
          allocate(H(maxnh,Enspin))
           do is = 1,Enspin
             do ih = 1,no_u
               read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
             enddo
           enddo

          if (.not.gamma) then
! Read interorbital vectors for K point phasing
! Allocate xij
            nullify(xij)
            allocate(xij(3,maxnh))
             do ih = 1,no_u
               read(iu) (xij(k,listhptr(ih)+1:listhptr(ih)+numh(ih)),k=1,3)
             enddo

          endif


! Close file
             call io_close( iu )

        else
            write(*,*) 'iohs: ERROR: file not found: ', fname
            call die('iohs: ERROR: file not found')
        endif
! FDN Temp Fim
! #################### WRITE ######################################
      elseif (task.eq.'write' .or. task.eq.'WRITE') then

! Find total numbers over all Nodes
#ifdef MPI
         call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum, &
         MPI_Comm_World,MPIerror)
#else
         maxnhtot = maxnh
#endif

        if (Node.eq.0) then
! Open file
          call io_assign( iu )
          open( iu, file=fname, form='unformatted', status='unknown' )      
! Write Dimensions Information
          write(iu) na_u, no_u, no_s, Enspin, maxnhtot

! Write Geometry information
          write(iu) xa(1:3,1:na_u)
          write(iu) isa(1:na_u)          
          write(iu) ucell

! Write k-point samplung information
          write(iu) gamma
          write(iu) onlyS
          write(iu) ts_gamma_scf_file
          write(iu) ts_kscell_file
          write(iu) ts_kdispl_file
	  write(iu) istep, ia1

! Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

! Write sparse listings
            write(iu) lasto(0:na_u)

          if (.not.gamma) then
            write(iu) (indxuo(ih),ih=1,no_s)
          endif

        endif

! Create globalised numh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            numhg(ih) = numh(hl)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(numh(hl),1,MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numhg(ih),1,MPI_integer, &
              BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

        if (Node.eq.0) then
! Write numh ( Sparse listings )
          maxhg = 0
          do ih = 1,no_u
            maxhg = max(maxhg,numhg(ih))
          enddo
	  write(iu) numhg(1:no_u)
#ifdef MPI
          allocate(buffer(maxhg))
          call memory('A','D',maxhg,'iohs')
          allocate(ibuffer(maxhg))
          call memory('A','I',maxhg,'iohs')
#endif



! Write Electronic Structure Information
          write(iu) qtot,temp
          write(iu) ef
        endif
! Write listh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
	    write(iu) listh(listhptr(hl)+1:listhptr(hl)+numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numhg(ih),MPI_integer,BNode,1, &
             MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(listh(listhptr(hl)+1),numh(hl),MPI_integer, &
             0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
                write(iu) ibuffer(1:numhg(ih))
            endif
          endif
#endif
        enddo
#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'iohs')
          deallocate(ibuffer)
        endif
#endif

! Write Overlap matrix
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
	    write(iu) S(listhptr(hl)+1:listhptr(hl)+numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
             BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(S(listhptr(hl)+1),numh(hl), &
             MPI_double_precision,0,1,MPI_Comm_World, &
             Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
	      write(iu) buffer(1:numhg(ih))
            endif
          endif
#endif
        enddo

	if (.not. onlyS) then
! Write Hamiltonian	 
        do is=1,Enspin	 
          do ih=1,no_u	 
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu) H(listhptr(hl)+1:listhptr(hl)+numh(hl),is)
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(H(listhptr(hl)+1,is),numh(hl), &
               MPI_double_precision,0,1,MPI_Comm_World, &
               Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
	        write(iu) buffer(1:numhg(ih))
              endif
            endif
#endif
          enddo
        enddo
	endif  ! onlyS

#ifdef MPI
          if (Node .eq. 0) then
! Free buffer array
             call memory('D','D',size(buffer),'iohs')
             deallocate(buffer)
          endif
#endif


        if (.not.gamma) then
#ifdef MPI
! Allocate buffer array
          if (Node .eq. 0) then
             allocate(buffer2(3,maxhg))
             call memory('A','D',3*maxhg,'iohs')
          endif
#endif
          do ih = 1,no_u
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu) (xij(k,listhptr(hl)+1:listhptr(hl)+numh(hl)),k=1,3)
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer2(1,1),3*numhg(ih), &
               MPI_double_precision,BNode,1,MPI_Comm_World, &
               Request,MPIerror) 
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(xij(1,listhptr(hl)+1),3*numh(hl), &
               MPI_double_precision,0,1,MPI_Comm_World, &
               Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                write(iu) (buffer2(k,1:numhg(ih)),k=1,3)
              endif
            endif
#endif
          enddo
#ifdef MPI
          if (Node .eq. 0) then
! Free buffer array
             call memory('D','D',size(buffer2),'iohs')
             deallocate(buffer2)
          endif
#endif
        endif   ! not gamma

        if (Node.eq.0) then
! Deallocate local array for global numh
          call memory('D','I',size(numhg),'iohs')
          deallocate(numhg)
! Close file
          call io_close( iu )
        endif

      endif


      end subroutine ts_iohs

!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


      subroutine TSiodm( task, maxnd, nbasis, nspin, numd, & 
                    listdptr, listd, dm, edm, ef, found )
! *******************************************************************
! Reads/writes density matrix from/to file
! Written by P.Ordejon and J.M.Soler. May 1997.
! ********* INPUT ***************************************************
! character task*3 : 'read' or 'write'
! integer   maxnd  : First dimension of listd and dm
! integer   nbasis : Number of atomic orbitals
! integer   nspin  : Number of spins (1 or 2)
! ********* INPUT OR OUTPUT (depending on task) *********************
! integer numd(nbasis)     : !ontrol vector of DM matrix
!                            (number of nonzero elements of each row)
! integer listdptr(nbasis) : !ontrol vector of DM matrix
!                            (pointer to the start of each row)
! integer listd(maxnd)     : !ontrol vector of DM matrix
!                            (list of nonzero elements of each row)
! real*8  dm(maxnd,nspin)  : Density matrix
! real*8  edm(maxnd,nspin) : Energy Density matrix
! real*8  ef,efs(nspin)    : Fermi energies
! ********* OUTPUT *************************************************
! logical found : Has DM been found in disk? (Only when task='read')
! ******************************************************************
      
!
!  Modules
!
      use precision, only : sp, dp
!      use parallel
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

! TSS Begin Added
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb, &
                              LocalToGlobalOrb
! TSS End

      implicit  none

      character task*(*), paste*33
      integer   maxnd, nbasis, nspin
      integer   listd(maxnd), numd(nbasis), listdptr(nbasis)
      real*8    dm(maxnd,nspin)
      real*8    edm(maxnd,nspin)
      real*8    ef
      logical, optional :: found

! Internal variables and arrays
      character fname*33, sname*30
      logical   exist1, exist2, exist3, frstme
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      integer   Node, Nodes, nbasistot, ml, ndmaxg
      integer, dimension(:), allocatable, save :: numdg
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real*8, dimension(:), allocatable, save :: buffer
      integer, dimension(:), allocatable, save :: ibuffer
#endif
      external          chkdim, io_assign, io_close, paste, timer, &
                        memory

      save      frstme, fname
      data      frstme /.true./

!     call timer( 'iodm', 1 )

! Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif
! Find file name
      if (frstme) then
        if (Node.eq.0) then
          sname = fdf_string('SystemLabel','siesta')
        endif
#ifdef MPI
        call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World, &
          MPIerror)
#endif
        fname = paste(sname,'.TSDE')
        frstme = .false.
      endif

! Find total number of basis functions over all Nodes
#ifdef MPI
      if (Nodes.gt.1) then
        call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum, &
          MPI_Comm_World,MPIerror)
      else
        nbasistot = nbasis
      endif
#else
      nbasistot = nbasis
#endif


! Allocate local buffer array for globalised numd
      allocate(numdg(nbasistot))
      call memory('A','I',nbasistot,'iodm')

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file='SAVE.TSDM',     exist=exist1)
          inquire (file='SAVE.ctrlTSDM', exist=exist2)
          inquire (file=fname,         exist=exist3)
        endif

#ifdef MPI
! Broadcast logicals so that all processors take the same route
        call MPI_Bcast(exist1,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist2,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

        if (exist1 .and. exist2) then
! Old-format files
         if (Node.eq.0) then
            write(6,'(/,a)')'TSiodmS: Reading Density Matrix from files'
            call io_assign(unit1)
            call io_assign(unit2)
            open( unit1, file='SAVE.DM', status='unknown')
            open( unit2, file='SAVE.ctrlDM', status='unknown')
            rewind(unit1)
            rewind(unit2)
          endif

          if (Node.eq.0) then
            read(unit2,*) (numdg(m),m=1,nbasistot)
          endif
#ifdef MPI
          call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World, &
            MPIerror)
#endif

! Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
            call LocalToGlobalOrb(m,Node,Nodes,ml)
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

! Check size of first dimension of dm
          call chkdim( 'TSiodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
! Create buffer arrays for transfering density matrix between nodes and lists
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'TSiodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'TSiodm')
#endif

          do m = 1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              read(unit2,*) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
            elseif (Node.eq.0) then
              do im = 1,numdg(m)
                read(unit2,*) ibuffer(im)
              enddo
              call MPI_ISend(ibuffer,numdg(m),MPI_integer, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_IRecv(listd(listdptr(ml)+1),numd(ml), &
                MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
          enddo

#ifdef MPI
          call memory('D','I',size(ibuffer),'TSiodm')
          deallocate(ibuffer)
#endif

          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                do im = 1,numd(ml)
                  read(unit1,*) dm(listdptr(ml)+im,is)
                enddo

#ifdef MPI
              elseif (Node.eq.0) then
                do im = 1,numdg(m)
                  read(unit1,*) buffer(im)
                enddo
                call MPI_ISend(buffer,numdg(m),DAT_double, &
                  BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml), &
                  DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo


#ifdef MPI
! Free buffer array
          call memory('D','D',size(buffer),'TSiodm')
          deallocate(buffer)
#endif
          if (Node.eq.0) then
            call io_close(unit1)
            call io_close(unit2)
          endif

          if(present(found)) found = .true.

        elseif (exist3) then
! New-format files
          if (Node.eq.0) then
            write(6,'(/,a)')'TSiodmS: Reading Density Matrix from files'
            call io_assign(unit1)
            open( unit1, file=fname, &
                form='unformatted', status='unknown' )
            rewind(unit1)
            read(unit1) nb, ns
          endif

! Communicate the values to all Nodes and adjust to allow for
! distributed memory before checking the dimensions
#ifdef MPI
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

! Check dimensions
          call chkdim( 'TSiodm', 'nbasis', nbasistot, nb, 0 )
          call chkdim( 'TSiodm', 'nspin',  nspin,  ns, 0 )

          if (Node.eq.0) then
            read(unit1) (numdg(m),m=1,nbasistot)
          endif
#ifdef MPI
          call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World, &
            MPIerror)
#endif

! Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
            call LocalToGlobalOrb(m,Node,Nodes,ml)
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

! Check size of first dimension of dm
          call chkdim( 'TSiodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
! Create buffer arrays for transfering density matrix between nodes and lists
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'TSiodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'TSiodm')
#endif

          do m = 1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
            elseif (Node.eq.0) then
              read(unit1) (ibuffer(im),im=1,numdg(m))
              call MPI_ISend(ibuffer,numdg(m),MPI_integer, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_IRecv(listd(listdptr(ml)+1),numd(ml), &
                MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
          enddo

#ifdef MPI
          call memory('D','I',size(ibuffer),'TSiodm')
          deallocate(ibuffer)
#endif

! read DM

          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                read(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
              elseif (Node.eq.0) then
                read(unit1) (buffer(im),im=1,numdg(m))
                call MPI_ISend(buffer,numdg(m),DAT_double, &
                  BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml), &
                  DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

! read EDM

          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
               read(unit1)(edm(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
              elseif (Node.eq.0) then
                read(unit1) (buffer(im),im=1,numdg(m))
                call MPI_ISend(buffer,numdg(m),DAT_double, &
                  BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(edm(listdptr(ml)+1,is),numd(ml), &
                  DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

#ifdef MPI
! Free buffer array
          call memory('D','D',size(buffer),'TSiodm')
          deallocate(buffer)
#endif
          if (Node.eq.0) then
            read(unit1) ef
            call io_close(unit1)
          endif

          if(present(found)) found = .true.

        else

          if(present(found)) found = .false.

        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file=fname, form='unformatted', &
            status='unknown' )
          rewind(unit1)
          write(unit1) nbasistot, nspin
        endif

! Create globalised numd
        do m = 1,nbasistot
#ifdef MPI
          call WhichNodeOrb(m,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            numdg(m) = numd(ml)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(numd(ml),1,MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numdg(m),1,MPI_integer, &
              BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

! Write out numd array
        if (Node.eq.0) then
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo
          write(unit1) (numdg(m),m=1,nbasistot)
#ifdef MPI
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'TSiodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'TSiodm')
#endif
        endif

! Write out listd array
        do m = 1,nbasistot
#ifdef MPI
          call WhichNodeOrb(m,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1, &
              MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
              write(unit1) (ibuffer(im),im=1,numdg(m))
            endif
          endif
#endif
        enddo

#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'TSiodm')
          deallocate(ibuffer)
        endif
#endif

! Write density matrix
        do is=1,nspin
          do m=1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              write(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numdg(m),DAT_double, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_ISend(dm(listdptr(ml)+1,is),numd(ml),DAT_double, &
                0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                write(unit1) (buffer(im),im=1,numdg(m))
              endif
            endif
#endif
          enddo
        enddo

! Write energy density matrix
        do is=1,nspin
          do m=1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              write(unit1) (edm(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numdg(m),DAT_double, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
             call MPI_ISend(edm(listdptr(ml)+1,is),numd(ml),DAT_double, &
               0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                write(unit1) (buffer(im),im=1,numdg(m))
              endif
            endif
#endif
          enddo
        enddo


        if (Node.eq.0) then
#ifdef MPI
          call memory('D','D',size(buffer),'TSiodm')
          deallocate(buffer)
#endif
          write(unit1) ef

          call io_close(unit1)
        endif

      else
        if (Node.eq.0) then
          stop 'TSiodmS: incorrect task'
        else
          stop
        endif
      endif

! Deallocate local buffer array for globalised numd
      call memory('D','I',size(numdg),'iodm')
      deallocate(numdg)


!     call timer( 'iodm', 2 )
      end subroutine TSiodm








END MODULE m_ts_io
