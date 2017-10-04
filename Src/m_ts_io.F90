module m_ts_io
!
! Routines that are used for Input and Output of files 
!
!=============================================================================
! CONTAINS:
!          1) ts_read_TSHS_na
!          2) ts_iohs


  implicit none

  public :: ts_iohs
  public :: ts_read_TSHS_na
  public :: ts_read_TSHS_lasto

  private

contains

  ! Reads in the number of atoms in the electrode. This is for easy operation
  ! In the options reading phase. We need the size of the electrodes to determine 
  ! the number of atoms in the electrode.
  subroutine ts_read_TSHS_na(TSHS,na_u)
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=200), intent(in) :: TSHS

! ***********************
! * OUTPUT variables    *
! ***********************    
    integer, intent(out)           :: na_u

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uTSHS,tmp(4)
#ifdef MPI
    integer :: MPIerror
#endif
    
    if ( IONode ) then
       call io_assign(uTSHS)
       open(file=TSHS,unit=uTSHS,form='unformatted')
       read(uTSHS) na_u, tmp(1:4) !na_u, no_u, no_s, Enspin, maxnh
       call io_close(uTSHS)
    end if

#ifdef MPI
    call MPI_Bcast(na_u,1,MPI_INTEGER,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_read_TSHS_na


  ! Reads in the orbitals in the electrode. This is for easy operation
  subroutine ts_read_TSHS_lasto(TSHS,na_u,lasto)
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=200), intent(in) :: TSHS
    integer,            intent(in) :: na_u

! ***********************
! * OUTPUT variables    *
! ***********************    
    integer,           intent(out) :: lasto(0:na_u)

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uTSHS
#ifdef MPI
    integer :: MPIerror
#endif
    
    if ( IONode ) then
       call io_assign(uTSHS)
       open(file=TSHS,unit=uTSHS,form='unformatted')
       read(uTSHS) ! na_u, no_u, no_s, Enspin, maxnh
       read(uTSHS) ! xa
       read(uTSHS) ! isa   
       read(uTSHS) ! ucell  
       read(uTSHS) ! gammaonfile   
       read(uTSHS) ! onlySfile
       read(uTSHS) ! ts_gamma_scf_file       
       read(uTSHS) ! ts_kscell_file
       read(uTSHS) ! ts_kdispl_file  
       read(uTSHS) ! istep, ia1
       read(uTSHS) lasto

       call io_close(uTSHS)
    end if

#ifdef MPI
    call MPI_Bcast(lasto,na_u+1,MPI_INTEGER,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_read_TSHS_lasto


! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
  subroutine ts_iohs(task, gamma, onlyS, no_u, no_s, Enspin,   &
       indxuo, maxnh, numh, listhptr, listh, H, S, qtot, &
       temp, xij, fnlength, fname, na_u, lasto, isa, ef, &
       ucell, ts_kscell_file, ts_kdispl_file,        	   &
       ts_gamma_scf_file, xa, istep, ia1, check_kcell)

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
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

!
!  Modules
!
    use precision,    only : dp
    use parallel,     only : Node, Nodes
#ifndef TBTRANS
    use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb, &
         GlobalToLocalOrb, GetNodeOrbs
#endif

    use files,        only : slabel, label_length
    use sys,          only : die
#ifdef TBTRANS
    use m_tbt_kpoints, only : ts_kscell => kscell
    use m_tbt_kpoints, only : ts_kdispl => kdispl
    use m_tbt_kpoints, only : ts_gamma_scf => Gamma
#else
    use m_ts_kpoints, only : ts_gamma_scf, ts_kscell, ts_kdispl
#endif

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
    external          io_assign, io_close
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
    logical,optional :: check_kcell
! TSS End

! Internal variables and arrays
    integer    is, iu, k
    integer    ih,hl, maxnhtot, maxhg
    integer, dimension(:), allocatable :: numhg
#ifdef MPI
    integer    MPIerror, Request, Status(MPI_Status_Size), BNode
    integer,  dimension(:),   allocatable :: ibuffer
    real(dp), dimension(:),   allocatable :: buffer
    real(dp), dimension(:,:), allocatable :: buffer2
#endif
    logical    found, gammaonfile, lcheck_kcell
! TSS Begin
! Aux Variables
    integer :: i, j
! TSS End

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io' )
#endif

    lcheck_kcell = .true.
    if ( present(check_kcell) ) then
       lcheck_kcell = check_kcell
    end if

! #################### READ ######################################
! Choose between read or write
    if (task.eq.'read' .or. task.eq.'READ') then

! Note that only the master node is suppossed to
! call this routine in "read" mode...

! Check if input file exists
       inquire( file=fname, exist=found)

       if (found) then

! Open file
          call io_assign( iu )
          open( iu, file=fname, form='unformatted', status='old' )      

! Read dimensions
          read(iu) na_u, no_u, no_s, Enspin, maxnh

! Check dimensions
          if (Enspin  .ne. nspin) then
             if (node == 0) then
                write(6,"(a,i10,a,i10)") "nspin mismatch: In "   &
                     // trim(fname) // " : ",                    &
                     Enspin, ". Current run: ", nspin
             endif
             call die('ts_iohs: TSHS file contains a wrong nspin')
          endif

! TSS Begin
! Allocate arrays that are going to be read now
          nullify(xa,isa)
          allocate(xa(3,na_u))
          allocate(isa(na_u)) 
          call memory('A','D',3*na_u,'iohs')
          call memory('A','I',na_u,'iohs')
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
             write(*,*) 'SIESTA system Gamma:',gamma
             write(*,*) 'SIESTA file Gamma:',gammaonfile
             call die('ts_iohs: Incompatible SIESTA k-sampling')
          endif
          if (.not.ts_gamma_scf .and. ts_gamma_scf_file) then
             write(*,*) 'TranSIESTA file Gamma:',ts_gamma_scf_file
             write(*,*) 'TranSIESTA system Gamma:',ts_gamma_scf
             call die('ts_iohs: Incompatible k-samplings!') 
          end if
          if ( lcheck_kcell ) then
             do i = 1,3
                if ( ts_kdispl_file(i) /= ts_kdispl(i) ) then
                   write(*,'(a,2F8.4)') 'TranSIESTA file k-displacements:', &
                        (ts_kdispl_file(j),j=1,2)
                   write(*,'(a,2F8.4)') 'TranSIESTA system k-displacements:', &
                        (ts_kdispl(j),j=1,2)
                   call die('ts_iohs: Incompatible k-displacements') 
                end if
             enddo
             if ( sum(ts_kscell_file - ts_kscell) == 0 ) then
                write(*,'(a)') 'TranSIESTA file k-grid:'
                do j = 1 , 3
                   write(*,'(t2,3f8.4)') (ts_kscell_file(i,j),i=1,3)
                end do
                write(*,'(a)') 'TranSIESTA system k-grid:'
                do j = 1 , 3
                   write(*,'(t2,3f8.4)') (ts_kscell(i,j),i=1,3)
                end do
                call die('ts_iohs: Incompatible k-grids') 
             end if
          end if

! Read sparse listings
          nullify(lasto)
          allocate(lasto(0:na_u))
          call memory('A','I',1+na_u,'iohs')
          read(iu) lasto(0:na_u)


          if (.not.gammaonfile) then
! Allocate arrays that are going to be read now
             nullify(indxuo)
             allocate(indxuo(1:no_s))
             call memory('A','I',no_s,'iohs')
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
          call memory('A','I',no_u,'iohs')
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
          call memory('A','I',maxnh,'iohs')
          do ih = 1,no_u
             read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih))
          enddo

! Read Overlap matrix
! Allocate S
          nullify(S)
          allocate(S(maxnh))
          call memory('A','D',maxnh,'iohs')
          do ih = 1,no_u
             read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
          enddo

! Read Hamiltonian
! Allocate H
          nullify(H)
          allocate(H(maxnh,Enspin))
          call memory('A','D',maxnh*Enspin,'iohs')
          do is = 1,Enspin
             do ih = 1,no_u
                read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
             enddo
          enddo

          if (.not.gammaonfile) then
! Read interorbital vectors for K point phasing
! Allocate xij
             nullify(xij)
             allocate(xij(3,maxnh))
             call memory('A','D',3*maxnh,'iohs')
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
#ifndef TBTRANS
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
#endif
    else
       call die("ts_iohs have not been given a correct task: "&
            //trim(task)//".")
    endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io' )
#endif

  end subroutine ts_iohs

end module m_ts_io
