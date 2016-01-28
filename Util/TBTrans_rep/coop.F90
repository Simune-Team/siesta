! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! #####################################################################
! ##           COOP curves resolved into atoms and orbitals          ##
! ##                              By                                 ##
! ##            Nick Papior Andersen, nickpapior@gmail.com           ##
! #####################################################################
subroutine COOP(uC,uCL,uCR, spin_F, &
     IsoAt1,IsoAt2, &
     noBufL,noL,noD, &
     na_u,lasto,GF,GFRGF,S, &
     iE,Energy)

  use precision,       only : dp
  use sys,             only : die
  use parallel,        only : Node, Nodes, IOnode
  use units,           only : Pi, eV
#ifdef MPI
  use mpi_siesta
#endif
  use m_tbt_out,     only : out_COOP, out_COOPLR

  implicit none

! **************************
! * INPUT variables        *
! **************************
  integer, intent(in) :: uC, uCL, uCR ! The units to write to
  real(dp), intent(in):: spin_F ! The spin factor
  integer, intent(in) :: IsoAt1, IsoAt2 ! Isolated atoms (index in lasto)
  integer, intent(in) :: noBufL, noL, noD ! Orbital counts, left buf, left elec, device
  integer, intent(in) :: na_u
  integer, intent(in) :: lasto(0:na_u) ! Last orbital index of each atom (full lasto)
  complex(dp), intent(in) :: GF   (noD,noD)
  complex(dp), intent(in) :: GFRGF(noD,noD)
  complex(dp), intent(in) :: S    (noD,noD)
  integer, intent(in)     :: iE
  real(dp), intent(in)    :: Energy

! **************************
! * LOCAL variables        *
! **************************
  real(dp), parameter :: r1dPi = 1.0_dp/Pi ! Local Pi factor
  integer             :: noShift ! The shift in orbitals duo to buffer and electrode
  integer             :: ia1, ia2
  integer             :: io1, io2 ! for loops
  real(dp)            :: coopTMP1, coopTMP2
  real(dp)            :: coopT,coopL,coopR          !
  real(dp)            :: coopL2L,coopR2L,coopL2R    !
  real(dp)            :: coopR2R

#ifdef MPI
  real(dp), allocatable :: buf_recv(:), buf_send(:)
  integer :: status(MPI_STATUS_SIZE), req
  integer :: MPIerror, iNode
#endif 

  call timer('coop',1)

!***********************************************************************
!     BEGIN

  noShift = noBufL + noL
#ifdef MPI
  allocate(buf_send(7))
  allocate(buf_recv(7))
#endif
  
  do ia1 = IsoAt1 , IsoAt2
     do ia2 = IsoAt1 , IsoAt2

        !coopL  = 0.0_dp
        coopT  = 0.0_dp
        coopR  = 0.0_dp

        do io2 = lasto(ia2-1) - noShift + 1 , lasto(ia2) - noShift
           do io1 = lasto(ia1-1) - noShift + 1 , lasto(ia1) - noShift
              coopT = coopT - r1dPi*DIMAG(GF   (io1,io2)*S(io1,io2))
              coopR = coopR - r1dPi*DIMAG(GFRGF(io1,io2)*S(io1,io2))
           end do
        end do
        coopT = spin_F * coopT 
        coopR = spin_F * coopR
        
        coopL = coopT - coopR

#ifdef MPI
        buf_send(1) = Energy
        buf_send(2) = coopT
        buf_send(3) = coopL
        buf_send(4) = coopR
        do iNode = 0 , Nodes-1
           if ( IONode ) then
              if ( iNode == Node ) then
                 buf_recv(:) = buf_send(:)
              else
                 call MPI_IRecv(buf_recv,4,MPI_Double_Precision, &
                      iNode,iNode,MPI_Comm_World,req,MPIerror) 
                 call MPI_Wait(req,status,MPIerror)
              end if
              call out_COOP(uC, ia1, ia2, &
                   buf_recv(1),buf_recv(2),buf_recv(3),buf_recv(4))
!              write(uC,12) ikpt, ' AO: ', &
!                   int(buf_recv(1)),int(buf_recv(2), &
!                   buf_recv(3:6)
           else if ( iNode == Node ) then
              call MPI_ISend(buf_send,4,MPI_Double_Precision, &
                   0,iNode,MPI_Comm_World,req,MPIerror) 
              call MPI_Wait(req,status,MPIerror)
           end if
        end do
#else
        call out_COOP(uC,ia1,ia2,Energy,coopT,coopL,coopR)
!        write(iC,12)ikxy,' AO: ',ia1,ia2,Energy/eV,coopT,coopR,coopL
#endif
        
     end do ! ia2
  end do  ! ia1


!=======================================================================
!     summedup COOP to all left/right atoms w.r.t. the given atom
!=======================================================================
  do ia1 = IsoAt1 , IsoAt2

     ! This initialization was done in the orbital loop,
     ! however, that did not make sense, as the calculation
     ! was only performed on the last orbital on each atom.
     !coopL2L = 0.0_dp
     coopL   = 0.0_dp
     !coopL2R = 0.0_dp
     coopR2R = 0.0_dp
     coopR   = 0.0_dp
     coopR2L = 0.0_dp
     

     do ia2 = IsoAt1 , IsoAt2
        do io2 = lasto(ia2-1) - noShift + 1 , lasto(ia2) - noShift

           ! for better performance
           do io1 = lasto(ia1-1) - noShift + 1 , lasto(ia1) - noShift

              ! Note that GFRGF= GF^dagger.(GammaR/2).GF,
              ! also that "GFRGF" + "GFIGF" = Im(G)
              coopTMP1 = - spin_F * r1dPi * DIMAG(GF   (io1,io2)*S(io1,io2))
              coopTMP2 = - spin_F * r1dPi * DIMAG(GFRGF(io1,io2)*S(io1,io2))
              if ( ia2 < ia1 ) then      !atom to the Left
                 coopL   = coopL   + coopTMP1
                 coopR2L = coopR2L + coopTMP2
              else if ( ia2 > ia1 ) then !atom to the Right
                 coopR   = coopR   + coopTMP1
                 coopR2R = coopR2R - coopTMP2
              end if ! atom 2 Right

           end do ! io1 -- orbital number on atom ia1

        end do ! io2 -- orbital number on atom ia2
     end do ! ia2
     
     coopL2L = coopL - coopR2L
     coopL2R = coopR - coopR2R

#ifdef MPI
     buf_send(1) = Energy
     buf_send(2) = coopL2L
     buf_send(3) = coopL
     buf_send(4) = coopL2R
     buf_send(5) = coopR2R
     buf_send(6) = coopR
     buf_send(7) = coopR2L
     do iNode = 0 , Nodes-1
        if ( IONode ) then
           if ( iNode == Node ) then
              buf_recv(:) = buf_send(:)
           else
              call MPI_IRecv(buf_recv,7,MPI_Double_Precision, &
                   iNode,iNode,MPI_Comm_World,req,MPIerror) 
              call MPI_Wait(req,status,MPIerror)
           end if
           call out_COOPLR(uCL, ia1, &
                buf_recv(1),buf_recv(2),buf_recv(3),buf_recv(4))
           call out_COOPLR(uCR, ia1, &
                buf_recv(1),buf_recv(5),buf_recv(6),buf_recv(7))

!           write(iu,11)int(bufN2(1,iNode)),' LR: ',
!           .                  int(bufN2(2,iNode)),bufN2(3:7,iNode)
        else if ( iNode == Node ) then
           call MPI_ISend(buf_send,7,MPI_Double_Precision, &
                0,iNode,MPI_Comm_World,req,MPIerror) 
           call MPI_Wait(req,status,MPIerror)
        end if
     end do
     
#else
     call out_COOPLR(uCL, ia1, &
          Energy,coopL2L,coopL,coopL2R)
     call out_COOPLR(uCR, ia1, &
          Energy,coopR2R,coopR,coopR2L)
!     write(iu,11)ikxy,' LR: ',ia1,Energy/eV,coopL2L,coopL2R,
!     .                                         coopR2L,coopR2R
#endif
  end do ! ia1

  call timer('coop',2)

end subroutine coop
