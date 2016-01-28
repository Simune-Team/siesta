! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_ts_scattering
!
! Routines that are used for the scattering system calculation
!
!=============================================================================
! CONTAINS:
!          1) getGF
!          2) getSFE
!          3) GFGammaGF
!          4) weightDM
!          5) weightDMC
!          6) get_scat_region
      implicit none

      public :: getSFE, getGF, GFGammaGF
      public :: weightDM, weightDMC
      public :: get_scat_region


      private

      CONTAINS

!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


C ##################################################################
C ## Calculating Full Greens functions of                         ## 
C ## contact coupled to electrodes using block tri-diagonal       ##
C ##                                                              ##          
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ##  Modified by F. D. Novaes                                    ##
C ##################################################################

      subroutine getGF(USEBULK,updatedmcr,VoltFDF,no,NGL,NGR,Hi,Si,
     &     HLL,HRR,ZEnergy,GF,ierr,useTri)

      use sys, only: die
      use m_ts_aux_rout, only : splittri, collectri, invmat, transpose,
     &                          csolveg, ts_alloc, ts_dealloc
      use precision, only: dp

      implicit none 


C ======================================================================
c     Note: Indices are starting at 0 in this routine
C ======================================================================
C INPUT:
      logical USEBULK           ! if true self-energy only is input else
c                                 z*S-H-Sigma for bulk is in sfe
      logical updatedmcr
      real(dp) VoltFDF
      integer no                !no : no. states
      complex(dp) Hi(0:(no*no)-1),Si(0:(no*no)-1)
      complex(dp) ZEnergy 
      integer NGL,NGR
      complex(dp) HLL(NGL,NGL)   ! Left electrode GF
      complex(dp) HRR(NGR,NGR)   ! Right electrode GF
C FDN 
C Use the tridiagonality ?
      logical :: useTri
C FDN 
      

C ======================================================================
C INPUT/OUTPUT:


C ======================================================================
C OUTPUT:
C      complex*16 GF(0:(no*no)-1)  ! 1/(Sc*Energy - H)
      complex(dp), intent(out) :: GF(no*no)
      integer ierr              !inversion err
C ======================================================================
c Helpers, tempos ...
      complex(dp), dimension(:), allocatable :: H
      integer, dimension (:),    allocatable :: ipvt

      integer i,j,ii,jj
      complex(dp) a,b
C --------------------------------------------------------------------

      integer n1,n2,n3
C FDN 
      integer :: allocStat
      complex(dp), dimension(:,:), allocatable ::  
     &                                             a11,a12,
     &                                             a21,a22,a23,
     &                                                 a32,a33
      complex(dp), dimension(:,:), allocatable ::
     &                                             b11,b12,
     &                                             b21,b22,b23,
     &                                                 b32,b33
      complex(dp), dimension(:,:), allocatable ::      x12
      complex(dp), dimension(:,:), allocatable ::  x21,    x23
      complex(dp), dimension(:,:), allocatable ::      x32
      logical :: contactGFonly   

c----------------------------------------------------------------------
      
      external timer, memory

      
    
c=======================================================================
C BEGIN
c=======================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE getGF' )
#endif

      call timer('GFT',1) 
      
c ================================================================
      ierr=0
! FDN
! Compute only GF(2,2) ?
! A difference could be made when updatedmcr(.false.) and volt=0.0 ...
! b11 and b33 are needed only when there is a bias, and because of the
! integration scheme, where both L and R electrodes are used as references.
      contactGFonly = updatedmcr .and. (voltfdf == 0.d0)      

c Setup full H (this part will be changed later..)

      call ts_alloc(H,no*no,"getGF")
 
c
c     Contact z*S-H without self-energies:
c
      if(useTri) then
         do j = 0, no-1
            do i = 0, no-1
               H(i + 1 + j*no) = 
     &              ZEnergy*Si(i+j*no)-Hi(i+j*no) 
            end do              !i
         end do                 !j
      else
         do j = 0, no-1
            do i = 0, no-1
               H(i + 1 + j*no) = 
     &              ZEnergy*Si(i+j*no)-Hi(i+j*no) 
               GF(i + 1 + no*j) = dcmplx(0d0,0d0)
            end do              !i
            GF(j + 1 + no*j) = dcmplx(1d0,0d0)
         end do                 !j
      end if 
      
c
c     Add Left Self-energy
c
      if (UseBulk) then 
         do j = 1, NGL
            do i = 1, NGL
               H(i + (j-1)*no) = HLL(i,j)
            end do              !i
         end do                 !j
      else
         do j = 1, NGL
            do i = 1, NGL
               H(i + (j-1)*no) = H(i + (j-1)*no) - HLL(i,j)
            end do              !i
         end do                 !j
      end if                    !USEBULK
c     
c     Add Right Self-energy
c     
      if (UseBulk) then
         do j = no-NGR+1,no
            jj=j-(no-NGR)
            do i = no-NGR+1,no
               ii=i-(no-NGR)
               H(i + (j-1)*no) = HRR(ii,jj)
            end do              !i
         end do                 !j
      else
         do j = no-NGR+1,no
            jj=j-(no-NGR)
            do i = no-NGR+1,no
               ii=i-(no-NGR)
               H(i + (j-1)*no) = H(i + (j-1)*no) - HRR(ii,jj)
            end do              !i
         end do                 !j
      end if
      
      if ( .not. useTri) then
! Invert directly
         allocate(ipvt(no),stat=allocStat)
         if(allocStat /=0) call die("Could not allocate in getGF") 
         call memory('A','I',size(ipvt),"getGF")      

         call zgesv(no,no,H,no,ipvt,GF,no,ierr)            

         call ts_dealloc(H,"getGF")
 
         call memory('D','I',size(ipvt),"getGF")
         deallocate(ipvt,stat=allocStat)
         if(allocStat /=0) call die("Could not deallocate in getGF")

      else

! Use the fact that H is tridiagonal       
c ================================================================
  
c allocate:
      n1 = NGL
      n2 = no - (NGL+NGR)
      n3 = NGR
              
      call ts_alloc(a11,n1,n1,"getGF")
      call ts_alloc(a12,n1,n2,"getGF")
      call ts_alloc(a21,n2,n1,"getGF")
      call ts_alloc(a22,n2,n2,"getGF")
      call ts_alloc(a23,n2,n3,"getGF")
      call ts_alloc(a32,n3,n2,"getGF")
      call ts_alloc(a33,n3,n3,"getGF")
      
      call splittri(n1,n2,n3,H,a11,a12,a21,a22,a23,a32,a33)
      call memory('D','Z',size(H),"getGF")
      deallocate(H,stat=allocStat)
      if(allocStat /=0) call die("Could not deallocate in getGF")
      
      

     
c====================================================================
       
      

c     inv(a11) -> store in b11
      call ts_alloc(b11,n1,n1,"getGF")
      call invmat(n1,a11,b11,ierr)
      call ts_dealloc(a11,"getGF")
      

c     inv(a33) -> store in b33
      call ts_alloc(b33,n3,n3,"getGF")
      call invmat(n3,a33,b33,ierr)
      call ts_dealloc(a33,"getGF")
      

c     x12 = inv(a11).a12 
      call ts_alloc(x12,n1,n2,"getGF")
      a=dcmplx(1.d0,0.d0)
      b=dcmplx(0.d0,0.d0) 
      call zgemm('N','N', n1, n2, n1 ,a ,b11, n1, a12,
     &     n1, b, x12, n1)     
      call ts_dealloc(a12,"getGF")
      if(contactGFonly) call ts_dealloc(b11,"getGF")
      

c     a22 - a21.x12 --> a22
      a=dcmplx(-1.d0,0.d0)
      b=dcmplx(1.d0,0.d0) 
      call zgemm('N','N', n2, n2, n1 ,a ,a21, n2, x12,
     &  n1, b, a22, n2)      
      if(contactGFonly) then
         call ts_dealloc(x12,"getGF")
         call ts_dealloc(a21,"getGF")
      end if
      
c     x32 = inv(a33).a32
      call ts_alloc(x32,n3,n2,"getGF")
      a=dcmplx(1.d0,0.d0)
      b=dcmplx(0.d0,0.d0) 
      call zgemm('N','N', n3, n2, n3 ,a ,b33, n3, a32,
     &     n3, b, x32, n3)   
      call ts_dealloc(a32,"getGF")
      if(contactGFonly) call ts_dealloc(b33,"getGF")
       

c     a22 - a23.x32 --> a22
      a=dcmplx(-1.d0,0.d0)
      b=dcmplx(1.d0,0.d0) 
      call zgemm('N','N', n2, n2, n3 ,a ,a23, n2, x32,
     &  n3, b, a22, n2)
      if(contactGFonly) then
         call ts_dealloc(x32,"getGF")
         call ts_dealloc(a23,"getGF")
      end if      
      
c now a22 contains (a22 - a21.x12 - a23.x32) 

c b22 = inv(a22 - a21.x12 - a23.x32) 
      call ts_alloc(b22,n2,n2,"getGF")
      call invmat(n2,a22,b22,ierr)   ! now we've got b22
      call ts_dealloc(a22,"getGF")

      
      if(.not.contactGFonly) then
         
c  b12 = -x12.b22, b21 = b12^T:
         call ts_alloc(b12,n1,n2,"getGF")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n1 , n2 , n2 ,a ,x12, n1, b22,
     &     n2, b, b12, n1)      ! got b12
         call ts_dealloc(x12,"getGF")
         
C FDN
C      call transpose(b12,n1,n2,b21) ! got b21
C b21 = -b22.a21.inv(a11)
C  a21.inv(a11) --> x21
         call ts_alloc(x21,n2,n1,"getGF")
         a=dcmplx(1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2, n1, n1 ,a ,a21, n2, b11, 
     &     n1, b, x21, n2)
         call ts_dealloc(a21,"getGF")
         
! b21 = -b22.x21
         call ts_alloc(b21,n2,n1,"getGF")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2 , n1 , n2 ,a ,b22, n2, x21, 
     &     n2, b, b21, n2)      ! got b21


! b11 is needed only when there is a bias.
! b11 = inv(a11) - b12.x21
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(1.d0,0.d0) 
         call zgemm('N','N', n1 , n1 , n2 ,a ,b12, n1, x21, 
     &     n2, b, b11, n1)
         call ts_dealloc(x21,"getGF")
C FDN 
 

c  b32 = -x32.b22, b23 = b32^T:
         call ts_alloc(b32,n3,n2,"getGF")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n3 , n2 , n2 ,a ,x32, n3, b22,
     &     n2, b, b32, n3)      ! got b32  
         call ts_dealloc(x32,"getGF")
     
C FDN
C      call transpose(b32,n3,n2,b23) ! got b23
C b23 = -b22.a23.inv(a33)
C   a23.inv(a33) --> x23
         call ts_alloc(x23,n2,n3,"getGF")
         a=dcmplx(1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2, n3, n3 ,a ,a23, n2, b33, 
     &     n3, b, x23, n2)
         call ts_dealloc(a23,"getGF")
         
C b23 = -b22.x23
         call ts_alloc(b23,n2,n3,"getGF")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2 , n3 , n2 ,a ,b22, n2, x23, 
     &     n2, b, b23, n2)      ! got b23


! b33 is needed only when there is a bias.
! b33 = inv(a33) - b32.x23
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(1.d0,0.d0) 
         call zgemm('N','N', n3 , n3 , n2 ,a ,b32, n3, x23, 
     &     n2, b, b33, n3)
         call ts_dealloc(x23,"getGF")

      end if ! .not.contactGFonly
         

c
c     now collect into GF
c     

 
      call collectri(n1,n2,n3,b11,b12,b21,b22,b23,b32,b33,GF)
 
      call ts_dealloc(b22,"getGF")
      
      if(.not.contactGFonly) then
         call ts_dealloc(b11,"getGF")
         call ts_dealloc(b12,"getGF")
         call ts_dealloc(b21,"getGF")
         call ts_dealloc(b23,"getGF")
         call ts_dealloc(b32,"getGF")
         call ts_dealloc(b33,"getGF")
      end if

      end if ! .not. useTri
      call timer('GFT',2)  

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS getGF' )
#endif
     
C ====================================================================
      END subroutine getGF
C ====================================================================



!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------



C ##################################################################
C ## Subroutine which read-in GF for the 1x1 surface cell         ##
C ## and returns the GF for the NA1xNA2 surface cell              ##
C ## by summing over the appropriate q-points                     ## 
C ##                                                              ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Changed to F90 by Jose-Luis Mozos, jlm@kanigo.icmab.es       ##
C ## Reimplemented expansion of the electrode by:                 ##
C ##    Nick Papior Andersen, nickpapior@gmail.com                ##
C ##################################################################

      subroutine getSFE(UseBulk,uGF,HAA,SAA,ZEnergy,ikpt,nq,qb,wq,
     &     no_GF,nua,lasto,NA1,NA2,no,SFE,GAMMA,points2read,
     &     errorgs)

      use m_ts_aux_rout, only : csolve 
      use precision, only : dp
      use units,     only : Pi
      use parallel,  only : IONode, Node, Nodes

#ifdef MPI
      use mpi_siesta, only : MPI_Bcast, MPI_AllReduce, 
     & MPI_Sum, MPI_Integer, MPI_Status_Size, 
     & MPI_Comm_World
      use mpi_siesta, only : DAT_dcomplex => MPI_double_complex
               
#endif

c-----------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------
      real(dp), parameter :: EPS = 1d-7
 
c ==================================================================
c INPUT

      logical frstme
      integer icall
      data      frstme /.true./
      save frstme,icall

#ifdef MPI
      integer  MPIerror, Request, Status(MPI_Status_Size)

#endif

      logical USEBULK           ! if true return self-energy only else
c                                 return z*S-H-Sigma for bulk

      integer uGF              !fileunit from which we read GF

      integer no                ! size of GF in NA1xNA2 supercell

      integer,    intent(in) :: nua ! no. uc (1x1) atoms
      integer,    intent(in) :: lasto(0:nua) ! index of last orbitals of uc atoms.
      integer,    intent(in) :: NA1,NA2 ! repetitions of surface UC.
      integer,    intent(in) :: nq !no. q-points <= NA1*NA2
      real(dp),   intent(in) :: qb(3,nq),wq(nq) !q-points (in reciprocal space) and weights 
      complex(dp) cphase,ccphase ! for use in transformation in Unit-cell
!     loop counters for Unit-cell conversion
      integer                :: ia1,ia2,iau,iuo 
      integer                :: ja1,ja2,jau,juo

      complex(dp) ZEnergy 
      integer, intent(in) :: ikpt

c ----------- 
c 1x1 surface cell
c
      integer no_GF               ! size of used part of GAA in UC.

      complex(dp) HAA(no_GF,no_GF,nq),SAA(no_GF,no_GF,nq)
      integer points2read

      integer ikpar
c======================================================================
c OUTPUT: 
      complex(dp) SFE(no,no) !Electrode atoms self-energy for USEBULK true
c                           !else z*S - H - Sigma for ideal surface
      complex(dp)  GAMMA(no,no)  ! 0.5*(Sigma - Sigma^dagger) 
      logical errorgs

c=====================================================================
c Helpers, tempos etc...
      integer ierr,itmp              ! inversion error      
      complex(dp) ZEi, wgfi
      integer iEni
      integer iq,inode
      integer io
      integer jo
      complex(dp) zenloc

      complex(dp), allocatable :: WRK1(:,:),WRK2(:,:),GAA(:,:,:)
      integer,dimension(:), allocatable :: ipvt
#ifdef MPI
      complex(dp), allocatable :: cbuffer(:,:,:)
#endif
      integer  ntot

      external memory

c=======================================================================
c     BEGIN:
c=======================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE getSFE' )
#endif

      call timer('getSFE',1)
      
      if(frstme) icall=1
      if(.not.frstme) icall=icall+1
      frstme=.false.

      errorgs =.false.
      
      ntot =  no_GF*no_GF*nq
      
      if(points2read .gt. Nodes) then
         if(IONode)
     &    write(6,*) 
     &    'ERROR GetSFE: points2read=',points2read,' > Nodes=',Nodes
          errorgs = .true.
          return
      end if
      if(points2read .le. 0) then
         if(IONode)
     &   write(6,*) 'ERROR GetSFE: points2read=',points2read,'< 0'  
         errorgs = .true.
         return
      end if

c
c     Init:
c
      allocate(WRK1(no,no))
      allocate(WRK2(no,no))
      allocate(ipvt(no))
      call memory('A','Z',2*no*no,'getsfe')
      call memory('A','I',no,'getsfe')
      allocate(GAA(no_GF,no_GF,nq))
      call memory('A','Z',no_GF*no_GF*nq,'getsfe')
      

#ifdef MPI
c In order for the inversion to work when inode > points2read
      GAA=dcmplx(0.d0,0.d0)
      do iq=1,nq
         do io=1,no_GF
            GAA(io,io,iq)=dcmplx(1.d0,1.d0)
         end do
      end do
#endif

  
c Loop over nodes. If root send stuff out to relevant nodes, if
c not root receive information from root.
c If there is more nodes than points left to read in file (=points2read)
c then just send out the data of the last point to the remaining nodes.


#ifdef MPI
      allocate(cbuffer(no_GF,no_GF,nq))
      call memory('A','Z',no_GF*no_GF*nq,'getsfe')
#endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      do inode = 0, Nodes-1
cccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ------ root read -------

         if(IONode .and. inode .lt. points2read) then
            read(uGF) iEni,ZEi,wgfi,ikpar !still enough points to feed to processors
         endif

#ifdef MPI
         if(Node.eq.inode) then
            zenloc = zenergy
         endif
        
         call MPI_Bcast(zenloc,1,DAT_dcomplex,inode,
     &        MPI_Comm_World,MPIerror)  
         
         call MPI_Bcast(iEni,1,MPI_integer,0,
     &        MPI_Comm_World,MPIerror)  
#else
         zenloc = zenergy
#endif

         if(IONode) then
            if(cdabs(ZEi-ZEnloc) .gt. EPS) then
               write(6,*)
     &              'ERROR: getSFE: ZEnergy=',ZEi,' expected:',ZEnloc,
     &              ' in processor ', inode
               errorgs = .true.
               return
            end if
C NPA Removed
            if(.not. (ikpt.eq.ikpar)) then
               write(6,*)
     &              'ERROR: getSFE: ikpar=',ikpar,' expected:',ikpt
               errorgs = .true.
               return
            end if
         end if
         
         if(iEni .eq. 1) then
            if(IONode) then
               read(uGF) HAA
               read(uGF) SAA
            endif
           
                     
#ifdef MPI
            call MPI_Bcast(HAA(1,1,1),no_GF*no_GF*nq,DAT_dcomplex,
     &           0,MPI_Comm_World,MPIerror)  
            call MPI_Bcast(SAA(1,1,1),no_GF*no_GF*nq,DAT_dcomplex,
     &           0,MPI_Comm_World,MPIerror)           
#endif
         end if                 !iEni==1
       


         if(inode .lt. points2read)  then

#ifdef MPI
            if(inode.eq.0.and.IONode)then
#endif

               read(uGF) GAA   !still enough points to feed to processors

#ifdef MPI
            elseif (IONode) then
               
               read(uGF) cbuffer
               
               call MPI_ISend(cbuffer(1,1,1),no_GF*no_GF*nq,DAT_dcomplex
     &              ,inode,1,MPI_Comm_World,Request,MPIerror) 
               call MPI_Wait(Request,Status,MPIerror)
               
            elseif (Node.eq.inode) then
               call MPI_IRecv(GAA(1,1,1),no_GF*no_GF*nq,DAT_dcomplex,
     &              0,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
               
            endif               !inode=0 & Node=0

#endif
         end if                 ! inode.lt.points2read                    
c     -----------------------

c     -----------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do                    !loop over nodes to distribute GF's etc ..
cccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
      call memory('D','Z',size(cbuffer),'getsfe')
      deallocate(cbuffer)
#endif

!-----------------------------------------------------------------------
!     Expansion of unit-cell: (1x1) -> (NA1 x NA2) Layer-wise ordering!
!-----------------------------------------------------------------------
      do io=1,no
         do jo=1,no
            WRK1(jo,io) = dcmplx(0d0,0d0)
            WRK2(jo,io) = dcmplx(0d0,0d0)
         end do                 ! jo
      end do                    ! io

      do iq=1,nq                !q-points
Cc----------------------------------------------
        io=0
        do iau=1,nua
         do ia2=0,NA2-1
          do ia1=0,NA1-1
           do iuo=1+lasto(iau-1),lasto(iau)  
            io=io+1
            jo=0
            do jau=1,nua
             do ja2=0,NA2-1
              do ja1=0,NA1-1
               do juo=1+lasto(jau-1),lasto(jau)  
                jo=jo+1
                cphase = cdexp(dcmplx(0.0_dp,-2.0_dp*Pi) *
     &               ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
!                ccphase = cdexp(dcmplx(0.0_dp,2.0_dp*Pi) *
!     &               ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
                
                WRK1(jo,io) = WRK1(jo,io) + dcmplx(wq(iq),0.0_dp)*
     &               cphase*GAA(juo,iuo,iq)

!     Here we employ the time-reversal symmetry: 
!          GF(1x1)(q)(i,j) = GF(1x1)(-q)(j,i)  
!                WRK1(jo,io) = WRK1(jo,io) + dcmplx(0.5_dp*wq(iq),0.0_dp)*
!     &               (cphase*GAA(juo,iuo,iq)+ccphase*GAA(iuo,juo,iq))
                
                WRK2(jo,io) = WRK2(jo,io) + dcmplx(wq(iq),0.0_dp)*
     &               cphase*(ZEnergy*SAA(juo,iuo,iq)-HAA(juo,iuo,iq))

!     Here we employ the time-reversal symmetry: 
!          GF(1x1)(q)(i,j) = GF(1x1)(-q)(j,i)  
!                WRK2(jo,io) = WRK2(jo,io) + dcmplx(0.5_dp*wq(iq),0.0_dp)*
!     &               (cphase*(ZEnergy*SAA(juo,iuo,iq)-HAA(juo,iuo,iq))
!     &               ccphase*(ZEnergy*SAA(iuo,juo,iq)-HAA(iuo,juo,iq)))
                
!     For debugging purposes...
!                ccphase  =  exp(dcmplx(0d0,2d0*Pi)  *  
!     &               ( (ia1-ja1)*q(1,iq) + (ia2-ja2)*q(2,iq)) )    

!!               WRK2(jo,io) = WRK2(jo,io) + 
!!     &              dcmplx(0.5d0*wq(iq),0d0)*
!!     &              (cphase*SAA(juo,iuo,iq)+ccphase*SAA(iuo,juo,iq))

               enddo !juo
              enddo !ja1
             enddo !ja2
            enddo !jau
           enddo !iuo
          enddo !ia1
         enddo !ia2
        enddo !iau
c----------------------------------------------
      enddo                     !q-points
Cc----------------------------------------------
      do io=1,no
         do jo=1,no
            SFE(jo,io) = dcmplx(0.d0,0.d0)
         end do                 ! jo
         SFE(io,io) = dcmplx(1.d0,0.d0) !to be used in inv.
      end do                    ! io

         
c     Find inverse of GAA --> SFE
      ierr=0
!the csolve has arrays of NG,NG in the interface...
      call csolve(no,WRK1,SFE,ipvt,ierr)

      if(ierr.ne.0) write(90+Node,*) 
     &     "err hola: ,Node,icall,ierr,zenloc,NG" 
      if(ierr.ne.0) write(90+Node,*) 
     &     "err hola: ",Node,icall,ierr,zenloc,no 

#ifdef MPI
      itmp=0
      call MPI_AllReduce(ierr,itmp,1,MPI_integer,MPI_sum,
     &     MPI_Comm_World,MPIerror)
      ierr=itmp      
#endif

      IF(IERR.NE.0) THEN
c         errorgs=.true.
         if(Node.eq.0) then
            write(6,*) 'ERROR: getsfe 1 MATRIX INVERSION FAILED'
            write(6,*) 'ERROR: ierr=',ierr
         endif 
      END IF


      do io=1,no
         do jo=1,no
c     SFE = Z*S - HAA - Sigma ---> WRK2 - SFE = Sigma -> WRK1:
            WRK1(jo,io) = WRK2(jo,io)-SFE(jo,io)
         end do                 !jo
      end do                    !io

c
c Gamma 
c
      do io=1,no
         do jo=1,no
            GAMMA(jo,io) = 
     &           dcmplx(0d0,0.5d0)*(WRK1(jo,io) - 
     &           dconjg(WRK1(io,jo)))
         end do                 !jo
      end do                    !io
      
      if(.not.USEBULK) then
c      
c    Sigma=(Z*S - H) - inv[GF] ---> SFE
c
         do io=1,no
            do jo=1,no
               SFE(jo,io) = WRK1(jo,io)
            end do              !jo
         end do                 !io

      end if                    !USEBULK


c     debug
c      do io=1,no
c         do jo=1,no
c            SFE(jo,io) = WRK2(jo,io)
c         end do                 !jo
c      end do                    !io


       call memory('D','Z',2*size(WRK1),'getsfe')
       call memory('D','Z',size(GAA),'getsfe')
       call memory('D','I',size(ipvt),'getsfe')
       deallocate(WRK1,WRK2,GAA,ipvt)

       call timer('getSFE',2)

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS getSFE' )
#endif

c======================================================================
      end subroutine getSFE
c======================================================================


!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      subroutine GFGammaGF(LEFTFLAG,nc,ng,GF,Gamma,GGG)

C ======================================================================
c  This routine returns GGG=(-i)*GF.Gamma.GF, where GF is a (nc)x(nc)
c  matrix and the first[LEFTFLAG true]/last[LEFTFLAG false] states
c  corresponds to the (ng) Left/Right electrode states.
c  Gamma is a (ng)x(ng) matrix.
C ======================================================================
      use precision, only : dp

      implicit none
C ======================================================================
C Edited by Nick Papior Andersen to limit memory use and do a lot fewer
C operations.
C However, this is at the expentiure of more cachemisses.
C It has however proved to be more than twice as fast.
c
C ======================================================================
C INPUT:

      
      logical, intent(in) :: LEFTFLAG !True if LEFT
      integer, intent(in) :: nc !no. states in contact region
      integer, intent(in) :: ng !if LeftFLAG the ng first states of nc
c                               !are left electrode states
c                               !else it is the ng last right electrode
c                               !states.
      complex(dp), intent(in) :: Gamma(ng,ng)   !i(Sigma-Sigma^dagger)/2
      complex(dp), intent(in) :: GF(nc,nc)

C ======================================================================
C OUTPUT:
C      
      complex(dp), intent(out) :: GGG(nc,nc)    !GF.GAMMA.GF

C ======================================================================
c Helpers, tempos ...
      complex(dp), allocatable :: wrk(:)
      integer :: i,j,ie,je

c=======================================================================
C BEGIN
c=======================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE GFGammaGF' )
#endif

      call timer("GFGGF",1)

      if(LEFTFLAG) then
         
         allocate(wrk(ng))
         call memory('A','Z',ng,'GFGammaGF')

         do i=1 , nc
            do ie=1 , ng
               wrk(ie) = dcmplx(0.0d0,0.0d0)
               do je=1 , ng
                  wrk(ie) = wrk(ie) +
     &                 dreal(Gamma(je,ie)) * GF(i,je)
               enddo
            enddo

            do j=1 , i
               GGG(j,i) = dcmplx(0.d0,0.d0)
               do ie=1 , ng
                  GGG(j,i) = GGG(j,i) +
     &                 dreal( wrk(ie)*dconjg(GF(j,ie)) )
               enddo
               GGG(j,i) = GGG(j,i) * dcmplx(0d0,-1.d0)
               GGG(i,j) = GGG(j,i)
            enddo

         enddo
         
      else ! Do right
 
         allocate(wrk(nc-ng+1:nc))
         call memory('A','Z',ng,'GFGammaGF')
         
         do i=1 , nc
            do ie= nc-ng+1 , nc
               wrk(ie) = dcmplx(0.0d0,0.0d0)
               do je= nc-ng+1 , nc
                  wrk(ie) = wrk(ie) +
     &                 dreal(Gamma(je-nc+ng,ie-nc+ng)) * GF(i,je)
               enddo
            enddo
            
            do j=1,i
               GGG(j,i) = dcmplx(0.d0,0.d0)
               do ie= nc-ng+1,nc
                  GGG(j,i) = GGG(j,i) +
     &                 dreal( wrk(ie)*dconjg(GF(j,ie)) )
                  
               enddo
               GGG(j,i) = GGG(j,i) * dcmplx(0d0,-1.d0)
               GGG(i,j) = GGG(j,i)
            enddo
            
         enddo
         
      endif
 
      call memory('D','Z',ng,'GFGammaGF')
      deallocate(wrk)

      call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS GFGammaGF' )
#endif

C ====================================================================
      END subroutine GFGammaGF
C ====================================================================





C ##################################################################
C ##   Mixing the Density matrixes according to the smallest      ##
C ##    realspace integral                                        ##
C ##                                                              ##
C ##  Version 011200  Kurt Stokbro, ks@mic.dtu.dk                 ##
C ##################################################################


      subroutine weightDM(nc,NGL,NGR,DML,DMR,DMneqL,DMneqR,
     &     EDML,EDMR)
C ======================================================================
c  This routine find weight for the DM integral. On output
C  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
C  EDML:= w EDML +(1-w) EDMR
C  In left electrode w=1 and in right electrode w=0
C ======================================================================

      use sys, only : die
      use precision, only : dp
      use parallel  , only: IONode

      implicit none

      logical PRINTALOT
      parameter(PRINTALOT=.FALSE.)
c      parameter(PRINTALOT=.TRUE.)

c
C ======================================================================
C INPUT:

      
      integer nc                !no. states in density matrix
      integer NGL               !no. states in left electrode
      integer NGR               !no. states in right electrode
      real(dp) DML(nc,nc),DMR(nc,nc)      !contour part of DM integration
      real(dp) DMneqL(nc,nc),DMneqR(nc,nc)!real axis part of DM integration
      real(dp) EDML(nc,nc),EDMR(nc,nc)    !L,R estimates of Edm

C ======================================================================
c Helpers, tempos ...
      real(dp) w,wl,wr,wsum
      real(dp) maxdif2            !largest integration error
      real(dp) maxndif2           !largest normalized integration error
      
      integer i,j


c=======================================================================
c BEGIN:
c=======================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE weightDM' )
#endif

      call timer('weightDM',1)
C Get Node number

      maxdif2=0.0d0
      maxndif2=0.0d0
      do j=1,nc
         do i=1,nc
! The values are real as we are in the Gamma point... No need
!     to complicate it...
            wr=DMneqR(i,j)*DMneqR(i,j)
            wl=DMneqL(i,j)*DMneqL(i,j)
            wsum=wr+wl
            if(wsum .gt. 0d0) then
                w = wl/wsum
             else
                w=0.5d0
             end if

c            write(*,*) i,j,DMneqR(i,j),DMneqL(i,j),DML(i,j),DMR(i,j)
c  In the leads we always use either left or right integration scheme.

            if (i .le. NGL .or. j .le. NGL) then
! FDN With TriDiag, it must be the opposite !! 
! GF(1,3) and GF(3,1) are not zero !, only in TriDiag ...
!               w = 1.0d0
               w = 0.d0
            endif
            if (i .gt. nc-NGR .or. j .gt. nc-NGR) then
! FDN With TriDiag, it must be the opposite !!
!               w = 0.0d0
               w = 1.d0
            endif
c weight DM and DML
            DML(i,j) = w *(DML(i,j)+DMneqR(i,j))
     &           +(1.0_dp-w)*(DMR(i,j)+DMneqL(i,j))
            EDML(i,j) = w*EDML(i,j)+(1.0_dp-w)*EDMR(i,j)

C make an integration error estimate
! FDN
!            dif=DML(i,j)+DMneqR(i,j)-(DMR(i,j)+DMneqL(i,j))
!            dif2=dif*dif
!            ndif2=dif2*wl*wr/(wsum*wsum)
!            if (dif2 .gt. maxdif2) then
!               maxdif2=dif2
!               idif2=i
!               jdif2=j
!            endif
!            if (ndif2 .gt. maxndif2) then
!               maxndif2=ndif2
!               indif2=i
!               jndif2=j
!            endif
         enddo                  !i
      enddo                     !j

! FDN
!      if(Node.eq.0) then
!
!      if (PRINTALOT) then
!         write(joutfile,*)
!     .    'Integration error estimate',indif2,jndif2, sqrt(maxndif2)
!         write(joutfile,*)
!     .    'unnormalized error ',idif2,jdif2, sqrt(maxdif2)
!      endif
!
!      endif

      call timer('weightDM',2)
C ====================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS weightDM' )
#endif

      END subroutine weightDM
C ====================================================================


C FDN Routine for complex matrices
      subroutine weightDMC(nc,NGL,NGR,DML,DMR,DMneqL,DMneqR,
     &     EDML,EDMR)
C ======================================================================
c  This routine find weight for the DM integral. On output
C  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
C  EDML:= w EDML +(1-w) EDMR
C  In left electrode w=1 and in right electrode w=0
C ======================================================================

      use sys, only : die
      use precision, only : dp
      use parallel, only : IONode
      implicit none

      logical PRINTALOT
      parameter(PRINTALOT=.FALSE.)
c      parameter(PRINTALOT=.TRUE.)

c
C ======================================================================
C INPUT:

      
      integer nc                !no. states in density matrix
      integer NGL               !no. states in left electrode
      integer NGR               !no. states in right electrode
      complex(dp) DML(nc,nc),DMR(nc,nc)      !contour part of DM integration
      complex(dp) DMneqL(nc,nc),DMneqR(nc,nc)!real axis part of DM integration
      complex(dp) EDML(nc,nc),EDMR(nc,nc)    !L,R estimates of Edm

C ======================================================================
c Helpers, tempos ...
      real(dp) w,wl,wr,wsum
      real(dp) maxdif2            !largest integration error
      real(dp) maxndif2           !largest normalized integration error
      
      integer i,j

c=======================================================================
c BEGIN:
c=======================================================================

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE weightDMC' )
#endif

      call timer('weightDMC',1)

      maxdif2=0.0d0
      maxndif2=0.0d0
      do j=1,nc
         do i=1,nc
            wr=cdabs(DMneqR(i,j))**2
            wl=cdabs(DMneqL(i,j))**2
            wsum=wr+wl
            if(wsum .gt. 0d0) then
                w = wl/wsum
             else
                w=0.5d0
             end if

c            write(*,*) i,j,DMneqR(i,j),DMneqL(i,j),DML(i,j),DMR(i,j)
c  In the leads we always use either left or right integration scheme.

            if (i .le. NGL .or. j .le. NGL) then
! FDN With TriDiag, it must be the opposite !!
!               w = 1.0d0
               w = 0.d0
            endif
            if (i .gt. nc-NGR .or. j .gt. nc-NGR) then
! FDN With TriDiag, it must be the opposite !!
!               w = 0.0d0
               w = 1.d0
            endif
c weight DM and DML
            DML(i,j) = w   *(DML(i,j)+DMneqR(i,j))
     &           +(1.0_dp-w)*(DMR(i,j)+DMneqL(i,j))
            EDML(i,j) = w*EDML(i,j)+(1.0_dp-w)*EDMR(i,j)
C make an integration error estimate
! FDN
!            dif=cdabs(DML(i,j)+DMneqR(i,j)-(DMR(i,j)+DMneqL(i,j)))
!            dif2=dif*dif
!            ndif2=dif2*wl*wr/(wsum*wsum)
!            if (dif2 .gt. maxdif2) then
!               maxdif2=dif2
!               idif2=i
!               jdif2=j
!            endif
!            if (ndif2 .gt. maxndif2) then
!               maxndif2=ndif2
!               indif2=i
!               jndif2=j
!            endif
         enddo                  !i
      enddo                     !j

! FDN
!      if(Node.eq.0) then

!      if (PRINTALOT) then
!         write(joutfile,*)
!     .    'Integration error estimate',indif2,jndif2, sqrt(maxndif2)
!         write(joutfile,*)
!     .    'unnormalized error ',idif2,jdif2, sqrt(maxdif2)
!      endif
!
!      endif

      call timer('weightDMC',2)

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS weightDMC' )
#endif

C ====================================================================
      END subroutine weightDMC
C ====================================================================


!     Function for calculating which region in the scattering matrix the
!     hamiltonian elements adhere.
      pure function get_scat_region(io,noL,jo,noR,no_u)
      integer, intent(in) :: io, noL, jo, noR, no_u
      integer :: get_scat_region
      
      
      if ( io < 1 .and. jo < 1 ) then
         get_scat_region = 1 ! Left buffer
      else if ( (io < 1 .and. jo <= noL) .or.
     &        (io <= noL .and. jo < 1) ) then
         get_scat_region = 2 ! Left buffer / left electrode
      else if ( (1 <= io .and. 1 <= jo) .and.
     &        (io <= noL .and. jo <= noL) ) then
         get_scat_region = 3 ! Left electrode
      else if ( ((1 <= io .and. noL < jo) .and.
     &        (io <= noL .and. jo <= no_u-noR)) .or.
     &        ((noL < io .and. 1 <= jo) .and.
     &        (io <= no_u-noR .and. jo <= noL)) ) then
         get_scat_region = 4 ! Left electrode / contact region
      else if ( (noL < io .and. noL < jo) .and.
     &        (io <= no_u-noR .and. jo <= no_u-noR) ) then
         get_scat_region = 5 ! Contact region
      else if ( ((noL < io .and. no_u-noR < jo) .and.
     &        (io <= no_u-noR .and. jo <= no_u)) .or.
     &        ((no_u-noR < io .and. noL < jo) .and.
     &        (io <= no_u .and. jo <= no_u-noR)) ) then
         get_scat_region = 6 ! Contact region / right electrode
      else if ( (no_u-noR < io .and. no_u-noR < jo) .and.
     &        (io <= no_u .and. jo <= no_u) ) then
         get_scat_region = 7 ! Right electrode
      else if ( ((no_u-noR < io .and. no_u < jo) .and.
     &        io <= no_u) .or.
     &        ((no_u < io .and. no_u-noR < jo) .and.
     &        jo <= no_u) ) then
         get_scat_region = 8 ! Right electrode / right buffer
      else if ( no_u < io .and. no_u < jo ) then
         get_scat_region = 9 ! Right buffer
      else
         get_scat_region = 0 ! Everything else (LE/RE ... etc.)
      end if
      
      end function get_scat_region

      END MODULE m_ts_scattering
