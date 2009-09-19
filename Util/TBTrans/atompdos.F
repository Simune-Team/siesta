C ##################################################################
C     ##         PDOS curves resolved into atoms          ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


      subroutine atompdos(ju,nc,nspin,ispin,is1,is2,nua,NGL,
     . NBUFL,lasto,
     .     energy,GF,GFR,S,slabel)

      implicit none
    
c     INPUT

      integer ju
      integer nc,is1,is2,nua,NGL,NBUFL
      integer lasto(0:nua)
      real*8 energy
      complex*16 GF(nc,nc),GFR(nc,nc),S(nc,nc)
      integer ispin,nspin

      character*20 slabel       ! System Label (to name output files)

c     LOCAL

      real*8 pdosTot,pdosL,pdosR
C FDN
      real*8 pdosAtTot    
C FDN
      real*8 pdosL2L,pdosR2L,pdosL2R,pdosR2R,pdos2L,pdos2R

      real*8 Pi
      parameter(Pi=3.141592653589793238462643383279503d0)
      real*8 eV
      parameter(eV=13.6058d0)

      complex*16 ztmp
      real*8 pdos(100)
c     real*8 pdosR(100),pdosL(100)
      integer iu,iu1
      integer i,ia1,ia2,ns1,ns2,il1,il2,i1,i2,io1,io2
      integer iastart,iaend,ia
      logical firsttime

      complex*16, dimension (:,:), allocatable:: GFS, GFRS

      data firsttime /.true./  
      save firsttime,iu,iu1,iastart,iaend
      character paste*33,filename*33
      external paste

c BEGIN

c     Open file

      if(is1.gt.is2) then
         i=is1
         is1=is2
         is2=i
      end if

      if(firsttime) then      
         call io_assign( iu )
         filename = paste(slabel,'.TOTDOS')
         if(nspin.gt.1) then
            if(ispin.eq.1) filename = paste(slabel,'.UP.TOTDOS')
            if(ispin.eq.1) filename = paste(slabel,'.DOWN.TOTDOS')
         end if
         open( iu, file=filename, status='unknown' )   
         write(iu,*) "#Total Mulliken pop. on atoms "
         write(iu,*) "#atom number,energy,Total, Left, Right "

         call io_assign( iu1 )
         filename = paste(slabel,'.ORBDOS') 
         if(nspin.gt.1) then
            if(ispin.eq.1) filename = paste(slabel,'.UP.ORBDOS')
            if(ispin.eq.1) filename = paste(slabel,'.DOWN.ORBDOS')
         end if

         open( iu1, file=filename, status='unknown' )   
         write(iu1,*) "#Mulliken population on orbitals "
         write(iu1,*) "#(orbs in same order as in SIESTA Mulliken)"
         write(iu1,*) "#atom number, energy, Mulliken pop. on orbital" 

         
         
         io1=is1 + NBUFL + NGL
         io2=is2 + NBUFL + NGL
         do ia=1,nua
            if(io1.le.lasto(ia) .and. io1.gt.lasto(ia-1)) iastart=ia !orbital is1 is on atom ia1
            if(io2.le.lasto(ia) .and. io2.gt.lasto(ia-1)) iaend=ia !orbital is1 is on atom ia1
         end do                 !ia

         write(ju,*) 'ATOMIC PDOS: starting on atom: ',iastart,
     .               ' ending on atom: ',iaend 
         firsttime=.false.
      end if

 12   format(i4,1x,F10.4,20(f8.4,1x))

      ztmp = dcmplx(-2d0/Pi,0d0)/eV


      allocate(GFS(nc,nc))
      allocate(GFRS(nc,nc))

      CALL zgemm('N','N',nc,nc,nc,dcmplx(1d0,0d0),
     &     GF,nc,S,nc,dcmplx(0d0,0d0),GFS,nc) 

      CALL zgemm('N','N',nc,nc,nc,dcmplx(1d0,0d0),
     &     GFR,nc,S,nc,dcmplx(0d0,0d0),GFRS,nc) 


      do ia1=iastart,iaend

         pdosL = 0d0
         pdosR = 0d0
         pdosTot  = 0d0
         
c     ---------------------------------            

         ns1 = lasto(ia1) - lasto(ia1-1)            

C FDN
         pdosAtTot=0.d0
C FDN
         do il1=1,ns1
            io1 = il1 + lasto(ia1-1)
            i1  = io1 - NBUFL - NGL
            
            pdos(il1) = DIMAG(ztmp*GFS(i1,i1))
C FDN
            pdosAtTot=pdosAtTot+pdos(il1)
C FDN 

            pdosTot   = pdosTot + DIMAG(ztmp*GFS(i1,i1))
            pdosR     = pdosR   + DIMAG(ztmp*GFRS(i1,i1))


         end do                 !il1


         pdosL = pdosTot  - pdosR            

C FDN         
         write(iu1,'(I4,1X,F9.5,1X,F12.8,1X,30(1X,e16.10))')
     .        ia1,energy*eV,pdosAtTot,(pdos(il1),il1=1,ns1)
C FDN     
    
c     write(6,'(I4,1X,F9.5,1X,30(1X,e16.10))')
c     .        ia1,energy*eV,(pdos(il1),il1=1,ns1)
         
         write(iu,'(I4,1X,F9.5,1X,3(1X,e16.10))')
     .        ia1,energy*eV,pdosTot,pdosL,pdosR
         
            
      end do                    !ia1


      deallocate(GFS)
      deallocate(GFRS)

c     call io_close( iu )

      return
      end





