C ##################################################################
C ##         COOP curves resolved into atoms and orbitals         ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


      subroutine coop(ju,nc,is1,is2,nua,NGL,NBUFL,lasto,
     .     energy,GF,GFR,S)

      implicit none
    
c     INPUT

      integer ju
      integer nc,is1,is2,nua,NGL,NBUFL
      integer lasto(0:nua)
      real*8 energy
      complex*16 GF(nc,nc),GFR(nc,nc),S(nc,nc)


c     LOCAL

      real*8 coopTot,coopL,coopR
      real*8 coopL2L,coopR2L,coopL2R,coopR2R,coop2L,coop2R

      real*8 Pi
      parameter(Pi=3.141592653589793238462643383279503d0)
      real*8 eV
      parameter(eV=13.6058d0)

      complex*16 ztmp
c      real*8 pdos(100),pdosR(100),pdosL(100)
      integer iu,iuL,iuR
      integer i,ia1,ia2,ns1,ns2,il1,il2,i1,i2,io1,io2
      integer iastart,iaend,ia
      logical firsttime
      data firsttime /.true./  
      save firsttime,iu,iuL,iuR,iastart,iaend


c BEGIN

c     Open file

      if(is1.gt.is2) then
         i=is1
         is1=is2
         is2=i
      end if

      if(firsttime) then      
         call io_assign( iu )
         call io_assign( iuL )
         call io_assign( iuR )
         open( iu, file='COOP.dat', status='unknown' )   
         open( iuL, file='COOPL.dat', status='unknown' )   
         open( iuR, file='COOPR.dat', status='unknown' )   
         
         io1=is1 + NBUFL + NGL
         io2=is2 + NBUFL + NGL
         do ia=1,nua
            if(io1.le.lasto(ia) .and. io1.gt.lasto(ia-1)) iastart=ia !orbital is1 is on atom ia1
            if(io2.le.lasto(ia) .and. io2.gt.lasto(ia-1)) iaend=ia !orbital is1 is on atom ia1
         end do                 !ia

         write(ju,*) 'ATOMIC COOP: starting on atom: ',iastart,
     .               ' ending on atom: ',iaend 
         firsttime=.false.
      end if


 10   format(i4,i4,1x,F10.4,20(f8.4,1x))
 11   format(a4,i4,1x,F10.4,20(f8.4,1x))
 12   format(a4,i4,i4,1x,F10.4,20(f8.4,1x))

      ztmp = dcmplx(-2d0/Pi,0d0)/eV


      do ia1=iastart,iaend
         do ia2=iastart,iaend

               coopL = 0d0
               coopR = 0d0
               coopTot  = 0d0
               
c     ---------------------------------            

            ns1 = lasto(ia1) - lasto(ia1-1)
            ns2 = lasto(ia2) - lasto(ia2-1)
            

            do il1=1,ns1
               io1 = il1 + lasto(ia1-1)
               i1  = io1 - NBUFL - NGL

               do il2=1,ns2
                  io2 = il2 + lasto(ia2-1)
                  i2  = io2 - NBUFL - NGL
            
                  coopTot  = coopTot  + DIMAG(ztmp*GF(i1,i2)*S(i1,i2))
                  coopR = coopR + DIMAG(ztmp*GFR(i1,i2)*S(i1,i2))

               end do           !il2 -- orbital number on atom ia2

               write(iuL,10)  ia1,ia2,energy*eV,
     .              ( DIMAG(ztmp*(GF(i1,i2)-GFR(i1,i2))*S(i1,i2)),
     .              i2 = 1 + lasto(ia2-1) - NBUFL - NGL,
     .              ns2 + lasto(ia2-1) - NBUFL - NGL )

               write(iuR,10)  ia1,ia2,energy*eV,
     .              ( DIMAG(ztmp*(GFR(i1,i2))*S(i1,i2)),
     .              i2 = 1 + lasto(ia2-1) - NBUFL - NGL,
     .              ns2 + lasto(ia2-1) - NBUFL - NGL )

            end do              !i1 -- orbital number on atom ia1
            
            coopL = coopTot  - coopR
            write(iu,12)  'AO: ',ia1,ia2,energy*eV,coopTot,coopR,coopL
            
         end do                 !ia2
      end do                    !ia1





c=============================================================
c summedup COOP to all left/right atoms w.r.t. the given atom
c=============================================================

      do ia1=iastart,iaend
         ns1 = lasto(ia1)-lasto(ia1-1)            
         do il1=1,ns1
            io1 = il1 + lasto(ia1-1)
            i1  = io1 - NBUFL - NGL
            
            coopL2L = 0d0
            coopL2R = 0d0
            coopR2L = 0d0
            coopR2R = 0d0
            coop2L =  0d0
            coop2R =  0d0
c     ---------------------------------            

            do ia2=1,nua
               ns2 = lasto(ia2)-lasto(ia2-1)
               do il2=1,ns2
                  io2 = il2 + lasto(ia2-1)
                  i2  = io2 - NBUFL - NGL
                  if(i2.gt.0 .and. i2.le.nc) then

                     if(ia2.lt.ia1) then !atom to the Left
                        coop2L = coop2L + 
     .                       DIMAG(ztmp*GF(i1,i2)*S(i1,i2))
                        coopR2L = coopR2L + 
     .                       DIMAG(ztmp*GFR(i1,i2)*S(i1,i2))
                     end if     !atom 2 Left

                     if(ia2.gt.ia1) then !atom to the Right
                        coop2R = coop2R + 
     .                       DIMAG(ztmp*GF(i1,i2)*S(i1,i2))
                        coopR2R = coopR2R + 
     .                       DIMAG(ztmp*GFR(i1,i2)*S(i1,i2))
                     end if     !atom 2 Right

                  end if        ! nc+1 > i2 > 0
               end do           !il2 -- orbital number on atom ia2
            end do              !ia2
c     ---------------------------------            
         end do                 !i1 -- orbital number on atom ia1
            
         coopL2L = coop2L - coopR2L
         coopL2R = coop2R - coopR2R
         write(iu,11) 'LR: ',
     .        ia1,energy*eV,coopL2L,coopL2R,coopR2L,coopR2R
      end do                    !ia1



c     call io_close( iu )

      return
      end





