








        subroutine phiatm(is,io,r,phi,grphi)
C***********************************************************************
C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C  Written by D,Sanchez-Portal. May, 1996.
C*************************INPUT*****************************************
C  INTEGER IS : Species index
C  INTEGER IO : Orbital index (within atom):
C              IO > 0 =>  Basis orbitals
C              IO = 0 =>  Local screened pseudopotential
C              IO < 0 =>  Kleynman-Bylander projectors
C  REAL*8 R(3): Point vector, relative to atom 
C*************************OUTPUT****************************************
C  REAL*8 PHI     : Value of the basis orbital, KB-projector or local
C                    pseudopotential
C  REAL*8 GRPHI(3): Gradient of the basis orbital, KB-projector or local
C                    pseudopotential
C*************************UNITS*****************************************
C Distances in Bohr
C************************BEHAVIOUR**************************************
C 0) Before using PHIATM, the pseudopotential must be initialized 
C    by calling ATOM for each atomic species required
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS and/or IO
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) PHIATM with IO = 0 is strictly equivalent to VLOCAL
C***********************************************************************
        implicit double precision (a-h,o-z)

        include 'atom.h'
        
        double precision
     .      r(3),table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .      tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .      phi,grphi(3),rly(lmx2), grly(3,lmx2)

        integer
     .      is,io,nzettb(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .      nkbmax(nsmax),loctab(nsmax)


        common/cmtab/table
        common/cmspline/tab2
        common/cmzeta/nzettb
        common/control/ismax,nomax,nkbmax
        common/cmloc/loctab


        if ((is.lt.1).or.(is.gt.ismax)) then
           write(6,*) 'PHIATM: THERE ARE NO DATA FOR IS=',IS
           write(6,*) 'PHIATM: ISMIN= 1, ISMAX= ',ismax
           STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'PHIATM: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'PHIATM: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif


       if (io.gt.0) then

          norb=0
          do  l=0,lmaxd
            do izeta=1,nzettb(l,is)
               norb=norb+(2*l+1)
               if(norb.ge.io) goto 20
            enddo
          enddo

  20      lorb=l
          izetorb=izeta
          morb=io-norb+lorb

          indx=0
          do l=0,lorb-1
             indx=indx+nzettb(l,is)
          enddo

          indx=indx+izetorb
  
         elseif(io.lt.0) then
         lloc=loctab(is) 
         nkb=0
         indx=0
         do 30 l=0,lmaxd
            if(l.eq.lloc) goto 30
            indx=indx+1
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 40
 30      continue


 40      lorb=l
         morb=-io+nkb+lorb 
         indx=-indx
      
         elseif (io.eq.0) then

             indx=0 
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif

        delt=table(1,indx,is) 
        rmax=delt*(ntbmax-1)
  
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1)=0.0d0
           grphi(2)=0.0d0
           grphi(3)=0.0d0

        else

           call splint(delt,table(3,indx,is),tab2(1,indx,is),ntbmax,
     .        rmod,phir,dphidr)

           if (indx.eq.0) then
              phi=phir
              do i=1,3 
                 grphi(i)=dphidr*r(i)/rmod
              enddo
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)
              
              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm)
              enddo


*             write(6,'(a,i4,2f12.6)')
*    .         'phiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                   
        endif


  
        return
        
        end
