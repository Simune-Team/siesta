! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine KSV_init( ucell, maxkpol, nkpol,kpol,wgthpol)
C *********************************************************************
C Finds polarization using the method of King-Smith and Vanderbilt
C ( Geometric Berry phase). Initialisation routine.
C Written by DSP, March 1999.
C **************************** INPUT **********************************
C real*8   ucell(3,3)         : Unit cell vectors
C integer  maxkpol            : Last dimension of kpoint 
C *************************** INPUT/OUTPUT ****************************
C integer  nkpol              : Maximum number of grid points for the
C                               bidimensional integrals
C *************************** OUTPUT **********************************
C real*8   kpol(3,maxkpol)    : Auxiliar array to store the kpoints
C                               for the bidimensional integrations.
C real*8   wgthpol(maxkpol)   : Auxiliar array to store the weights of
C                               the kpoints.            
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C *********************************************************************
      implicit          none
      integer           maxkpol, nkpol
      double precision  kpol(3,maxkpol), ucell(3,3),
     .                  wgthpol(maxkpol)
C *********************************************************************

C Internal variables 
      integer
     .  ix, iy, kscell(3,3), igrd, nk, nmeshk(3,3)
         
      double precision
     .  rcell(3,3), displ(3), dsp(3), cutoff

      external
     .  reclat, timer

C Start time counter 
      call timer( 'KSV_init', 1 )

C Reading unit cell and calculate the reciprocal cell
      call reclat( ucell, rcell, 1 )

C Find the integration grids in reciprocal space   
C In the first call this is dome to set the rigth 
C dimension for the arrays which depend on parameter
C nkpol
      call repol(nmeshk,dsp) 
      
      nk=nmeshk(1,1)*nmeshk(2,1)*nmeshk(3,1)+
     .  nmeshk(1,2)*nmeshk(2,2)*nmeshk(3,2)+
     .  nmeshk(1,3)*nmeshk(2,3)*nmeshk(3,3) 

      if (nk.ne.0) then  
        nkpol=1
        do igrd=1,3
          cutoff=0.0d0
          do ix=1,3
            do iy= 1,3
              kscell(ix,iy)=0
            enddo
          enddo
          do ix=1,3
            if (ix.ne.igrd) then
              kscell(ix,ix)= nmeshk(ix,igrd)
              displ(ix)=dsp(igrd)
            else 
              kscell(ix,ix)=1
              displ(ix)=0.0d0
            endif
          enddo 
          nk=kscell(1,1)*kscell(2,2)*kscell(3,3) 
          if (nk.gt.0) then 
            call kgridinit( ucell, kscell, displ, cutoff, nk )
            if (maxkpol.gt.nk) then
              call kgrid( ucell, kscell,  displ,
     .                    nk, kpol, wgthpol )
            endif
          endif 
          nkpol=max(nk,nkpol)       
        enddo  
      else
        nkpol=0
      endif  

      call timer( 'KSV_init', 2 )

      end
