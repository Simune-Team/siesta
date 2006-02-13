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
      subroutine obc(polxyz, polR, ucell, dx, nspin, node)
C *******************************************************************
C Writes the born effective charges to file 
C by Tom Archer
C ********* INPUT ***************************************************
C real*8 polxyz(3,nspin): polarization along cartisian coordinates(Bohr)
C real*8 dx             : atomic displacements (in Bohr)
C real*8 polR(3,nspin)  : polarization along lattice vectors(Bohr)
C real*8 ucell(3,3)     : cell vectors
C integer nspin         : spin polarized calculation flag
C integer node          : node on which process is running
C ********** BEHAVIOUR **********************************************
C On the first call (undisplaced coordinates), initial polarization  
C is calculated.
C The born charges are calculated as the diffrence between the 
C polarization of the undisplaced and the displaced coordinates 
C divided by the displacment
C 
C a phase is introduced in the KSV calculation which is removed before
C the BC matrix is saved.  
C *******************************************************************

      use fdf

      implicit          none

c External
      integer           
     .     node, nspin 
      double precision  
     .     dot,polxyz(3,nspin), polR(3,nspin), ucell(3,3), dx
      external          
     .     chkdim, io_assign, io_close, paste, timer,memory, dot


c Intrnal

      integer           
     .     igrd, ispin,  i, ix, unit1, nwritten, n

      character 
     .     fname*33, sname*30, line*132, paste*33

      logical   frstme

      double precision 
     .     rdummy, dmod(3),
     .     phaseR(3), uR(3,3), pR(3), pxyz(3)

      double precision, dimension(:,:), allocatable, save :: pres

      save      frstme, fname, nwritten
      data      frstme /.true./
      data      nwritten / 0 /
c###################################################################

c     Allocate local array for storing residual forces
      if (.not.allocated(pres)) then
        allocate(pres(3,nspin))
        call memory('A','D',3*nspin,'obc')
      endif

c     Find file name
      if (frstme) then
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.BC')
      endif

      call io_assign(unit1)
      open( unit1, file=fname, status='unknown' )
      rewind(unit1)

c **************matrix to convert from lattice vectors to xyz
      do igrd=1,3
        dmod(igrd)=dot(ucell(1,igrd),ucell(1,igrd),3)
        dmod(igrd)=dsqrt(dmod(igrd))
        do ix=1,3
          uR(ix,igrd)=ucell(ix,igrd)/dmod(igrd)
        enddo
      enddo
c ***********end matrix to convert from lattice vectors to xyz


      if (node.eq.0) then
c **************print polarization***********************************
         do ispin=1,nspin
            if (nspin.gt.1) then
               if (ispin.eq.1) write(6,'(/,a)')
     .              'obc: Macroscopic polarization for spin Up:'
               if (ispin.eq.2) write(6,'(/,a)')
     .              'obc: Macroscopic polarization for spin Down:'
            endif
            write(6,'(/,a,a)')
     .       'obc: Macroscopic polarization per unit cell',
     .       ' along the lattice vectors (Bohr): '
            write(6,'(a,3f12.6)') 'obc:',(polR(ix,ispin),ix=1,3)
c            write(6,'(/,a)')'obc: Along cartesian directions '
c            write(6,'(3f12.6)') (polxyz(ix,ispin),ix=1,3)
         enddo 
c **************end print polarization********************************
      endif

      if (node.eq.0) then
c write BC matrix to file *******************************************
         if (frstme) then
c     Write header message if frstme
            if (nspin.eq.1) then
               write(unit1,'(a)') 'BC matrix'
            else
               write(unit1,'(a)') 'BC matrix for spin-up + spin-down'
            endif  
c     Set values of residual polarization
            do ispin=1,nspin
               do igrd=1,3
                  pres(igrd,ispin) = polR(igrd,ispin)
               enddo
            enddo
            frstme = .false.
            call io_close(unit1)
            return
         endif

c     Read file written so far to put pointer for write in the correct place
         read(unit1,'(a)') line
         do n = 1,nwritten            
           read(unit1,'(3f15.7)') (rdummy, ix=1,3)
         enddo
     
c     write BC matrix in lattice vectors
         pR(:)=0.0d0
         pxyz(:)=0.0d0
         do ispin=1,nspin
            do igrd=1,3
               pR(igrd)=pR(igrd)+(polR(igrd,ispin)-pres(igrd,ispin))
            enddo
         enddo
         
c *******remove extra phase******************************************
      do igrd=1,3
         phaseR(igrd) = nint((pR(igrd)/dmod(igrd))) * dmod(igrd)
         pR(igrd) = pR(igrd) - phaseR(igrd) 
      enddo

c     get polarization in cartesian coordinates
      do ix=1,3
         pxyz(ix) = uR(ix,1)*pR(1)+ uR(ix,2)*pR(2)+uR(ix,3)*pR(3)
      enddo
c*******end remove phase*********************************************
         write(unit1,'(3f15.7)') (pxyz(ix)/dx,ix=1,3)

         nwritten = nwritten + 1
         
         call io_close(unit1)
c**  End write BC matrix to file  **************************
      endif
      
      return
      end
