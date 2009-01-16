C FDN kscell and kdispl added as dummy
      subroutine TSiohs2( task, finname,gamma,nua,nuotot,notot,nspin,
     .                 maxnh,numh, listhptr, listh, H, S, qtot, temp, 
     .                 xij, indxuo, efs, cell, isa, lasto, xa,
     .                 kscell,kdispl  )
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C Note because of the new more compact method of storing H and S
C this routine is NOT backwards compatible
C *************************** INPUT **********************************
C character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
C logical       gamma         : Is only gamma point used?
C ******************** INPUT or OUTPUT (depending on task) ***********
C integer nuotot              : Number of basis orbitals per unit cell
C integer nspin               : Spin polarization (1 or 2)
C integer maxnh               : First dimension of listh, H, S and
C                               second of xij
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row (-1)
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  qtot                : Total number of electrons
C real*8  temp                : Electronic temperature for Fermi smearing
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not read/written if only gamma point)
C integer nuotot              : Number of basis orbitals per supercell
C integer nua                 : No. atoms in unitcell
C integer isa(nua)            : Species index of each atom
C integer lasto(0:nua)        : Index of last orbital of each atom
C real*8 xa(3,nua)            : Atomic coordinates (Bohr)
C real*8 efs(nspin)           : Fermi energy for each spin
c real*8  cell(3,3)           : Unit cell vectors in real space cell(ixyz,ivec)
C integer indxuo(notot)       : Index of equivalent orbital in unit cell
C                             Unit cell orbitals must be the first in
C                             orbital lists, i.e. indxuo.le.nuo, with
C                             nuo the number of orbitals in unit cell


C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

      implicit          none
      character         task*(*),finname*33
      logical           gamma
      integer           maxnh, nuotot, nspin
c      integer           listh(maxnh),  listhptr(*)
c      double precision  H(maxnh,nspin), S(maxnh),
c     .                  xij(3,maxnh)
      double precision  qtot, temp
      external          io_assign, io_close

      integer, dimension (:), pointer:: listh, listhptr,
     .                        numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs

c Internal variables and arrays
      integer    im, is, iu, j, ju, k, mnh, ns, i,ic
      integer    ih,hl,nut
      logical    found
cccc
      integer nua                
c      integer isa(nua)            
c      integer lasto(0:nua)    
c      real*8 xa(3,nua)    
      real*8  cell(3,3)       
c      integer indxuo(notot)     
      integer notot

C FDN
      integer kscell(3,3)
      real*8  kdispl(3)
      double precision   ef
C FDN

cccc


C Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then

C Check if input file exists
        inquire( file=finname, exist=found )

        if (found) then
C Open file
            call io_assign( iu )
            open( iu, file=finname, status='old' ,
     .       form = 'unformatted')

C FDN kscell and kdispl writen and cell
            read(iu)((kscell(i,j),i=1,3),j=1,3)
            read(iu)(kdispl(i),i=1,3)
            read(iu) cell
C FDN            



C Read dimensions
            read(iu) nua, nuotot, nspin, notot, maxnh
c            write(6,*) 'dim:', nua, nuotot, nspin, notot, maxnh
c  allocations 

            
            allocate(isa(nua))            
            allocate(lasto(0:nua))
            allocate(xa(3,nua))
            allocate(indxuo(notot))
    
            read(iu)(isa(i),i=1,nua)
            read(iu)(lasto(i),i=0,nua)
            read(iu)((xa(j,i),j=1,3),i=1,nua)
            read(iu)(indxuo(i),i=1,notot)

C Allocate local array for global numh
            allocate(numh(nuotot))

           nut = 0
          do ih = 1,nuotot
            read(iu) numh(ih)
             nut = nut + numh(ih)
          enddo

          allocate(listhptr(nuotot))

C Create listhptr
          listhptr(1) = 0
          do hl = 2,nuotot
            listhptr(hl) = listhptr(hl-1) + numh(hl-1)
          enddo

          
          allocate(listh(maxnh))

            ic = 0
          do ih = 1,nuotot
              hl = ih
              do im = 1,numh(hl)
                 ic = ic +1
                read(iu) listh(listhptr(hl)+im)
              enddo
          enddo


C Read Hamiltonian

          allocate(H(maxnh,nspin))

          do is = 1,nspin
            do ih = 1,nuotot
                hl = ih
                do im = 1,numh(hl)
                  read(iu) H(listhptr(hl)+im,is)
                enddo
            enddo
          enddo

C Read Overlap matrix

          allocate(S(maxnh))

          do ih = 1,nuotot
              hl = ih
              do im = 1,numh(hl)
                read(iu) S(listhptr(hl)+im)
              enddo
          enddo
          
            allocate(efs(nspin)) 
           
            read(iu) qtot,temp
C FDN No fixed spin
            read(iu) ef
            efs(1:nspin)=ef
C FDN
C FDN
C            read(iu) cell
C FDN
            read(iu) gamma  

          if (.not.gamma) then
C Read interorbital vectors for K point phasing

            allocate(xij(3,maxnh)) 

            do ih = 1,nuotot
                hl = ih
                do im = 1,numh(hl)
                  read(iu) (xij(k,listhptr(hl)+im),k=1,3)
                enddo
            enddo
          endif

            
c Close file
          call io_close( iu )
        else
            write(6,*) 'iohs: ERROR: file not found: ', finname
            stop 'iohs: ERROR: file not found'
        endif
 
       endif


      end
