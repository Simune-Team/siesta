      subroutine detover(psiprev,psi,maxpsi, S, Sr, 
     .               numh, listh, indxuo,
     .               no, nuo, xij, maxno, maxuo, nocc,
     .               kpoint, dk, detr,deti)
C *********************************************************************
C Finds the determinant of the overlap matrix
C between the periodic Bloch functions corresponding to neighboring  
C k points
C Written by DSP. March 1999
C **************************** INPUT ********************************** 
C real*8  psiprev(maxpsi)     : Wavefunctions in previous k point
C real*8  psi(maxpsi)         : Wavefunctions in current k point
C real*8  maxpsi              : Array dimensions for psi and psiprev
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C real*8  Sr(maxno,maxuo)     : Position operator matrix elements (sparse)
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listh(maxno,nuo)    : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row 
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuo                 : Number of basis orbitals in the unit cell
C integer maxno               : Maximum number of orbitals interacting
C                               with any orbital 
C integer maxuo               : Third dimension of xij
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)  
C integer nocc                : number of occupied states
C real*8  kpoint(3)           : Current kpoint
C real*8  dk(3)               : Vector joining the previous and current 
C                               kpoint
C *************************** INPUT/OUTPUT ****************************
C real*8  detr                  : Real part of the determinant
C real*8  deti                  : Imaginary part of the determinant
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      implicit none
 
      integer nuo, maxno, maxuo, maxpsi, no,
     .  listh(maxno,nuo), numh(nuo), indxuo(no),nocc,
     .  info, job
      parameter(job=10)

      double precision psiprev(maxpsi), dk(3),detr, deti,
     .  xij(3,maxno,maxuo), S(maxno,maxuo),  Sr(maxno,maxuo),
     .  psi(maxpsi), kpoint(3)

      complex*16  det(2), pipj, determ

C**** Internal variables ***********************************************

      integer iuo, juo, j, io, ie, je, iu, ju, jo
      double precision  kxij, skxij, ckxij,pipj1, pipj2
      integer i

      complex*16, dimension(:,:), allocatable ::
     .  Aux, Aux2

      integer, dimension(:), allocatable ::
     .  Auxint

C Allocate local memory
      allocate(Aux(nuo,nuo))
      call memory('A','Z',nuo*nuo,'detover')
      allocate(Aux2(nuo,nuo))
      call memory('A','Z',nuo*nuo,'detover')
      allocate(Auxint(nuo))
      call memory('A','I',nuo,'detover')

      do iuo = 1,nuo
        do juo = 1,nuo
          Aux2(juo,iuo) = dcmplx(0.d0,0.d0)
          Aux(juo,iuo) = dcmplx(0.d0,0.d0)
        enddo
      enddo 

      do io = 1,nuo
        do j = 1,numh(io)
          jo = listh(j,io)
          iuo = indxuo(io)
          juo = indxuo(jo)
          kxij = (kpoint(1)-0.5d0*dk(1)) * xij(1,j,io) +
     .           (kpoint(2)-0.5d0*dk(2)) * xij(2,j,io) +
     .           (kpoint(3)-0.5d0*dk(3)) * xij(3,j,io) 
          ckxij = dcos(kxij)
          skxij = dsin(kxij)
          Aux2(iuo,juo)=Aux2(iuo,juo)+
     .    dcmplx(  S(j,io)*ckxij + Sr(j,io)*skxij, 
     .        S(j,io)*skxij - Sr(j,io)*ckxij  )  
 
        enddo 
      enddo 
 
      do ie = 1,nocc
        do je = 1, nocc
              
          do iuo=1,nuo
            do juo=1,nuo    

              iu=2*nuo*(ie-1)+2*(iuo-1)
              ju=2*nuo*(je-1)+2*(juo-1) 

              pipj1 = psiprev(iu+1) * psi(ju+1) +
     .                psiprev(iu+2) * psi(ju+2)
              pipj2 = psiprev(iu+1) * psi(ju+2) -
     .                psiprev(iu+2) * psi(ju+1) 
              pipj=dcmplx(pipj1,pipj2)

              Aux(ie,je) = Aux(ie,je) + 
     .                 pipj*Aux2(iuo,juo)

            enddo 
          enddo  

        enddo 
      enddo 

      call zgefa(Aux,nuo,nocc,Auxint,info)
      call zgedi(Aux,nuo,nocc,Auxint,det,Aux2,job)

      determ=det(1)
      deti=dimag(determ)
      detr=dreal(determ)   

C Deallocate local memory
      call memory('D','Z',size(Aux),'detover')
      deallocate(Aux)
      call memory('D','Z',size(Aux2),'detover')
      deallocate(Aux2)
      call memory('D','I',size(Auxint),'detover')
      deallocate(Auxint)

      end
