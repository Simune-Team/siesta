      subroutine mulliken(iopt,nspin,natoms,nbasistot,maxnh,numh,
     .                    listhptr,listh,s,dm,isa,lasto,iaorb,
     .                    iphorb)
C ********************************************************************
C Subroutine to perform Mulliken population analysis.
C (Overlap and total populations, both for orbitals and for atoms)
C The density matrix (d.m.) and overlap matrix are passed in sparse form
C (both with the same sparse structure)
C There is no output. The populations are printed to the output.
C
C Written by P.Ordejon, October'96
C Non-collinear spin added by J.M.Soler, May 1998. 
C Symmetry label for each orbital included by DSP, Oct. 1998.
C Label with the principal quantum number introduced by DSP, July 1999.
C ************************** INPUT ************************************
C integer iopt                : Work option: 1 = atomic and orbital charges
C                                            2 = 1 + atomic overlap pop.
C                                            3 = 2 + orbital overlap pop.
C integer nspin               : Number of spin components
C integer natoms              : Number of atoms
C integer nbasistot           : Number of basis orbitals over all nodes
C integer maxnh               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements
C integer numh(nbasis)        : First Control vector of d.m. and overlap
C integer listhptr(nbasis)    : Second Control vector of d.m. and overlap
C integer listh(maxnh)        : Third Control vector of d.m. and overlap
C real*8  s(maxnh)            : Overlap matrix in sparse form
C real*8  dm(maxnh,nspin)     : Density matrix in sparse form 
C integer isa(natoms)         : Species index of each atom
C integer lasto(0:maxa)       : Index of last orbital of each atom
C                               (lasto(0) = 0) 
C integer iaorb(nbasis)       : Atomic index of each orbital
C integer iphorb(nbasis)      : Orbital index of each orbital in its atom
C ************************* OUTPUT *************************************
C No output. The results are printed to standard output
C **********************************************************************
C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs, GlobalToLocalOrb
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer
     .  iopt,natoms,nbasistot,maxnh,nspin

      integer
     .  numh(*),lasto(0:natoms),listh(maxnh),listhptr(*),
     .  iphorb(*), isa(natoms), iaorb(*)

      real(dp)
     .  dm(maxnh,nspin),s(maxnh)

      external
     .  memory
C Internal parameters ..................................................
C Number of culumns in printout.  Must be smaller than 20
      integer ncol, nbasis
      parameter (ncol = 8)

#ifdef MPI
      integer MPIerror
      real(dp)
     .  pb(ncol)
      real(dp), dimension(:), allocatable, save :: qb
#endif
      integer i,ia,ib,ii,imax,in,io,ior,ip,is,j,ja,jja,jo,jor,
     .  ind, nao, nblock, ns, ispec, irow, nrow, nres, 
     .  config(ncol), config2(8), iNode, itot

      real(dp)
     .  p(ncol),qa,qas(4),qts(4),qtot,stot,svec(3) 

      real(dp), dimension(:),   allocatable, save :: qo
      real(dp), dimension(:,:), allocatable, save :: qos

      character sym_label(ncol)*7, atm_label*20,
     .     sym_label2(8)*7   

C ......................

C Get Node number
#ifdef MPI
C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(nbasistot,Node,Nodes,nbasis)
#else
      nbasis = nbasistot
#endif

      if (iopt.eq.0) then
C iopt = 0 implies no analysis
        return
      elseif (iopt.lt.0 .or. iopt.gt.3) then
        if (Node.eq.0) then
          write(6,"(a)") 'mulliken: ERROR: Wrong iopt'
        endif
        return
      endif 

C Allocate local memory
      allocate(qos(nspin,nbasis))
      call memory('A','D',nspin*nbasis,'mulliken')
      allocate(qo(nbasis))
      call memory('A','D',nbasis,'mulliken')
#ifdef MPI
      allocate(qb(nbasis))
      call memory('A','D',nbasis,'mulliken')
#endif

      ns=0
      do i = 1,natoms
        ns=max(ns,isa(i))
      enddo 

C Compute and print Overlap Populations for Orbitals ....................
      if (iopt .eq. 3) then
        if (Node.eq.0) then
          write(6,*) 
          write(6,"(a)")'mulliken: Overlap Populations between Orbitals'
        endif
        nblock = nbasistot / ncol
        ip=1
        if (nblock*ncol .eq. nbasistot) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = nbasistot - nblock * ncol  
          do ii=1,imax 
             sym_label(ii)=symfio(isa(iaorb((ib-1)*ncol+ii)),
     .                iphorb((ib-1)*ncol+ii))
             config(ii)=cnfigfio(isa(iaorb((ib-1)*ncol+ii)),
     .                iphorb((ib-1)*ncol+ii))
          enddo 
          if (Node.eq.0) then
            write(6,*) 
            write(6,'(14x,20(2x,i4,3x))')((ib-1)*ncol+ii,ii=1,imax) 
            write(6,'(17x,20(1x,i1,a7))') 
     .              (config(ii),sym_label(ii),ii=1,imax)
          endif

C Loop over the following section so that each Node can print
C the terms relating to the local orbitals in sequence. Not
C the most efficient algorithm, but preserves order and keeps
C memory usage the same as in the serial version.
          do itot = 1,nbasistot
            call GlobalToLocalOrb(itot,Node,Nodes,i)
            if (i.gt.0) then
              sym_label(1)=symfio(isa(iaorb(itot)),iphorb(itot))
              config(1)=cnfigfio(isa(iaorb(itot)),iphorb(itot))
              do ii = 1,imax
                p(ii) = 0.0d0
              enddo
              do in = 1,numh(i)
                ind = listhptr(i)+in
                j = listh(ind)
                ii = j - (ib - 1) * ncol
                if (ii .ge. 1 .and. ii .le. imax) then
                  p(ii) = 0.d0
                  do is = 1,min(nspin,2)
                    p(ii) = p(ii) + dm(ind,is) * s(ind)
                  enddo
                endif
              enddo
              write(6,15) itot,config(1), sym_label(1),(p(ii),ii=1,imax)
            endif
#ifdef MPI
            call MPI_Barrier( MPI_Comm_World, MPIerror )
#endif
          enddo
        enddo
      endif
C ...................

C Compute and print Overlap Populations for Atoms ....................
      if (iopt .ge. 2) then
        if (Node.eq.0) then
          write(6,*) 
          write(6,"(a)")'mulliken: Overlap Populations between Atoms'
        endif
        nblock = natoms / ncol
        ip=1
        if (nblock*ncol .eq. natoms) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = natoms - nblock * ncol
          if (Node.eq.0) then
            write(6,*) 
            write(6,10) ((ib-1)*ncol+ii,ii=1,imax)
          endif
          do i = 1,natoms
            do ii = 1,imax
              p(ii) = 0.0
            enddo
            do ior = lasto(i-1)+1,lasto(i)
              call GlobalToLocalOrb(ior,Node,Nodes,itot)
              if (itot.gt.0) then
                do in = 1,numh(itot)
                  ind = listhptr(itot)+in
                  jor = listh(ind)
                  do jja = 1,natoms
                    if (lasto(jja) .ge. jor) then
                      ja = jja
                      goto 100
                    endif
                  enddo
                  goto 110
 100              ii = ja - (ib - 1) * ncol
                  if (ii .ge. 1 .and. ii .le. imax) then
                    do is = 1,min(nspin,2)
                      p(ii) = p(ii) + dm(ind,is) * s(ind)
                    enddo
                  endif
 110              continue
                enddo
              endif
            enddo

#ifdef MPI
C Global sum of values stored in p
            pb(1:imax)=p(1:imax)
#ifdef NODAT
            call MPI_Reduce(pb,p,imax,MPI_double_precision,MPI_sum,
     .        0,MPI_Comm_World,MPIerror)
#else
            call MPI_Reduce(pb,p,imax,DAT_double,MPI_sum,
     .        0,MPI_Comm_World,MPIerror)
#endif
#endif

            if (Node.eq.0) then
              write(6,11) i,(p(ii),ii=1,imax)
            endif
          enddo
        enddo
      endif
C ....................

C Compute and print Mulliken Orbital and Atomic Populations ..........
      if (iopt .ge. 1) then
        if (Node.eq.0) then
          write(6,*) 
          write(6,"(a)")'mulliken: Atomic and Orbital Populations:'
        endif
        if (nspin .le. 2) then
          do is = 1,nspin
            if (nspin .eq. 2.and.Node.eq.0) then
              if(is .eq. 1) write(6,'(/,a)') 'mulliken: Spin UP '
              if(is .eq. 2) write(6,'(/,a)') 'mulliken: Spin DOWN '
            endif
            qtot = 0.0d0
            do ispec = 1,ns  

             atm_label = labelfis(ispec)
             if (Node.eq.0) then
               write(6,'(/2a)')'Species: ', atm_label 
               write(6,'(a4,a7,a6)') 'Atom', 'Qatom', 'Qorb'
             endif
C DSP, Writing symmetries for each orbital. 
C DSP, Orbitals with a 'P' belong to the polarization shell
             nao = nofis(ispec)
             nrow=nao/8  
             nres=nao-8*nrow 
             nao=0
             do irow=1,nrow 
               do io=1,8  
                 nao=nao+1
                 sym_label2(io)=symfio(ispec,nao) 
                 config2(io)=cnfigfio(ispec,nao)
               enddo 
               if (Node.eq.0) then
                 write(6,'(15x,8(i1,a7))')
     .             (config2(io),sym_label2(io),io=1,8)
               endif
             enddo 
             do io=1,nres 
               nao=nao+1 
               sym_label2(io)=symfio(ispec,nao) 
               config2(io)=cnfigfio(ispec,nao)
             enddo 
             if (Node.eq.0) then
               write(6,'(15x,8(i1,a7))')
     .           (config2(io),sym_label2(io),io=1,nres) 
             endif

             do ia = 1,natoms
               if (isa(ia).eq.ispec) then
C             Compute charge in each orbital of atom ia
                 qa = 0.0d0
                 do io = lasto(ia-1)+1,lasto(ia)
                   nao = io - lasto(ia-1)
                   qo(nao) = 0.0d0
                   call GlobalToLocalOrb(io,Node,Nodes,itot)
                   if (itot.gt.0) then
                     do in = 1,numh(itot)
                       ind = listhptr(itot) + in
                       qo(nao) = qo(nao) + dm(ind,is) * s(ind)
                     enddo
                     qa = qa + qo(nao)
                   endif
                 enddo
                 qtot = qtot + qa
#ifdef MPI
C Global sum of values stored in p
                 qb(1:nao)=qo(1:nao)
#ifdef NODAT
                 call MPI_Reduce(qb,qo,nao,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
#else
                 call MPI_Reduce(qb,qo,nao,DAT_double,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
#endif
                 qb(1)=qa
#ifdef NODAT
                 call MPI_Reduce(qb(1),qa,1,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
#else
                 call MPI_Reduce(qb(1),qa,1,DAT_double,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
#endif
#endif
                 if (Node.eq.0) then
                   write(6,'(i4,f7.3,8f8.3,(/11x,8f8.3))')
     .               ia, qa, (qo(io),io=1,nao) 
                 endif
               endif 
             enddo 
            enddo
#ifdef MPI
C Global sum of total charge
            qb(1)=qtot
#ifdef NODAT
            call MPI_Reduce(qb(1),qtot,1,MPI_double_precision,MPI_sum,
     .        0,MPI_Comm_World,MPIerror)
#else
            call MPI_Reduce(qb(1),qtot,1,DAT_double,MPI_sum,
     .        0,MPI_Comm_World,MPIerror)
#endif
#endif
            if (Node.eq.0) then
              write(6,"(/a,f8.3)") 'mulliken: Qtot = ', qtot
            endif
          enddo
        elseif (nspin .eq. 4) then
          do is = 1,nspin
            qts(is) = 0.0d0
            do io = 1,nbasis
              qos(is,io) = 0.0d0
            enddo
          enddo
          do is = 1,nspin
            do io = 1,nbasis
              do in = 1,numh(io)
                ind = listhptr(io)+in
                jo = listh(ind)
                qos(is,io) = qos(is,io) + dm(ind,is)*s(ind)/2
                if (jo.le.nbasis)
     .            qos(is,jo) = qos(is,jo) + dm(ind,is)*s(ind)/2
              enddo
            enddo
          enddo 
          do ispec=1,ns
            atm_label=labelfis(ispec)
            if (Node.eq.0) then
              write(6,'(/2a)')'Species: ', atm_label
              write(6,'(/,a4,a9,4x,2a10,3x,a8,/,64(1h-))')
     .        'Atom', 'Orb', 'Charge', 'Spin', 'Svec'
            endif

           do ia = 1,natoms 
            if (isa(ia).eq.ispec) then 
            qas(1:nspin) = 0.d0
            do iNode = 0,Nodes-1
              if (iNode.eq.Node) then
                do io = lasto(ia-1)+1,lasto(ia)  
                  call GlobalToLocalOrb(io,Node,Nodes,itot)
                  if (itot.gt.0) then
C DSP, Writing symmetries for each orbital.
C DSP, Orbitals with a 'P' belong to the polarization shell

                    sym_label(1)=symfio(ispec,iphorb(io))
                    config(1)=cnfigfio(ispec,iphorb(io))
                    do is = 1,nspin
                      qas(is) = qas(is) + qos(is,itot)
                      qts(is) = qts(is) + qos(is,itot)
                    enddo
                    call spnvec( nspin, qos(1,itot), qtot, stot, svec)
                    write(6,'(i4,i5,i1,a7,2f10.5,3x,3f8.3)')
     .                ia, io,config(1),sym_label(1), qtot, stot, svec
                  endif
                enddo
              endif
#ifdef MPI
              call MPI_Barrier( MPI_Comm_World, MPIerror )
#endif
            enddo
            call spnvec( nspin, qas, qtot, stot, svec )
#ifdef MPI
C Global reduction of terms
            qb(1)=qtot
            qb(2)=stot
            qb(3)=svec(1)
            qb(4)=svec(2)
            qb(5)=svec(3)
#ifdef NODAT
            call MPI_Reduce(qb,pb,5,MPI_double_precision,MPI_sum,0,
     .        MPI_Comm_World,MPIerror)
#else
            call MPI_Reduce(qb,pb,5,DAT_double,MPI_sum,0,
     .        MPI_Comm_World,MPIerror)
#endif
            qtot=pb(1)
            stot=pb(2)
            svec(1)=pb(3)
            svec(2)=pb(4)
            svec(3)=pb(5)
#endif
            if (Node.eq.0) then
              write(6,'(i4,4x,a6,3x,2f10.5,3x,3f8.3,/)')
     .          ia, 'Total', qtot, stot, svec 
            endif
            endif
           enddo
          enddo 
          call spnvec( nspin, qts, qtot, stot, svec )
#ifdef MPI
C Global reduction of terms
          qb(1)=qtot
          qb(2)=stot
          qb(3)=svec(1)
          qb(4)=svec(2)
          qb(5)=svec(3)
#ifdef NODAT
          call MPI_Reduce(qb,pb,5,MPI_double_precision,MPI_sum,0,
     .        MPI_Comm_World,MPIerror)
#else
          call MPI_Reduce(qb,pb,5,DAT_double,MPI_sum,0,
     .        MPI_Comm_World,MPIerror)
#endif
          qtot=pb(1)
          stot=pb(2)
          svec(1)=pb(3)
          svec(2)=pb(4)
          svec(3)=pb(5)
#endif
          if (Node.eq.0) then
            write(6,'(64(1h-),/,2a8,1x,2f10.5,3x,3f8.3,/)')
     .        'Total', 'Total', qtot, stot, svec
          endif
        endif
      endif

C Deallocate local memory
      call memory('D','D',size(qos),'mulliken')
      deallocate(qos)
      call memory('D','D',size(qo),'mulliken')
      deallocate(qo)
#ifdef MPI
      call memory('D','D',size(qb),'mulliken')
      deallocate(qb)
#endif

C ...................

10    format(12x,20(2x,i4,2x))
11    format(i12,20(1x,f7.3))
12    format('i = ',i4,'   q = ',f6.3,'   q_orb = ',9f6.3)
13    format(32x,9f6.3) 
15    format(i4,1x,i1,a7,20(1x,f8.3))
      return
      end



      subroutine spnvec( ns, qs, qt, st, sv )
c ********************************************************************
c Finds the spin vector components from the spin density matrix
c Written by J.M.Soler, May 1998.
c ******* Input ******************************************************
c integer ns     : Number of components in spin density matrix
c real*8  qs(ns) : Spin density matrix elements with the convention
c                  is=1 => Q11; is=2 => Q22; is=3 => Real(Q12);
c                  is=4 => Imag(Q12)
c ******* Output *****************************************************
c real*8  qt    : Total charge
c real*8  st    : Total spin
c real*8  sv(3) : Spin vector
c ********************************************************************

      use precision

      implicit          none
      integer           ns
      real(dp)          qs(ns), qt, st, sv(3)
      real(dp)          cosph, costh, sinph, sinth, tiny
      parameter ( tiny = 1.d-12 )

      if (ns .eq. 1) then
        qt = qs(1)
        st = 0.0d0
        sv(1) = 0.0d0
        sv(2) = 0.0d0
        sv(3) = 0.0d0
      elseif (ns .eq. 2) then
        qt = qs(1) + qs(2)
        st = qs(1) - qs(2)
        sv(1) = 0.0d0
        sv(2) = 0.0d0
        sv(3) = st
      elseif (ns .eq. 4) then
        qt = qs(1) + qs(2)
        st = sqrt( (qs(1)-qs(2))**2 + 4.d0*(qs(3)**2+qs(4)**2) )
        costh = ( qs(1) - qs(2) ) / ( st + tiny )
        sinth = sqrt( 1.0d0 - costh**2 )
        cosph =  qs(3) / ( sqrt( qs(3)**2 + qs(4)**2 ) + tiny )
        sinph = -qs(4) / ( sqrt( qs(3)**2 + qs(4)**2 ) + tiny )
        sv(1) = st * sinth * cosph
        sv(2) = st * sinth * sinph
        sv(3) = st * costh
      else
        write(6,*) 'spnvec: ERROR: invalid argument ns =', ns
        return
      endif
      end
