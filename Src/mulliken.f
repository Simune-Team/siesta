C $Id: mulliken.f,v 1.11 1999/03/11 19:26:36 ordejon Exp $

      subroutine mulliken( iopt,nspin,natoms,nbasis,nbmax,nhmax,
     .              numh,listh,s,dm,isa,lasto,iaorb,iphorb )
C ********************************************************************
C Subroutine to perform Mulliken population analisys.
C (Overlap and total populations, both for orbitals and for atoms)
C The density matrix (d.m.) and overlap matrix are passed in sparse form
C (both with the same sparse structure)
C There is no output. The populations are printed to the output.
C
C Written by P.Ordejon, October'96
C Non-collinear spin added by J.M.Soler, May 1998. 
C Symmetry label for each orbital included by DSP, Oct. 1998.
C ************************** INPUT ************************************
C integer iopt                : Work option: 1 = atomic and orbital charges
C                                            2 = 1 + atomic overlap pop.
C                                            3 = 2 + orbital overlap pop.
C integer nspin               : Number of spin components
C integer natoms              : Number of atoms
C integer nbasis              : Number of basis orbitals
C integer nbmax               : Second dimension of dm
C integer nhmax               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements of each row
C integer numh(nbasis)        : First Control vector of d.m. and overlap
C integer listh(nhmax,nbasis) : Second Control vector of d.m. and overlap
C real*8 s(nhmax,nbasis)      : Overlap matrix in sparse form
C real*8 dm(nhmax,nbmax,nspin): Density matrix in sparse form 
C integer isa(natoms)         : Species index of each atom
C integer lasto(0:maxa)       : Index of last orbital of each atom
C                               (lasto(0) = 0) 
C integer iaorb(nbasis)       : Atomic index of each orbital
C integer iphorb(nbasis)      : Orbital index of each orbital in its atom
C ************************* OUTPUT *************************************
C No output. The results are printed to standard output
C **********************************************************************
      implicit none
      integer
     .   iopt,natoms,nbasis,nbmax,nhmax,nspin

      integer
     .   numh(nbasis),lasto(0:natoms),listh(nhmax,nbasis),
     .   iphorb(nbasis), isa(natoms), iaorb(nbasis)

      double precision
     .   dm(nhmax,nbmax,nspin),s(nhmax,nbasis)

C Internal parameters ..................................................
C Number of culumns in printout.  Must be smaller than 20
      integer ncol
      parameter (ncol = 8)

C maxo  : Maximum number of orbitals (total)
C normax: Maximum number of orbitals per atom
      integer maxo, normax
      parameter ( maxo   = 2000 ) 
      parameter ( normax =   60 )

      integer i,ia,ib,ii,imax,in,io,ior,ip,is,j,ja,jja,jo,jor,
     .        nao,nblock, ns, ispec, irow, nrow, nres, nofis

      double precision
     .  ap(ncol),p(ncol),qo(normax),qa,qas(4),qos(4,maxo),
     .  qts(4),qtot,stot,svec(3) 

      character symfio*7, sym_label(ncol)*7, atm_label*20,
     .     labelfis*20, sym_label2(8)*7   

      external symfio, labelfis, nofis
C ......................

      if (iopt.eq.0) then
        return
      elseif (iopt.lt.0 .or. iopt.gt.3) then
        write(6,"(a)") 'mulliken: ERROR: Wrong iopt'
        return
      endif 

      ns=0
      do i = 1,natoms
         ns=max(ns,isa(i))
      enddo 
      do i = 1,ncol
        ap(i) = 0.0d0
      enddo

C Compute and print Overlap Populations for Orbitals ....................
      if (iopt .eq. 3) then
        write(6,*) 
        write(6,"(a)")'mulliken: Overlap Populations between Orbitals'
        write(6,*) 
        nblock = nbasis / ncol
        ip=1
        if (nblock*ncol .eq. nbasis) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = nbasis - nblock * ncol  
          do ii=1,imax 
             sym_label(ii)=symfio(isa(iaorb((ib-1)*ncol+ii)),
     .                iphorb((ib-1)*ncol+ii))
          enddo 
          write(6,*) 
          write(6,10) ((ib-1)*ncol+ii,ii=1,imax) 
          write(6,'(15x,20(x,a7))') (sym_label(ii),ii=1,imax)
          write(6,*) 
          do i = 1,nbasis 
            sym_label(1)=symfio(isa(iaorb(i)),iphorb(i))
            do ii = 1,imax
              p(ii) = 0.0
            enddo
            do in = 1,numh(i)
              j = listh(in,i)
              ii = j - (ib - 1) * ncol
              if (ii .ge. 1 .and. ii .le. imax) then
                p(ii) = 0.d0
                do is = 1,min(nspin,2)
                  p(ii) = p(ii) + dm(in,i,is) * s(in,i)
                enddo
              endif
            enddo
            write(6,15) i,sym_label(1),(p(ii),ii=1,imax)
          enddo
        enddo
      endif
C ...................

C Compute and print Overlap Populations for Atoms ....................
      if (iopt .ge. 2) then
        write(6,*) 
        write(6,"(a)")'mulliken: Overlap Populations between Atoms'
        write(6,*) 
        nblock = natoms / ncol
        ip=1
        if (nblock*ncol .eq. natoms) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = natoms - nblock * ncol
          write(6,*) 
          write(6,10) ((ib-1)*ncol+ii,ii=1,imax)
          write(6,*) 
          do i = 1,natoms
            do ii = 1,imax
              p(ii) = 0.0
            enddo
            do ior = lasto(i-1)+1,lasto(i)
              do in = 1,numh(ior)
                jor = listh(in,ior)
                ja = 0
                do jja = 1,natoms
                  if (lasto(jja) .ge. jor) then
                    ja = jja
                    goto 100
                  endif
                enddo
                if (ja .eq. 0) stop 'mulliken: ERROR:  ja = 0'
100             ii = ja - (ib - 1) * ncol
                if (ii .ge. 1 .and. ii .le. imax) then
                  do is = 1,min(nspin,2)
                    p(ii) = p(ii) + dm(in,ior,is) * s(in,ior)
                  enddo
                endif
              enddo
            enddo
            write(6,11) i,(p(ii),ii=1,imax)
          enddo
        enddo
      endif
C ....................

C Compute and print Mulliken Orbital and Atomic Populations ..........
      if (iopt .ge. 1) then
        write(6,*) 
        write(6,"(a)")'mulliken: Atomic and Orbital Populations:'
        if (nspin .le. 2) then
          do is = 1,nspin
            if (nspin .eq. 2) then
              if(is .eq. 1) write(6,'(/,a)') 'mulliken: Spin UP '
              if(is .eq. 2) write(6,'(/,a)') 'mulliken: Spin DOWN '
            endif
            qtot = 0.0 
            do ispec =1, ns  

             atm_label=labelfis(ispec)
             write(6,'(/2a)')'Specie: ', atm_label 
             write(6,'(a4,a7,a6)') 'Atom', 'Qatom', 'Qorb'
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
                 enddo 
                  write(6,'(14x,8a8))')
     .           (sym_label2(io),io=1,8)
              enddo 
              do io=1,nres 
                 nao=nao+1 
                  sym_label2(io)=symfio(ispec,nao) 
               enddo 
               write(6,'(14x,8a8)')
     .         (sym_label2(io),io=1,nres) 


             do ia = 1,natoms
              if(isa(ia).eq.ispec) then
C             Compute charge in each orbital of atom ia
              qa = 0.0
              nao = 0
              do io = lasto(ia-1)+1,lasto(ia)
                nao = nao + 1
                qo(nao) = 0.0
                do in = 1,numh(io)
                  qo(nao) = qo(nao) + dm(in,io,is) * s(in,io)
                enddo
                qa = qa + qo(nao)
              enddo
              qtot = qtot + qa
              write(6,'(i4,f7.3,8f8.3,(/11x,8f8.3))')
     .          ia, qa, (qo(io),io=1,nao) 
              endif 
             enddo 
            enddo
            write(6,"(a,f8.3)") 'mulliken: Qtot = ', qtot
          enddo
        elseif (nspin .eq. 4) then
          call chkdim( 'mulliken', 'normax', normax, nbasis, 1 )
          do is = 1,nspin
            qts(is) = 0.d0
            do io = 1,nbasis
              qos(is,io) = 0.d0
            enddo
          enddo
          do is = 1,nspin
            do io = 1,nbasis
              do in = 1,numh(io)
                jo = listh(in,io)
                qos(is,io) = qos(is,io) + dm(in,io,is)*s(in,io)/2
                qos(is,jo) = qos(is,jo) + dm(in,io,is)*s(in,io)/2
              enddo
            enddo
          enddo 
          do ispec=1,ns
            atm_label=labelfis(ispec)
            write(6,'(/2a)')'Specie: ', atm_label
            write(6,'(/,a4,a9,4x,2a10,3x,a8,/,64(1h-))')
     .      'Atom', 'Orb', 'Charge', 'Spin', 'Svec'

           do ia = 1,natoms 
            if(isa(ia).eq.ispec) then 
            do is = 1,nspin
              qas(is) = 0.d0
            enddo
            do io = lasto(ia-1)+1,lasto(ia)  
C DSP, Writing symmetries for each orbital.
C DSP, Orbitals with a 'P' belong to the polarization shell

              sym_label(1)=symfio(ispec,iphorb(io))
              do is = 1,nspin
                qas(is) = qas(is) + qos(is,io)
                qts(is) = qts(is) + qos(is,io)
              enddo
              call spnvec( nspin, qos(1,io), qtot, stot, svec )
              write(6,'(i4,i5,a8,2f10.5,3x,3f8.3)')
     .          ia, io, sym_label(1), qtot, stot, svec
            enddo
            call spnvec( nspin, qas, qtot, stot, svec )
            write(6,'(i4,4x,a6,3x,2f10.5,3x,3f8.3,/)')
     .        ia, 'Total', qtot, stot, svec 
            endif
           enddo
          enddo 
          call spnvec( nspin, qts, qtot, stot, svec )
          write(6,'(64(1h-),/,2a8,x,2f10.5,3x,3f8.3,/)')
     .      'Total', 'Total', qtot, stot, svec
        endif
      endif
C ...................

10    format(12x,20(2x,i4,2x))
11    format(i12,20(1x,f7.3))
12    format('i = ',i4,'   q = ',f6.3,'   q_orb = ',9f6.3)
13    format(32x,9f6.3) 
15    format(i4,1x,a7,20(1x,f7.3))
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

      implicit          none
      integer           ns
      double precision  qs(ns), qt, st, sv(3)
      double precision  cosph, costh, sinph, sinth, tiny
      parameter ( tiny = 1.d-12 )

      if (ns .eq. 1) then
        qt = qs(1)
        st = 0.d0
        sv(1) = 0.d0
        sv(2) = 0.d0
        sv(3) = 0.d0
      elseif (ns .eq. 2) then
        qt = qs(1) + qs(2)
        st = qs(1) - qs(2)
        sv(1) = 0.d0
        sv(2) = 0.d0
        sv(3) = st
      elseif (ns .eq. 4) then
        qt = qs(1) + qs(2)
        st = sqrt( (qs(1)-qs(2))**2 + 4.d0*(qs(3)**2+qs(4)**2) )
        costh = ( qs(1) - qs(2) ) / ( st + tiny )
        sinth = sqrt( 1.d0 - costh**2 )
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


