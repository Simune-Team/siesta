      subroutine mulliken(iopt,natoms,nbasis,nhmax,numh,listh,s,
     .   dm,lasto)
C ********************************************************************
C Subroutine to perform Mulliken population analisys.
C (Overlap and total populations, both for orbitals and for atoms)
C The density matrix (d.m.) and overlap matrix are passed in sparse form
C (both with the same sparse structure)
C There is no output. The populations are printed to the output.
C
C Written by P.Ordejon, October'96
C ************************** INPUT ************************************
C integer iopt                : Work option: 1 = atomic and orbital charges
C                                            2 = 1 + atomic overlap pop.
C                                            3 = 2 + orbital overlap pop.
C integer natoms              : Number of atoms
C integer nbasis              : Number of basis orbitals
C integer nhmax               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements of each row
C integer numh(nbasis)        : First Control vector of d.m. and overlap
C integer listh(nhmax,nbasis) : Second Control vector of d.m. and overlap
C real*8 s(nhmax,nbasis)      : Overlap matrix in sparse form
C real*8 dm(nhmax,nbasis)     : Density matrix in sparse form
C integer lasto(0:maxa)       : Index of last orbital of each atom
C                               (lasto(0) = 0)
C ************************* OUTPUT *************************************
C No output. The results are printed to standard output
C **********************************************************************
      implicit none
      integer
     .   iopt,natoms,nbasis,nhmax

      integer
     .   numh(nbasis),lasto(0:natoms),listh(nhmax,nbasis)

      double precision
     .   dm(nhmax,nbasis),s(nhmax,nbasis)

C Internal parameters ..................................................
C Number of culumns in printout.  Must be smaller than 20
      integer ncol
      parameter (ncol = 10)

C Maximum number of orbitals per atom
      integer normax
      parameter (normax = 60)

      integer i,ia,ib,ii,imax,in,ior,iior,ip,j,ja,jja,jor,nblock

      double precision
     .  ap(ncol),p(ncol),qo(normax),qa,qtot
C ......................

      if ((iopt .ne. 1) .and. (iopt .ne. 2) .and. (iopt .ne. 3)) then
        write(6,"(a)")'mulliken: Wrong iopt'
        stop
      endif
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
          write(6,*) 
          write(6,10) ((ib-1)*ncol+ii,ii=1,imax)
          write(6,*) 
          do i = 1,nbasis
            do ii = 1,imax
              p(ii) = 0.0
            enddo
            do in = 1,numh(i)
              j = listh(in,i)
              ii = j - (ib - 1) * ncol
              if (ii .ge. 1 .and. ii .le. imax) then
                p(ii) = dm(in,i) * s(in,i)
              endif
            enddo
            write(6,11) i,(p(ii),ii=1,imax)
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
                do jja = 1,natoms
                  if (lasto(jja) .ge. jor) then
                    ja = jja
                    goto 100
                  endif
                enddo
100             ii = ja - (ib - 1) * ncol
                if (ii .ge. 1 .and. ii .le. imax) then
                  p(ii) = p(ii) + dm(in,ior) * s(in,ior)
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
        write(6,"(a)")'mulliken: Atomic and Orbital Populations'
        write(6,*) 
	qtot = 0.0
        do ia = 1,natoms
C Compute charge in each orbital of atom ia
          qa = 0.0
          ior = 0
          do iior = lasto(ia-1)+1,lasto(ia)
            ior = ior + 1
            qo(ior) = 0.0
            do in = 1,numh(iior)
              qo(ior) = qo(ior) + dm(in,iior) * s(in,iior)
            enddo
            qa = qa + qo(ior)
          enddo
	  qtot = qtot + qa
          if (ior .le. 9) then
            write(6,12) ia,qa,(qo(i), i=1,ior)
          else
            write(6,12) ia,qa,(qo(i), i=1,9)
            nblock = ior/9
            ip=1
            if (nblock*9 .eq. ior) ip=0
            do j = 2,nblock+ip
              if (j .ne. nblock+ip) then
                write(6,13) (qo(i), i=(j-1)*9+1,j*9)
              else
                write(6,13) (qo(i), i=(j-1)*9+1,ior)
              endif
            enddo
          endif
        enddo
	write(6,"(a,f8.3)")'mulliken: Qtot = ', qtot
      endif
C ...................
10    format(5x,20(3x,i4))
11    format(i5,20f7.3)
12    format('i = ',i4,'   q = ',f6.3,'   q_orb = ',9f6.3)
13    format(32x,9f6.3)


      return
      end

