      subroutine rhooda( NO, indxuo, endCt, listCt, CtoCt, 
     .                   C, maxC, NSP, NP, Datm, RHOatm, 
     .                   iaorb, iphorb, isa, listP2, dxa, 
     .                   xdop, xdsp, DirectPhi )
C ********************************************************************
C Finds the Harris density at the mesh points from the atomic
C occupations.
C Written by P.Ordejon and J.M.Soler. May'95.
C Inverted so that grid points are the outer loop, J.D. Gale, Jan'99
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer indxuo(NO)      : Index of equivalent atom in unit cell
C integer endCt(NO)       : Accumulated umber of nonzero elements in
C                           each row of C
C integer listCt(*)       : List of nonzero elements in each row of C
C integer CtoCt(*)        : Maps C to transpose of C
C real*4  C(NSP,maxC)     : Values of nonzero elements in C
C integer maxC            : second dimension of C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of mesh points
C real*8  Datm(NO)        : Occupations of basis orbitals in free atom
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C integer listP2(*)       : Mesh point pointer for orbital
C real*8  dxa(3,*)        : Atom position within mesh-cell
C real*8  xdop(3,*)       : Vector to mesh points within rmax
C real*8  xdsp(3,*)       : Vector to mesh sub-points
C logical DirectPhi       : if .true. then C is not to be used and phi
C                         : values must be calculated on the fly
C *********************** OUTPUT **************************************
C real*4  RHOatm(NSP,NP)  : Harris (sum of atoms) density at mesh points
C *********************************************************************
C
C  Modules
C
      use atmfuncs, only: rcut, phiatm
C
      implicit none

      integer          NO, NP, NSP, maxC
      integer          endCt(0:NP), indxuo(NO), listCt(*),
     .                 CtoCt(*), iaorb(*), iphorb(*), isa(*),
     .                 listP2(*)
      real             C(NSP,maxC),RHOatm(NSP,NP)
      double precision Datm(NO), dxa(3,*), xdop(3,*), xdsp(3,*), phip
      logical          DirectPhi

      integer          i, in, ip, isp, iu, kn, iop, is, iphi, ia, ix
      double precision Ci, gradCi(3), r2o, r2sp, dxsp(3) 

C  Loop on mesh points
      do ip = 1,np

C  Initialise RHOatm
        do isp = 1, nsp
          RHOatm(isp,ip) = 0.d0
        enddo

C  Loop on orbitals of mesh point
        do kn = 1+endCt(ip-1), endCt(ip)
          i = listCt(kn)
          in = CtoCt(kn)
          iu = indxuo(i)

          if (DirectPhi) then
C  Generate phi value and loop on subpoints
            iphi = iphorb(i)
            ia = iaorb(i)
            is = isa(ia)
            r2o = rcut(is,iphi)**2
            iop = listP2(in)
            do isp = 1,nsp
              do ix = 1,3
                dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
              enddo
              r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
              if (r2sp.lt.r2o) then
                call phiatm(is,iphi,dxsp,phip,gradCi)
                Ci = phip
                RHOatm(isp,ip) = RHOatm(isp,ip) + Datm(iu) * Ci * Ci
              endif
            enddo
          else
C  Loop on sub-points
            do isp = 1,nsp
              Ci = C(isp,in)
              RHOatm(isp,ip) = RHOatm(isp,ip) + Datm(iu) * Ci * Ci
            enddo
          endif

        enddo

      enddo

      end
