      subroutine rhooda( no, indxuo, np, Datm, rhoatm, 
     .                   iaorb, iphorb, isa )
C ********************************************************************
C Finds the Harris density at the mesh points from the atomic
C occupations.
C Written by P.Ordejon and J.M.Soler. May'95.
C Inverted so that grid points are the outer loop, J.D. Gale, Jan'99
C *********************** InpUT **************************************
C integer no              : Number of basis orbitals
C integer indxuo(no)      : Index of equivalent atom in unit cell
C integer np              : Number of mesh points
C real*8  Datm(no)        : Occupations of basis orbitals in free atom
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C *********************** OUTPUT **************************************
C real    rhoatm(nsp,np)  : Harris (sum of atoms) density at mesh points
C *********************************************************************
C
C  Modules
C
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshphi
C
      implicit none

      integer          no, np
      integer          indxuo(no), 
     .                 iaorb(*), iphorb(*), isa(*)
      real             rhoatm(nsp,np)
      real*8           Datm(no), phip

      integer          i, in, ip, isp, iu, kn, iop, is, iphi, ia, ix
      real*8           Ci, gradCi(3), r2o, r2sp, dxsp(3) 

C  Loop on mesh points
      do ip = 1,np

C  Initialise rhoatm
        do isp = 1, nsp
          rhoatm(isp,ip) = 0.0
        enddo

C  Loop on orbitals of mesh point
        do kn = 1+endpht(ip-1), endpht(ip)
          i = lstpht(kn)
          iu = indxuo(i)

          if (DirectPhi) then
C  Generate phi value and loop on subpoints
            iphi = iphorb(i)
            ia = iaorb(i)
            is = isa(ia)
            r2o = rcut(is,iphi)**2
            iop = listp2(kn)
            do isp = 1,nsp
              do ix = 1,3
                dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
              enddo
              r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
              if (r2sp.lt.r2o) then
                call phiatm(is,iphi,dxsp,phip,gradCi)
                Ci = phip
                rhoatm(isp,ip) = rhoatm(isp,ip) + Datm(iu) * Ci * Ci
              endif
            enddo
          else
C  Loop on sub-points
            do isp = 1,nsp
              Ci = phi(isp,kn)
              rhoatm(isp,ip) = rhoatm(isp,ip) + Datm(iu) * Ci * Ci
            enddo
          endif

        enddo

      enddo

      end
