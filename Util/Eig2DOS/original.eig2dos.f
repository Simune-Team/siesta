! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
c $Id: eig2dos.f,v 1.1 1999/04/08 10:20:52 emilio Exp $

      program eig2dos

***********************************************************************
* Utility for obtaining the total density of states out of the .EIG 
* file generated by SIESTA. Energies in eV.
*
* Written by Emilio Artacho, April 1999.
***********************************************************************
* Read from standard input "systemlabel.EIG" generated by SIESTA adding
*   - peak width (in eV) for broadening (gaussian or lorentzian).
*   - number of points in the energy window
*   - energy window: Emin and Emax (wrt E_F shifted to zero)
* in the first line of the file (after the fermi energy present there)
*
* Density of states in standard output (energy shifted, zero at E_F)
*   plus the number of electrons (states) in the energy window, 
*   with and without broadening (number of eigenvalues and 
*   DOS integral, with simple sum)
*
* Equal weight for all k-points is assumed. This condition can be 
* removed by reading the weights from "systemlabel.KP"
***********************************************************************

      implicit none

      integer           maxne, maxb, maxs
      parameter         (maxne = 10000, maxb = 10000, maxs = 2)

      integer           nk, nspin, nband, ne, ie, ik, ika, is, ib, nel
      double precision  e, emin, emax, eincr, ef, eta, pi, x, sum, sta,
     .                  eig(maxb,maxs), dos(maxne), norm

      logical           overflow, loren
      parameter         (loren = .false.)
c ---------------------------------------------------------------------

c reading and initializing --------------------------------------------

      pi = dacos(-1.d0)

      read(5,*) ef, eta, ne, emin, emax
      read(5,*) nband, nspin, nk

      overflow = (ne .gt. maxne) .or. (nband .gt. maxb)
      if (overflow) stop 'Dimensions in eig2dos too small'

      if (ne .lt. 2) ne = 2
      eincr = (emax - emin) / dfloat(ne - 1)

      nel = 0
      do ie = 1, ne
         dos(ie) = 0.d0
      enddo

c big loop: for each eigenvalue a lorentzian is added -----------------

      do ik = 1, nk

         read(5,*) ika, ((eig(ib,is), ib = 1, nband), is = 1, nspin)

         do is = 1, nspin
            do ib = 1, nband
               e = eig(ib,is) - ef
               if ( (e.ge.emin) .and. (e.le.emax) ) nel = nel + 1
               do ie = 1, ne
                  x = emin + (ie-1)*eincr - e
                  if (loren) then
                     dos(ie) = dos(ie) + eta / (eta*eta + x*x) 
                  else
                     dos(ie) = dos(ie) + dexp( - x*x/(eta*eta) )
                  endif
               enddo
            enddo
         enddo
      enddo

      if (loren) then
         norm = pi * dfloat(nk)
      else
         norm = dsqrt(pi) * eta * dfloat(nk)
      endif

      do ie = 1, ne
         dos(ie) = dos(ie)/norm
      enddo

c integral, extremely sophisticated -----------------------------------

      sum = 0.d0
      do ie = 1, ne
         sum = sum + dos(ie)*eincr
      enddo

c number of electrons -------------------------------------------------

      if (nspin .eq. 1) then
         sum = sum * 2.d0
         nel = nel * 2
      endif
      sta = nel / dfloat(nk)

c output, prepared for gnuplot ----------------------------------------

      write(6,"(2a)") '# EIG2DOS: Utility for SIESTA to obtain the ',
     .                'electronic density of states'
      write(6,"(a)") '#'
      write(6,"(2a)") '#                                           ',
     .                '       Emilio Artacho, Apr. 1999'
      write(6,"(2a)") '# ------------------------------------------',
     .                '--------------------------------'
      write(6,"(a,3i6)") '# Nbands, Nspin, Nk   = ', nband, nspin, nk
      write(6,"(a,f10.4,a)") '# E_F                 = ', ef , ' eV'
      write(6,"(a,f10.4,a)") '# Broadening          = ', eta, ' eV'
      write(6,"(a,2f10.4,a)") '# Number of electrons = ', sta, sum,
     .                        ' per cell'
      write(6,"(a)") '#'
      write(6,"(a)") '#        E            N'
      write(6,"(2a)") '# ------------------------------------------',
     .                '--------------------------------'

      write(6,"(2f14.6)") (emin + (ie-1)*eincr, dos(ie), ie = 1, ne)

      stop
      end
