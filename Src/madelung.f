c $Id: madelung.f,v 1.2 1999/05/05 17:25:32 emilio Exp $

      subroutine madelung(cell, shape, charnet, Emad)

c *******************************************************************
c Finds Madelung correction for charged systems
c (Makov & Payne, PRB 51, 4014 (1995))
c (Values of Madelung constant: 
C  Coldwell-Horsfall & Maradudin, J. Math. Phys. 1, 395 (1960))
c
c Written by P. Ordejon, March 1999.
c ********* INPUT ***************************************************
c real*8       cell(3,3): Lattice vectors (in Bohr)
c character*10 shape    : Cell shape (atom,molecule,chain,slab,bulk)
c real*8       charnet  : Net charge of the system
c ********* OUTPUT **************************************************
c real*8       Emad     : Madelung correction energy (in Ry)
c *******************************************************************
      implicit none

      double precision
     .  cell(3,3), charnet, Emad

      character
     .  shape*10

C ÑLocal variables.........

      double precision
     .  alpha, lv
 
      character
     .  ctype*4


C If system is charged, find if energy correction terms can be applied
C (cluster or molecule, with SC, FCC or BCC cell) .....................
      if (shape .ne. 'molecule' .and. shape .ne. 'atom') then
        write(6,'(2a)')
     .  'madelung: WARNING: Charged system, but not'
     .  ,' an atom or molecule.'
        write(6,'(2a)')
     .  'madelung: WARNING: Energy correction terms'
     .  ,' can not be applied. '
        write(6,'(2a)')
     .  'madelung: WARNING: Continue only if you really know'
     .  ,' what you are doing.'
        Emad = 0.0d0
        return
      else
        call typecell(cell,ctype,lv)
        if (ctype .eq. 'none') then
          write(6,'(2a)')
     .    'madelung: WARNING: Charged system, but not'
     .    ,' SC, FCC or BCC cell.'
          write(6,'(2a)')
     .    'madelung: WARNING: Energy correction terms'
     .    ,' can not be applied. '
          write(6,'(2a)')
     .    'madelung: WARNING: Continue only if you really know'
     .    ,' what you are doing.'
          Emad = 0.0d0
          return
        else if (ctype .eq. 'sc') then
          alpha = 5.6745947
          write(6,'(a)')
     .    'madelung: Charged system with SC cell'
        else if (ctype .eq. 'bcc') then
          alpha = 7.278473
          write(6,'(a)')
     .    'madelung: Charged system with BCC cell'
        else if (ctype .eq. 'fcc') then
          alpha = 9.1697536
          write(6,'(a)')
     .    'madelung: Charged system with FCC cell'
        else
          write(6,'(a)')
     .    'madelung: ERROR: Wrong type of cell'
          stop
        endif
      endif

      Emad = charnet**2 * alpha / (2.0d0 * lv)
C ...................
      return
      end

