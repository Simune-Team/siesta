      module chemical

      use sys

      implicit none

      private
      public number_of_species, atomic_number, species_label
      public is_floating, is_bessel
      public read_chemical_types
!
!     Species information
!
      type chemical_types
         integer                       :: no_of_species
         character(len=20), pointer    :: spec_label(:)
         integer, pointer              :: z(:)
      end type chemical_types

      type(chemical_types), save       :: chemical_list


      CONTAINS

      subroutine check(i)
      integer, intent(in) :: i
      if (i.lt.0 .or. i.gt.chemical_list%no_of_species)
     $     call die("Wrong species number requested")
      end subroutine check

      function number_of_species()
      integer number_of_species
      number_of_species = chemical_list%no_of_species
      end function number_of_species

      function species_label(i)
      character(len=20) species_label
      integer, intent(in)  :: i

      call check(i)
      species_label = chemical_list%spec_label(i)
      end function species_label

      function atomic_number(i)
      integer atomic_number
      integer, intent(in)  :: i

      call check(i)
      atomic_number = chemical_list%z(i)
      end function atomic_number

      function is_floating(i)
      logical is_floating
      integer, intent(in)  :: i

      call check(i)
      is_floating = (chemical_list%z(i) .le. 0)
      end function is_floating

      function is_bessel(i)
      logical is_bessel
      integer, intent(in)  :: i

      call check(i)
      is_bessel = (chemical_list%z(i) .eq. -100)
      end function is_bessel

!---
      subroutine read_chemical_types

      use parse
      use fdf

      integer nsp, isp
      integer ns_read

      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      character(len=20) label
      integer z
      logical floating, bessel

      nsp = fdf_get('Number_of_species',0)
      if (nsp.eq.0) call die("No species found!!!")

      allocate(chemical_list%spec_label(nsp))
      allocate(chemical_list%z(nsp))
      chemical_list%no_of_species = nsp
     
      nullify(bp)
      if (.not. fdf_block('Chemical_species_label',bp) )
     $     call die("Block Chemical_species_label does not exist.")

      ns_read = 0
      loop: DO
        if (.not. fdf_bline(bp,line)) exit loop
        ns_read = ns_read + 1
        p => digest(line)
        if (.not. match(p,"iin"))
     $       call die("Wrong format in Chemical_species_label")
        isp = integers(p,1)
        if (isp .gt. nsp .or. isp .lt. 1)
     $       call die("Wrong specnum in Chemical_species_label")
        label = names(p,1)
        z = integers(p,2)
        floating = (z.le.0)
        bessel = (z.eq.-100)

C        if (.not.floating) then
C         write(6,*) 'Species number: ', isp,
C     $             ' Label: ', trim(label),
C     $             ' Atomic number:',  z
C        elseif (bessel) then
C         write(6,*) 'Species number: ', isp,
C     $             ' Label: ', trim(label),
C     $             ' (floating Bessel functions)'
C        else
C         write(6,*) 'Species number: ', isp,
C     $             ' Label: ', trim(label),
C     $             ' Atomic number:',  z,
C     $             ' (floating PAOs)'
C        endif
!
        chemical_list%z(isp) = z
        chemical_list%spec_label(isp) = label
 
       call destroy(p)
      enddo loop
      if (ns_read .ne. nsp) call die("Not enough species in block")
      call destroy(bp)

      end subroutine read_chemical_types

      end module chemical






