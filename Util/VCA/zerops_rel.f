      program zerops

!
!     Pseudopotential for a "zero" atom.
!
      use precision,       only: dp
      use pseudopotential, only: pseudopotential_t,
     $                           pseudo_write_formatted
      use sys,             only: die

      implicit none

      type(pseudopotential_t) :: p
      integer :: i
      
      p%name = "ZR"
      p%nr = 1029
      p%nrval = p%nr + 1
      p%zval = 0.0_dp
      p%gen_zval = 0.0_dp
      p%icorr = "ca"
      p%irel = "rel"
      p%nicore = "nc"
      p%a = 0.0125_dp
      p%b = 0.309844e-03_dp

      p%method(1) = "ZEROPSEUDO"
      p%method(2) = " 4down-3up"
      p%method(3:6) = ""

      p%text = "Zero pseudo, rel"

      p%npotu = 3
      p%npotd = 4
      allocate(p%r(p%nrval))
      p%r = 1.1_dp
      allocate(p%chcore(p%nrval))
      allocate(p%chval(p%nrval))
      p%chcore = 0.0_dp
      p%chval = 0.0_dp

      allocate(p%vdown(p%npotd,p%nrval))
      allocate(p%ldown(p%npotd))
      do i = 1, p%npotd
         p%vdown(i,:) = 0.0_dp
         p%ldown(i) = i - 1
      enddo

      allocate(p%vup(p%npotu,p%nrval))
      allocate(p%lup(p%npotu))
      do i = 1, p%npotu
         p%vup(i,:) = 0.0_dp
         p%lup(i) = i   ! Note this !!
      enddo

      call pseudo_write_formatted("Zero.rel.psf",p)

      end program zerops


