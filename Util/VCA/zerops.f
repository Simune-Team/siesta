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
      
      p%name = "ZN"
      p%nr = 1029
      p%nrval = p%nr + 1
      p%zval = 0.0_dp
      p%gen_zval = 0.0_dp
      p%icorr = "ca"
      p%irel = "nrl"
      p%nicore = "nc"
      p%a = 0.0125_dp
      p%b = 0.309844e-03_dp

      p%method(1) = "ZEROPSEUDO"
      p%method(1) = " 4down-0up"
      p%method(3:6) = ""

      p%text = "Zero pseudo, nrl"

      p%npotu = 0
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

      call pseudo_write_formatted("Zero.nrl.psf",p)

      end program zerops


