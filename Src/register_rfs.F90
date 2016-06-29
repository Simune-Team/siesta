! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
  subroutine register_rfs()

!
!   Installs a record of the PAO, KB projector, and Vna
!   radial functions in the global registry used by the
!   new version of matel.
!   The registration returns a global index which should be
!   kept for the invokation of matel.
!
    use m_matel_registry, only: register_in_rf_pool, show_pool
    use atm_types, only: species_info, species, nspecies
    use radial, only: rad_func

    implicit none 
!
    type(species_info), pointer        :: spp
    type(rad_func), pointer            :: func   
!
    integer :: is, io, ko, l, m, gindex

    do is = 1, nspecies
       spp => species(is)
       do io=1,spp%norbs
          func => spp%orbnl(spp%orb_index(io))
          l = spp%orb_l(io)
          m = spp%orb_m(io)
          call register_in_rf_pool(func,l,m,"orb",(/is,io/),gindex)
          spp%orb_gindex(io) = gindex
!!         For debugging
!          write(6,*)'Atomic orbitals'
!          write(6,*)'is, io, gindex = ', is, io, gindex 
!!         End debugging
       enddo
    enddo

    ! KB projectors
    do is = 1, nspecies
       spp => species(is)
       do ko=1,spp%nprojs
          func => spp%pjnl(spp%pj_index(ko))
          l = spp%pj_l(ko)
          m = spp%pj_m(ko)
          io = -ko
          call register_in_rf_pool(func,l,m,"kbproj",(/is,io/),gindex)
          spp%pj_gindex(ko) = gindex
!!         For debugging
!          write(6,*)'KB projectors'
!          write(6,*)'is, ko, gindex = ', is, ko, gindex 
!!         End debugging
       enddo
    enddo
    
    ! Vna
    do is = 1, nspecies
       spp => species(is)
       func => spp%vna
       l = 0
       m = 0
       call register_in_rf_pool(func,l,m,"vna",(/is/),gindex)
       spp%vna_gindex = gindex
!!         For debugging
!          write(6,*)'VNA'
!          write(6,*)'is, gindex = ', is, gindex 
!!         End debugging
    enddo

!!    call show_pool()
    
  end subroutine register_rfs
!
!   Test
!
!!$  subroutine test_register()
!!$
!!$    use precision, only: dp
!!$    use m_matel_registry, only: cutoff=>rcut, evaluate
!!$    use atm_types, only: species_info, species, nspecies
!!$    use atmfuncs, only: rcut, phiatm
!!$    use atmfuncs, only: orb_gindex, kbproj_gindex, vna_gindex
!!$
!!$    implicit none 
!!$!
!!$    type(species_info), pointer        :: spp
!!$!
!!$    integer :: pass, is, io, ko, gindex
!!$    real(dp) :: r(3) = (/0.5_dp, 0.5_dp, 0.5_dp/)
!!$    real(dp) :: grad(3), phi
!!$    logical, external :: io_node
!!$
!!$   do pass = 0, 1
!!$    if (pass==1) r(:) = 0.0_dp
!!$    if (io_node()) print *, "Evaluation for r:",r(:)
!!$    do is = 1, nspecies
!!$       if (io_node()) print *, "---IS: ", is
!!$       spp => species(is)
!!$       do io=1,spp%norbs
!!$          call phiatm(is,io,r,phi,grad)
!!$          if (io_node()) print "(a,i3,f12.6,g20.10)", "io, rcut, phiatm_h:",   &
!!$                    io, rcut(is,io), phi
!!$          gindex = orb_gindex(is,io)
!!$          call evaluate(gindex,r,phi,grad)
!!$          if (io_node()) print "(a,i3,f12.6,g20.10)", "ig, rcut, phiatm_h:",   &
!!$                    gindex, cutoff(gindex), phi
!!$       enddo
!!$    enddo
!!$
!!$    do is = 1, nspecies
!!$       if (io_node()) print *, "---IS projs: ", is
!!$       spp => species(is)
!!$       do io=1,spp%nprojs
!!$          call phiatm(is,-io,r,phi,grad)
!!$          if (io_node()) print "(a,i3,f12.6,g20.10)", "io, rcut, phiatm_h:",   &
!!$                    io, rcut(is,-io), phi
!!$          gindex = kbproj_gindex(is,-io)
!!$          call evaluate(gindex,r,phi,grad)
!!$          if (io_node()) print "(a,i3,f12.6,g20.10)", "ig, rcut, phiatm_h:",   &
!!$                    gindex, cutoff(gindex), phi
!!$       enddo
!!$    enddo
!!$    
!!$    ! Vna
!!$    do is = 1, nspecies
!!$       if (io_node()) print *, "---IS vna: ", is
!!$       spp => species(is)
!!$       call phiatm(is,0,r,phi,grad)
!!$       if (io_node()) print "(a,f12.6,g20.10)", "rcut, phiatm_h:",   &
!!$                 rcut(is,0), phi
!!$       gindex = vna_gindex(is)
!!$       call evaluate(gindex,r,phi,grad)
!!$       if (io_node()) print "(a,i3,f12.6,g20.10)", "ig, rcut, phiatm_h:",   &
!!$                    gindex, cutoff(gindex), phi
!!$    enddo
!!$ enddo
!!$    
!!$  end subroutine test_register
