  subroutine register_rfs()

    use m_radfunc_registry, only: register_in_rf_pool
    use atm_types, only: species_info, species, nspecies
    use radial, only: rad_func

    implicit none 
!
    type(species_info), pointer        :: spp
    type(rad_func), pointer            :: func   
!
    integer :: is, io, l, m, gindex

    do is = 1, nspecies
       spp => species(is)
       do io=1,spp%norbs
          func => spp%orbnl(spp%orb_index(io))
          l = spp%orb_l(io)
          m = spp%orb_m(io)
          call register_in_rf_pool(func,l,m,"orb",(/is,io/),gindex)
          spp%orb_gindex(io) = gindex
       enddo
    enddo

    ! KB projectors
    do is = 1, nspecies
       spp => species(is)
       do io=1,spp%nprojs
          func => spp%pjnl(spp%pj_index(io))
          l = spp%pj_l(io)
          m = spp%pj_m(io)
          call register_in_rf_pool(func,l,m,"kbproj",(/is,io/),gindex)
          spp%pj_gindex(io) = gindex
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
    enddo
    
  end subroutine register_rfs
!
!   Test
!
  subroutine test_register()

    use precision, only: dp
    use m_radfunc_registry, only: cutoff=>rcut, evaluate
    use atm_types, only: species_info, species, nspecies
    use atmfuncs, only: rcut, phiatm
    use atmfuncs, only: orb_gindex, kbproj_gindex, vna_gindex

    implicit none 
!
    type(species_info), pointer        :: spp
!
    integer :: is, io, gindex
    real(dp) :: r(3) = (/0.5_dp, 0.5_dp, 0.5_dp/)
    real(dp) :: grad(3), phi

    do is = 1, nspecies
       print *, "---IS: ", is
       spp => species(is)
       do io=1,spp%norbs
          call phiatm(is,io,r,phi,grad)
          print *, "io, rcut, phiatm_h:",   &
                    io, rcut(is,io), phi
          gindex = orb_gindex(is,io)
          call evaluate(gindex,r,phi,grad)
          print *, "ig, rcut, phiatm_h:",   &
                    gindex, cutoff(gindex), phi
       enddo
    enddo

    do is = 1, nspecies
       print *, "---IS projs: ", is
       spp => species(is)
       do io=1,spp%nprojs
          call phiatm(is,-io,r,phi,grad)
          print *, "io, rcut, phiatm_h:",   &
                    io, rcut(is,-io), phi
          gindex = kbproj_gindex(is,io)
          call evaluate(gindex,r,phi,grad)
          print *, "ig, rcut, phiatm_h:",   &
                    gindex, cutoff(gindex), phi
       enddo
    enddo
    
    ! Vna
    do is = 1, nspecies
       print *, "---IS vna: ", is
       spp => species(is)
       call phiatm(is,0,r,phi,grad)
       print *, "rcut, phiatm_h:",   &
                 rcut(is,0), phi
       gindex = vna_gindex(is)
       call evaluate(gindex,r,phi,grad)
       print *, "ig, rcut, phiatm_h:",   &
                    gindex, cutoff(gindex), phi
    enddo
    
  end subroutine test_register
