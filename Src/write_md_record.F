      subroutine write_md_record( istep )
        USE SIESTA_TODO
        Use siesta_options
        use siesta_geom
        use atomlist, only: iza
        use md_out,     only: md_v_format
        use parallel, only: IOnode
        use m_energies
        use m_steps

        implicit none

        integer, intent(in) ::  istep

        real(dp) :: Etot_output
        logical:: lastst  

        ! Save atomic positions and velocities accumulatively
        if (writmd.and.IOnode) then
          if (.not. harrisfun) then
            Etot_output = Etot
          else
            Etot_output = Eharrs1
          endif
          getot = Etot_output + Ekinion + kn + kpr + vn + vpr
          call iomd( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, istep, inicoor, 
     .        fincoor, tempion, Etot_output, getot,
     .        volume/Ang**3, Psol/kbar)
          call md_v_format(na_u,isa,xa,ucell)
#         ifdef CDF
            call md_netcdf( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, 
     .        tempion, Etot_output, getot,
     .        volume/Ang**3, Psol/kbar)
#         endif

        endif

        ! Accumulate coor in Xmol file for animation 
        lastst = fincoor .le. istep
        if (writpx.and.IOnode) then
          call pixmol(iza, xa, na_u, lastst)
        endif
      END
