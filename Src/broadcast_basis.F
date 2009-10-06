! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine broadcast_basis
!
!     Globalizes the basis and potentials data structures
!
!     Alberto Garcia, June 2000--

      use atm_types
      use radial
      use parallel,   only : Node, Nodes

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifndef MPI
!
!     Do nothing...
!
      end subroutine broadcast_basis
#else
      integer MPIerror

      integer is, i
      type(species_info), pointer  :: spp
      type(rad_func), pointer      :: func

      if (Nodes.eq.1) return

      call MPI_Bcast(nspecies,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(npairs,1,MPI_integer,0,MPI_Comm_World,MPIerror)

      if (Node.ne.0) allocate (species(nspecies))

      do is=1,nspecies
         spp => species(is)
         call MPI_Bcast(spp%symbol,2,MPI_character,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%label,20,MPI_character,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%z,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%mass,1,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%zval,1,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%self_energy,1,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)

         call MPI_Bcast(spp%n_orbnl,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%lmax_basis,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orbnl_l,maxn_orbnl,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orbnl_n,maxn_orbnl,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orbnl_z,maxn_orbnl,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orbnl_ispol,maxn_orbnl,MPI_logical,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orbnl_pop,maxn_orbnl,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)



         call MPI_Bcast(spp%n_pjnl,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%lmax_projs,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pjnl_l,maxn_pjnl,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pjnl_n,maxn_pjnl,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pjnl_ekb,maxn_pjnl,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)

         call MPI_Bcast(spp%norbs,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)

         call MPI_Bcast(spp%orb_index,maxnorbs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orb_n,maxnorbs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orb_l,maxnorbs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orb_m,maxnorbs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%orb_pop,maxnorbs,MPI_double_precision,
     $                  0,MPI_Comm_World,MPIerror)

         call MPI_Bcast(spp%nprojs,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)

         call MPI_Bcast(spp%pj_index,maxnprojs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pj_n,maxnprojs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pj_l,maxnprojs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(spp%pj_m,maxnprojs,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)


         if (Node.ne.0) then
            allocate(spp%orbnl(spp%n_orbnl))
            allocate(spp%pjnl(spp%n_pjnl))
         endif
         
         do i=1,spp%n_orbnl
            func=>spp%orbnl(i)
            call broadcast_radial(func,Node)
         enddo

         do i=1,spp%n_pjnl
            func=>spp%pjnl(i)
            call broadcast_radial(func,Node)
         enddo

         func=>spp%vna
         call broadcast_radial(func,Node)

         func=>spp%chlocal
         call broadcast_radial(func,Node)

!!       !! Note that reduced_vlocal might not be initialized
!!       !! if we have read the basis information from a .ion file.
!!       !! reduced_vlocal is not needed by Siesta, anyway, and
!!       !! at this point is has already been written out to
!!       !! any plotting files (or ion.nc files, which for some 
!!       !! historical reason are more complete
!!
!!         func=>spp%reduced_vlocal
!!         call broadcast_radial(func,Node)

         call MPI_Bcast(spp%there_is_core,1,MPI_logical,
     $                  0,MPI_Comm_World,MPIerror)
         
         func=>spp%core
         call broadcast_radial(func,Node)

      enddo
!
!     Electrostatic correction tables
!
      if (Node.ne.0)   allocate (elec_corr(npairs))

      do i = 1, npairs
         func=>elec_corr(i)
         call broadcast_radial(func,Node)
      enddo
         
!
!     Done
!
      CONTAINS

      subroutine broadcast_radial(func,Node)

      type(rad_func), pointer   :: func
      integer, intent(in)       :: Node

      call MPI_Bcast(func%n,1,MPI_integer,
     $                  0,MPI_Comm_World,MPIerror)
      if (func%n .eq. 0) then
         call rad_zero(func)
      else
         if (Node.ne.0) call rad_alloc(func,func%n)
         call MPI_Bcast(func%cutoff,1,MPI_double_precision,
     $        0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(func%delta,1,MPI_double_precision,
     $        0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(func%f,func%n,MPI_double_precision,
     $        0,MPI_Comm_World,MPIerror)
         call MPI_Bcast(func%d2,func%n,MPI_double_precision,
     $        0,MPI_Comm_World,MPIerror)
      endif
      end subroutine broadcast_radial

      end subroutine broadcast_basis


#endif



