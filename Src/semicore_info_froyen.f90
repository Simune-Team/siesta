module m_semicore_info_froyen

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)
  
  public :: get_n_semicore_shells

CONTAINS
subroutine get_n_semicore_shells(p,nsemic)
  use m_ncps, only: pseudopotential_t => froyen_ps_t
  
    type(pseudopotential_t), intent(in) :: p
    integer, intent(out)       :: nsemic(0:3)

  ! Returns an array with the number of semicore shells
  ! on channel l (0..3)
  ! It uses the information in Froyen-style pseudopotential files

  ! If no pseudo is available for a given l, a 0 is returned.
  ! If for some reason the generation step did not pseudize
  ! a proper valence state (e.g. Cu 3d), the number of semicore
  ! shells is returned as -1. This should be checked by the caller.
!
!

  integer   :: lmax, inp_lun, l, n, z
  integer    :: nval_gs(0:3)

  character(len=2), allocatable :: orb_arr(:)
  real(dp), allocatable         :: zdown_arr(:)
  real(dp), allocatable         :: zup_arr(:)
  real(dp), allocatable         :: rc_arr(:)
  integer,  allocatable         :: gen_n(:)

  real(dp) :: chgvps

  z = atomic_number(p%name)
  call cnfig(z,nval_gs)

  lmax = p%npotd-1
  allocate(orb_arr(0:lmax))
  allocate(zdown_arr(0:lmax))
  allocate(zup_arr(0:lmax))
  allocate(rc_arr(0:lmax))
  allocate(gen_n(0:lmax))

  call get_ps_conf(p%irel,lmax,p%text,chgvps, &
                   orb_arr,zdown_arr,zup_arr,rc_arr)

  nsemic(:) = 0
  do l = 0, lmax
     read(orb_arr(l),"(i1)") gen_n(l)
     if (gen_n(l) < nval_gs(l)) then
        nsemic(l) =  nval_gs(l) - gen_n(l)
     endif
  enddo
  if (sum(nsemic) > 0) then
     print *, " -- Semicore analysis"
     do l = 0, lmax
        if (gen_n(l) < nval_gs(l)) then
           print "(a,i1,a,2i2)", "For l=",l, ", gen_n, nval_gs:", gen_n(l), nval_gs(l)
           nsemic(l) =  nval_gs(l) - gen_n(l)
           print "(a,i2,1x,a)", "* There are ", nsemic(l), " semicore shells"
        endif
     enddo
  endif

end subroutine get_n_semicore_shells

subroutine get_ps_conf(irel,lmax,text,chgvps, &
                       orb_arr,zdown_arr,zup_arr,rc_arr)
!
!     Attempt to decode the valence configuration used for
!     the generation of the pseudopotential
!     (At least, the valence charge)

      character(len=3), intent(in)  :: irel
      integer, intent(in)           :: lmax
      character(len=70), intent(in) :: text
      real(dp), intent(out)         :: chgvps
      character(len=2), intent(out) :: orb_arr(0:)
      real(dp), intent(out)         :: zdown_arr(0:)
      real(dp), intent(out)         :: zup_arr(0:)
      real(dp), intent(out)         :: rc_arr(0:)

      integer  :: l, itext
      real(dp) :: ztot, zup, zdown, rc_read
      character(len=2) :: orb

      chgvps=0.0_dp

            if(irel.eq.'isp') then
               write(6,'(/,2a)') &
               'Pseudopotential generated from an ', &
               'atomic spin-polarized calculation'

               write(6,'(/,a)') 'Valence configuration '// &
                      'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8080) &
                      orb, zdown, zup, rc_read
 8080             format(a2,f4.2,1x,f4.2,1x,f4.2)

                  orb_arr(l) = orb
                  zdown_arr(l) = zdown
                  zup_arr(l) = zup
                  rc_arr(l) = rc_read

                  chgvps = chgvps + zdown + zup
                  write(6,8085) orb, zdown, zup, rc_read
 8085             format(a2,'(',f4.2,',',f4.2,') rc: ',f4.2)
               enddo

            else
               if(irel.eq.'rel') then
                  write(6,'(/,2a)')  &
               'Pseudopotential generated from a ', &
                      'relativistic atomic calculation'
                  write(6,'(2a)')   &
               'There are spin-orbit pseudopotentials available'
                  write(6,'(2a)')  &
               'Spin-orbit interaction is not included in ', &
                      'this calculation'
               endif
 
               write(6,'(/,a)') 'Valence configuration '// &
                      'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8090) &
                      orb, ztot, rc_read
 8090             format(a2,f5.2,4x,f5.2)

                  orb_arr(l) = orb
                  zdown_arr(l) = ztot
                  zup_arr(l) = 0.0_dp
                  rc_arr(l) = rc_read

                  chgvps = chgvps + ztot
                  write(6,8095) orb, ztot, rc_read
 8095             format(a2,'(',f5.2,') rc: ',f4.2)
               enddo

           endif
           return

 5000    continue       ! Error return: set chgvps to zero
         call die("Error in get_ps_conf")

         end subroutine get_ps_conf

!
! This routine encodes some choices regarding the core-valence split,
! which might not be universal.
! 
      SUBROUTINE CNFIG( Z, NVAL ) 
! Returns the valence configuration for atomic ground state, i.e.
! the principal quantum number NVAL of the valence orbilas for each L
! Originally written by A.R.Williams. Modified by J.M.Soler

      integer,intent(in) :: Z        ! Atomic number
      integer,intent(out):: NVAL(0:3) ! Valence electrons for each L

      integer, parameter :: LMAX=3, NCHNG=15
      integer :: ICHNG, L, LCHNG(NCHNG), ZCHNG(NCHNG)

      ! Originally: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
!!     DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/

      ! Changed to: s valence orbital switched for full p occupation
      !           Li,Na,Na, K, K,Ga,Rb,Rb,In,Cs,Cs,Hf,Tl,Fr,Fr
      DATA ZCHNG / 3,11,11,19,19,31,37,37,49,55,55,72,81,87,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
      DO L=0,LMAX
         NVAL(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         NVAL(L)=NVAL(L)+1
      END DO

      END subroutine cnfig

!
      FUNCTION atomic_number(SYMBOL) result(z)

! Given the atomic symbol, it returns the atomic number
! Based on code by J. Soler

      character(len=2), intent(in)    :: SYMBOL  ! Atomic symbol
      integer                         :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) =  &
               (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
                 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
                 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
                 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                 'Md','No','Lr'/)

     do z = 1, NZ
        if (SYMBOL == NAME(Z)) then
           RETURN
        endif
     enddo
     call die("Cannot find atomic number for " // symbol)
        
   end FUNCTION atomic_number

end module m_semicore_info_froyen

