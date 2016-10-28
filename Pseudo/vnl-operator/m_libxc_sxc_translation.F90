module m_libxc_sxc_translation
!
!  Support for XC functional name translation
!  Version appropriate for ATOM-SiestaXC-Libxc
!
!  Alberto Garcia, June-Sep 2014
!

! Get the symbolic handles from the specific modules
use m_libxc_list
use m_siestaxc_list

implicit none

private

type, public :: xc_id_t
   type(siestaxc_t)   :: siestaxc_id
   type(libxc_t)      :: libxc_id(2)
   character(len=2)   :: atom_id      ! simple two-char code is enough
end type xc_id_t


public :: get_xc_id_from_atom_id, print_xc_id, xc_id_to_string
public :: get_xc_id_from_siestaxc, xc_is_not_lda, xc_nfuncs_libxc
public :: get_xc_id_from_libxc

! These names are placeholders, not currently implemented in libxc
type(libxc_t), parameter :: XC_GGA_C_PBE_XXX = XC_NOT_IMPL
type(libxc_t), parameter :: XC_GGA_X_PBE_XXX = XC_NOT_IMPL
type(libxc_t), parameter :: XC_VDW_C_DF1 = XC_NOT_IMPL
type(libxc_t), parameter :: XC_VDW_C_DF2 = XC_NOT_IMPL
type(libxc_t), parameter :: XC_VDW_C_VV10 = XC_NOT_IMPL
type(libxc_t), parameter :: XC_GGA_X_CX_VDW = XC_NOT_IMPL

type(xc_id_t), dimension(26) :: xct  = (/   &

! LDA ------
xc_id_t(SXC_LDA_PZ, (/XC_LDA_X, XC_LDA_C_PZ/), "ca"),  &
xc_id_t(SXC_LDA_CA, (/XC_LDA_X, XC_LDA_C_PZ/), "ca"),     &  !alias
xc_id_t(SXC_LDA_PW92, (/XC_LDA_X, XC_LDA_C_PW/), "pw"),  &

! not in the current SiestaXC; only in old excorr in atom.
xc_id_t(SXC_LDA_WIGNER, (/XC_LDA_X, XC_LDA_C_WIGNER/), "wi"),  &
xc_id_t(SXC_LDA_HL, (/XC_LDA_X, XC_LDA_C_HL/), "hl"),  &
xc_id_t(SXC_LDA_GL, (/XC_LDA_X, XC_LDA_C_GL/), "gl"),  &
xc_id_t(SXC_LDA_VBH, (/XC_LDA_X, XC_LDA_C_vBH/), "bh"),  &

! GGA ------
xc_id_t(SXC_GGA_PW91, (/XC_GGA_X_PW91, XC_GGA_C_PW91/), "wp"), &
xc_id_t(SXC_GGA_PBE, (/XC_GGA_X_PBE, XC_GGA_C_PBE/), "pb"), &
!"RPBE - Hammer et al"
xc_id_t(SXC_GGA_RPBE, (/XC_GGA_X_RPBE, XC_GGA_C_PBE/), "rp"), &
!"revPBE Zhang+Yang"
xc_id_t(SXC_GGA_revPBE, (/XC_GGA_X_PBE_R, XC_GGA_C_PBE/), "rv"), &
!"Becke-Lee-Yang-Parr"
xc_id_t(SXC_GGA_LYP, (/XC_GGA_X_B88, XC_GGA_C_LYP/), "bl"), &
!"Wu-Cohen"
xc_id_t(SXC_GGA_WC, (/XC_GGA_X_WC, XC_GGA_C_PBE/), "wc"), & ! ???
!"Perdew-Burke-Ernzerhof-solid"
xc_id_t(SXC_GGA_PBEsol, (/XC_GGA_X_PBE_SOL, XC_GGA_C_PBE_SOL/), "ps"), &
!"Armiento-Mattsson-05"
xc_id_t(SXC_GGA_AM05, (/XC_GGA_X_AM05, XC_GGA_C_AM05/), "am"), &

! not yet implemented in libxc
xc_id_t(SXC_GGA_PBEJsJrLO, (/XC_GGA_X_PBE_JSJR, XC_GGA_C_PBE_XXX/), "jo"), & 
xc_id_t(SXC_GGA_PBEJsJrHEG, (/XC_GGA_X_PBE_XXX, XC_GGA_C_PBE_XXX/), "jh"), & 
xc_id_t(SXC_GGA_PBEGcGxLO, (/XC_GGA_X_PBE_XXX, XC_GGA_C_PBE_XXX/), "go"), &
xc_id_t(SXC_GGA_PBEGcGxHEG, (/XC_GGA_X_PBE_XXX, XC_GGA_C_PBE_XXX/), "gh"), &

! VDW -----
xc_id_t(SXC_VDW_DRSLL, (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_DF1/), "vw"), &
xc_id_t(SXC_VDW_DRSLL, (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_DF1/), "vf"), & !alias
xc_id_t(SXC_VDW_LMKLL, (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_DF2/), "vl"), &
xc_id_t(SXC_VDW_KKBM,  (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_DF1/), "vk"), &
xc_id_t(SXC_VDW_C09,   (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_DF1/), "vc"), &
xc_id_t(SXC_VDW_BH,    (/XC_GGA_X_CX_VDW, XC_VDW_C_DF1/), "vb"), &
xc_id_t(SXC_VDW_VV,    (/XC_GGA_X_OPTB88_VDW, XC_VDW_C_VV10/), "vv") &
                                      /)


CONTAINS

  subroutine get_xc_id_from_atom_id(atom_id,xc_id,stat)
    character(len=2), intent(in) :: atom_id
    type(xc_id_t), intent(out)   :: xc_id
    integer, intent(out)         :: stat

    integer :: i

    stat = -1
    do i = 1, size(xct)
       if (xct(i)%atom_id == atom_id) then
          xc_id = xct(i)
          stat = 0
       endif
    enddo
  end subroutine get_xc_id_from_atom_id

  subroutine get_xc_id_from_siestaxc(xc_type,xc_authors,xc_id,stat)
    character(len=*), intent(in) :: xc_type
    character(len=*), intent(in) :: xc_authors
    type(xc_id_t), intent(out)   :: xc_id
    integer, intent(out)         :: stat

    integer :: i

    stat = -1
    do i = 1, size(xct)
       if ((xct(i)%siestaxc_id%family == xc_type)  .and. &
           (xct(i)%siestaxc_id%authors == xc_authors))  then
          xc_id = xct(i)
          stat = 0
       endif
    enddo
  end subroutine get_xc_id_from_siestaxc

  function xc_is_not_lda(xc_id) result (p)
    type(xc_id_t), intent(in)   :: xc_id
    logical :: p

    p = .not. (xc_id%siestaxc_id%family(1:3) == "LDA")

  end function xc_is_not_lda

  function xc_nfuncs_libxc(xc_id) result (n)
    type(xc_id_t), intent(in)   :: xc_id
    integer :: n

    integer :: i

    n = 2   ! Always
!    n = 0
!    do i = 1, 2
!       if (trim(xc_id%libxc_id(i)%name)  == "XC_EMPTY") cycle
!       n = n + 1
!    enddo

  end function xc_nfuncs_libxc

  subroutine get_xc_id_from_libxc(libxc_ids,xc_id,stat)
    integer, intent(in)          :: libxc_ids(2)
    type(xc_id_t), intent(out)   :: xc_id
    integer, intent(out)         :: stat

    integer :: i, c1, c2

    stat = -1
    do i = 1, size(xct)
       c1 = xct(i)%libxc_id(1)%code
       c2 = xct(i)%libxc_id(2)%code
       if ( (c1 == libxc_ids(1) .and. c2 == libxc_ids(2)) .or.  &
            (c2 == libxc_ids(1) .and. c1 == libxc_ids(2)) ) then
          xc_id = xct(i)
          stat = 0
       endif
    enddo
  end subroutine get_xc_id_from_libxc
         
  subroutine print_xc_id(xc_id)
    type(xc_id_t), intent(in)   :: xc_id
    print "(a)", xc_id_to_string(xc_id)
  end subroutine print_xc_id
    
  function xc_id_to_string(xc_id) result (s)
    type(xc_id_t), intent(in)   :: xc_id
    character(len=72)           :: s
    write(s,"(a,'--',a,1x,a,'--',a,1x,a2)")  &
          trim(xc_id%siestaxc_id%family), &
          trim(xc_id%siestaxc_id%authors), &
          trim(xc_id%libxc_id(1)%name), &
          trim(xc_id%libxc_id(2)%name), &
          trim(xc_id%atom_id)
  end function xc_id_to_string
    
end module m_libxc_sxc_translation

#ifdef __TEST__
program xcid_test

use m_libxc_sxc_translation

type(xc_id_t) :: xc_id
integer       :: stat
character(len=40) :: id, xc_authors, xc_type

do
 write(*,fmt="(a)",advance="no") "Enter string: "
 read(*,"(a)") id
 if (len_trim(id) == 2) then
    call get_xc_id_from_atom_id(trim(id),xc_id,stat)
    if (stat ==0) then
       call print_xc_id(xc_id)
    else
       print *, "UNKNOWN atom_id"
    endif
 else
    write(*,fmt="(a,a)") "String considered as XC type:", trim(id)
    xc_type = id
    write(*,fmt="(a)",advance="no") "Enter XC authors: "
    read(*,"(a)") xc_authors
    call get_xc_id_from_siestaxc(xc_type,xc_authors,xc_id,stat)
    if (stat ==0) then
       call print_xc_id(xc_id)
    else
       print *, "UNKNOWN SiestaXC id"
    endif
 endif
enddo

end program xcid_test
#endif

#ifdef __TEST_LIBXC__
program xcid_test

use m_libxc_sxc_translation

type(xc_id_t) :: xc_id
integer       :: stat, c1, c2

do
 write(*,fmt="(a)",advance="no") "Enter codes: "
 read(*,*) c1, c2
    call get_xc_id_from_libxc((/c1, c2/),xc_id,stat)
    if (stat ==0) then
       call print_xc_id(xc_id)
    else
       print *, "Cannot find functional with those libxc codes"
    endif
enddo

end program xcid_test
#endif
 

   
