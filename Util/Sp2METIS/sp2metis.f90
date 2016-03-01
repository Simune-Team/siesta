! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
!  Nick Papior, 2015

! This program will print a ndmetis
! prepared graph file
program sp2metis

  use class_Sparsity
  use geom_helper, only : UCORB
  use m_io_s
  use m_os, only : file_exist
  use create_sparsity_SC

  character(len=250) :: fname
  character(len=20) :: fmt
  
  ! The sparsity
  type(Sparsity) :: sp_full, sp_uc
  integer :: no_u, io, o, j, o2, jo
  integer, pointer :: ncol(:), l_ptr(:), l_col(:)
  integer, pointer :: cur_col(:), col(:)
  integer, allocatable :: cw(:)
  integer :: N_arg
  logical :: has_weight

  has_weight = .false.

  N_arg = command_argument_count()
  ! Get the file-name from the command-line
  if ( N_arg < 1 ) then
     stop "Call this program with a file-name *.DM, *.TSDE"
  end if
  if ( N_arg == 1 ) then
     call get_command_argument(1, value = fname)
  else 
     call get_command_argument(1, value = fname)
     if ( fname == '-w' ) has_weight = .true.
     call get_command_argument(2, value = fname)
  end if

  ! open file
  if ( .not. file_exist(fname) ) then
     print *,trim(fname)
     stop "File does not exist. Please provide file that exists..."
  end if

  j = 10
  open(j, file = trim(fname) , status = 'old', form = 'unformatted' )

  ! read the first line
  read(j) no_u, io ! nspin

  call io_read_Sp(j, no_u, sp_full, 'sp')

  close(j)

  ! Convert to UC
  call crtSparsity_SC(sp_full,sp_uc, UC = .true. )
  call delete(sp_full)

  ! Convert sparsity pattern to graph
  call attach(sp_uc, n_col = ncol, list_ptr = l_ptr, &
       list_col = l_col, nnzs = n_nzs )

  ! Start writing it out
  if ( has_weight ) then
     write(*,'(2(tr1,i0),tr1,a)') no_u , (n_nzs - no_u)/2,'001'
  else
     write(*,'(2(tr1,i0))') no_u , (n_nzs - no_u)/2
  end if

  ! Write out each graph point
  write(fmt,'(i0,a)') 2*no_u,'(tr1,i0)'
  do io = 1 , no_u

     cur_col => l_col(l_ptr(io)+1:l_ptr(io)+ncol(io))

     ! allocate graph and weight
     allocate(cw((ncol(io)-1) * 2))
     
     ! Create the weights graph
     o2 = 0
     do o = 1 , ncol(io)

        jo = cur_col(o)
        if ( jo == io ) cycle

        ! step position
        o2 = o2 + 2

        ! copy over position
        cw(o2-1) = jo

        ! create weights
        col => l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo))
        cw(o2) = 0

        do j = 1 , ncol(jo)
           if ( any(col(j) == cur_col) ) cw(o2) = cw(o2) + 1
        end do

     end do
     
     if ( has_weight ) then
        write(*,'('//trim(fmt)//')') cw
     else
        write(*,'('//trim(fmt)//')') pack(cur_col, cur_col /= io)
     end if

     deallocate(cw)
  end do

  call delete(sp_uc)

end program sp2metis
