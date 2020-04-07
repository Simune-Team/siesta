! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!< Bloch expansion module implementing a way to expand a Hermitian matrix obeying Bloch's theorem.
!!
!! Bloch's theorem states that only a constant phase-factor is accummulated upon a translation
!! vector `T`.
!!
!! This module implements a type `bloch_unfold_t` implementing the required algorithms
!! for unfolding a matrix.
module bloch_unfold_m

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(p=15)

  type bloch_unfold_t
    
    !< Number of times a matrix is unfolded
    integer, public :: B(3) = 1

    !< Local variable for quick query of product(this%B)
    integer, private :: prod_B = 1

  contains

    procedure, pass :: initialize => bloch_unfold_init
    procedure, pass :: size => bloch_unfold_size
    procedure, pass :: matrix_unfold => bloch_unfold_matrix_z
    procedure, pass :: matrix_unfold_HS_G => bloch_unfold_matrix_HS_G_z
    procedure, pass :: matrix_unfold_HS => bloch_unfold_matrix_HS_z
    procedure, pass :: get_k => bloch_unfold_get_k

  end type bloch_unfold_t
  public :: bloch_unfold_t

contains

  subroutine bloch_unfold_init(this, B)
    class(bloch_unfold_t), intent(inout) :: this
    integer, intent(in) :: B(3)

    this%B(:) = B(:)

    if ( any(this%B < 1) ) then
      call die('bloch_unfold requires unfolding to be >= 1!')
    end if

    this%prod_B = product(this%B)

  end subroutine bloch_unfold_init

  pure function bloch_unfold_size(this) result(size)
    class(bloch_unfold_t), intent(in) :: this
    integer :: size
    
    size = this%prod_B
    
  end function bloch_unfold_size

  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_matrix_z(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return
      
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_matrix_z_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_z_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_matrix_z_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_z_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_matrix_z_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_matrix_z_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_matrix_z

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_z_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: ph, phc, phA_step

    ! Initialize un-folded matrix
    uM(:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,phc)

    w = 1._dp / real(NA, dp)
    
    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)
      
      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
!$OMP do
      do j = 1, N
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        phc = conjg(ph)
!$OMP do
        do j = 1, N
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * phc
          end do
        end do
!$OMP end do nowait
        ph = ph * phA_step
      end do

    end do

!$OMP barrier

    ! At this point the following has been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uM(:,TMA,:,iMA) = uM(:,TMA-1,:,iMA-1)
      end do
!$OMP end do
      ! we need to wait because we copy the previous columns
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_z_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_z_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, phc

    uM(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,phc)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)
    
    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
!$OMP do
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA
          ph = phA
          phc = conjg(phA)
!$OMP do
          do j = 1, N
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait
          phA = phA * phA_step
        end do
        
        
        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          phc = conjg(phB)
!$OMP do
          do j = 1, N
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA
            ph = phB * phA
            phc = phB * conjg(phA)
!$OMP do
            do j = 1, N
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait

            ph = conjg(phB) * phA
            phc = conjg(phB) * conjg(phA)
!$OMP do
            do j = 1, N
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait
            
            phA = phA * phA_step
          end do
          
          phB = phB * phB_step
        end do
      end do
      
    end do

!$OMP barrier

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uM(:,TMA,1,:,iMA,1) = uM(:,TMA-1,1,:,iMA-1,1)
      end do
!$OMP end do
      ! we need a wait since it copies from the previous iMA
    end do
    do iMB = 2, NB
      do iMA = 2, NA
!$OMP do schedule(static)
        do TMA = 2, NA
          uM(:,TMA,iMB,:,iMA,1) = uM(:,TMA-1,iMB,:,iMA-1,1)
          uM(:,TMA,1,:,iMA,iMB) = uM(:,TMA-1,1,:,iMA-1,iMB)
        end do
!$OMP end do
        ! we need a wait since it copies from the previous iMA
      end do
    end do

    ! Now copy all duplicated data
    do iMB = 2, NB
!$OMP do schedule(static)
      do TMB = 2, NB
        do TMA = 1, NA
          uM(:,TMA,TMB,:,TMA,iMB) = uM(:,TMA,TMB-1,:,TMA,iMB-1)
        end do
      end do
!$OMP end do
      ! we need a wait since it copies from the previous iMB
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_z_2
  
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_matrix_HS_G_z(this, bk, N, H, S, G, Z, uSZmH, uG)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH, uG

    if ( this%prod_B == 1 ) then
      
      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)
      call zcopy(N*N, G(1,1,1), 1, uG(1,1,1,1), 1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_matrix_HS_G_z_1(N, H, S, G, Z, bk(3), this%B(3), uSZmH, uG)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_HS_G_z_1(N, H, S, G, Z, bk(2), this%B(2), uSZmH, uG)
      else
        call bloch_unfold_matrix_HS_G_z_2(N, H, S, G, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_HS_G_z_1(N, H, S, G, Z, bk(1), this%B(1), uSZmH, uG)
      else
        call bloch_unfold_matrix_HS_G_z_2(N, H, S, G, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_matrix_HS_G_z_2(N, H, S, G, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH, uG)
    else
      call die('currently not implemented')
    end if
    
  end subroutine bloch_unfold_matrix_HS_G_z

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_HS_G_z_1(N, H, S, G, Z, kA, NA, uSZmH, uG)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), dimension(N,NA,N,NA), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: phA_step
    complex(dp) :: ph, phc, phZ, phcZ

    ! Initialize un-folded matrix
    uSZmH(:,:,:,:) = cmplx(0._dp, 0._dp, dp)
    uG(:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,phc,phZ,phcZ) firstprivate(Z)

    w = 1._dp / real(NA, dp)
    
    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)
      
      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
      phZ = Z * w
!$OMP do
      do j = 1, N
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * w
        end do
        do i = 1, N
          uG(i,1,j,1) = uG(i,1,j,1) + G(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        phc = conjg(ph)
        phZ = ph * Z
        phcZ = phc * Z
!$OMP do
        do j = 1, N
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + S(i,j,TMA) * phcZ - H(i,j,TMA) * phc
          end do
          do i = 1, N
            uG(i,iMA,j,1) = uG(i,iMA,j,1) + G(i,j,TMA) * ph
          end do
          do i = 1, N
            uG(i,1,j,iMA) = uG(i,1,j,iMA) + G(i,j,TMA) * phc
          end do
        end do
!$OMP end do nowait
        ph = ph * phA_step
      end do

    end do

!$OMP barrier

    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uSZmH(:,TMA,:,iMA) = uSZmH(:,TMA-1,:,iMA-1)
        uG(:,TMA,:,iMA) = uG(:,TMA-1,:,iMA-1)
      end do
!$OMP end do
      ! we need to wait because we copy the previous columns
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_HS_G_z_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_HS_G_z_2(N, H, S, G, Z, kA, NA, kB, NB, uSZmH, uG)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), dimension(N,NA,NB,N,NA,NB), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, phc, phZ, phcZ

    uSZmH(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)
    uG(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,phc,phZ,phcZ) firstprivate(Z)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
        phZ = Z * w
!$OMP do
        do j = 1, N
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
          do i = 1, N
            uG(i,1,1,j,1,1) = uG(i,1,1,j,1,1) + G(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA          
          ph = phA
          phc = conjg(phA)
          phZ = ph * Z
          phcZ = phc * Z
!$OMP do
          do j = 1, N
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
            end do

            do i = 1, N
              uG(i,iMA,1,j,1,1) = uG(i,iMA,1,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,iMA,1) = uG(i,1,1,j,iMA,1) + G(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait
          phA = phA * phA_step
        end do
        
        
        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          phc = conjg(phB)
          phZ = ph * Z
          phcZ = phc * Z
!$OMP do
          do j = 1, N
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
            end do
            
            do i = 1, N
              uG(i,1,iMB,j,1,1) = uG(i,1,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,1,iMB) = uG(i,1,1,j,1,iMB) + G(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA
            ph = phB * phA
            phc = phB * conjg(phA)
            phZ = ph * Z
            phcZ = phc * Z
!$OMP do
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
              end do
              
              do i = 1, N
                uG(i,iMA,iMB,j,1,1) = uG(i,iMA,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,iMB,j,iMA,1) = uG(i,1,iMB,j,iMA,1) + G(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait

            ph = conjg(phB) * phA
            phc = conjg(phB) * conjg(phA)
            phZ = ph * Z
            phcZ = phc * Z
!$OMP do
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
              end do
              
              do i = 1, N
                uG(i,iMA,1,j,1,iMB) = uG(i,iMA,1,j,1,iMB) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,1,j,iMA,iMB) = uG(i,1,1,j,iMA,iMB) + G(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait
            
            phA = phA * phA_step
          end do
          
          phB = phB * phB_step
        end do
      end do
      
    end do

!$OMP barrier

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uSZmH(:,TMA,1,:,iMA,1) = uSZmH(:,TMA-1,1,:,iMA-1,1)
        uG(:,TMA,1,:,iMA,1) = uG(:,TMA-1,1,:,iMA-1,1)
      end do
!$OMP end do
        ! we need a wait since it copies from the previous iMA
    end do
    do iMB = 2, NB
      do iMA = 2, NA
!$OMP do schedule(static)
        do TMA = 2, NA
          uSZmH(:,TMA,1,:,iMA,iMB) = uSZmH(:,TMA-1,1,:,iMA-1,iMB)
          uSZmH(:,TMA,iMB,:,iMA,1) = uSZmH(:,TMA-1,iMB,:,iMA-1,1)
          uG(:,TMA,1,:,iMA,iMB) = uG(:,TMA-1,1,:,iMA-1,iMB)
          uG(:,TMA,iMB,:,iMA,1) = uG(:,TMA-1,iMB,:,iMA-1,1)
        end do
!$OMP end do
        ! we need a wait since it copies from the previous iMA
      end do
    end do

    ! Now copy all duplicated data
    do iMB = 2, NB
!$OMP do schedule(static)
      do TMB = 2, NB
        do TMA = 1, NA
          uSZmH(:,TMA,TMB,:,TMA,iMB) = uSZmH(:,TMA,TMB-1,:,TMA,iMB-1)
          uG(:,TMA,TMB,:,TMA,iMB) = uG(:,TMA,TMB-1,:,TMA,iMB-1)
        end do
      end do
!$OMP end do
      ! we need a wait since it copies from the previous iMB
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_HS_G_z_2


  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_matrix_HS_z(this, bk, N, H, S, Z, uSZmH)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH

    if ( this%prod_B == 1 ) then
      
      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_matrix_HS_z_1(N, H, S, Z, bk(3), this%B(3), uSZmH)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_HS_z_1(N, H, S, Z, bk(2), this%B(2), uSZmH)
      else
        call bloch_unfold_matrix_HS_z_2(N, H, S, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_matrix_HS_z_1(N, H, S, Z, bk(1), this%B(1), uSZmH)
      else
        call bloch_unfold_matrix_HS_z_2(N, H, S, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_matrix_HS_z_2(N, H, S, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH)
    else
      call die('currently not implemented')
    end if
    
  end subroutine bloch_unfold_matrix_HS_z

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_HS_z_1(N, H, S, Z, kA, NA, uSZmH)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uSZmH(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: phA_step
    complex(dp) :: ph, phc, phZ, phcZ

    ! Initialize un-folded matrix
    uSZmH(:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,phc,phZ,phcZ) firstprivate(Z)

    w = 1._dp / real(NA, dp)
    
    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)
      
      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
      phZ = Z * w
!$OMP do
      do j = 1, N
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        phc = conjg(ph)
        phZ = Z * ph
        phcZ = Z * phc
!$OMP do
        do j = 1, N
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + S(i,j,TMA) * phcZ - H(i,j,TMA) * phc
          end do
        end do
!$OMP end do nowait
        ph = ph * phA_step
      end do

    end do

!$OMP barrier

    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uSZmH(:,TMA,:,iMA) = uSZmH(:,TMA-1,:,iMA-1)
      end do
!$OMP end do
      ! we need to wait because we copy the previous columns
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_HS_z_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_matrix_HS_z_2(N, H, S, Z, kA, NA, kB, NB, uSZmH)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uSZmH(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, phc, phZ, phcZ

    uSZmH(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,phc,phZ,phcZ) firstprivate(Z)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
        phZ = Z * w
!$OMP do
        do j = 1, N
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA          
          ph = phA
          phc = conjg(phA)
          phZ = Z * ph
          phcZ = Z * phc
!$OMP do
          do j = 1, N
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait
          phA = phA * phA_step
        end do
        
        
        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          phc = conjg(phB)
          phZ = Z * ph
          phcZ = Z * phc
!$OMP do
          do j = 1, N
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
            end do
          end do
!$OMP end do nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA
            ph = phB * phA
            phc = phB * conjg(phA)
            phZ = Z * ph
            phcZ = Z * phc
!$OMP do
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait

            ph = conjg(phB) * phA
            phc = conjg(phB) * conjg(phA)
            phZ = Z * ph
            phcZ = Z * phc
!$OMP do
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * phcZ - H(i,j,TMA,TMB) * phc
              end do
            end do
!$OMP end do nowait
            
            phA = phA * phA_step
          end do
          
          phB = phB * phB_step
        end do
      end do
      
    end do

!$OMP barrier

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
!$OMP do schedule(static)
      do TMA = 2, NA
        uSZmH(:,TMA,1,:,iMA,1) = uSZmH(:,TMA-1,1,:,iMA-1,1)
      end do
!$OMP end do
      ! we need a wait since it copies from the previous iMA
    end do
    do iMB = 2, NB
      do iMA = 2, NA
!$OMP do schedule(static)
        do TMA = 2, NA
          uSZmH(:,TMA,iMB,:,iMA,1) = uSZmH(:,TMA-1,iMB,:,iMA-1,1)
          uSZmH(:,TMA,1,:,iMA,iMB) = uSZmH(:,TMA-1,1,:,iMA-1,iMB)
        end do
!$OMP end do
        ! we need a wait since it copies from the previous iMA
      end do
    end do

    ! Now copy all duplicated data
    do iMB = 2, NB
!$OMP do schedule(static)
      do TMB = 2, NB
        do TMA = 1, NA
          uSZmH(:,TMA,TMB,:,TMA,iMB) = uSZmH(:,TMA,TMB-1,:,TMA,iMB-1)
        end do
      end do
!$OMP end do
      ! we need a wait since it copies from the previous iMB
    end do

!$OMP end parallel

  end subroutine bloch_unfold_matrix_HS_z_2

  pure function bloch_unfold_get_k(this, iA, iB, iC) result(bk)
    class(bloch_unfold_t), intent(in) :: this
    integer, intent(in) :: iA, iB, iC
    real(dp) :: bk(3)

    ! TODO, the current implementation assumes k-symmetry
    ! of the electrode electronic structure.
    ! Using Bloch expansion with non-symmetry will, likely, produce
    ! wrong results.
    ! Luckily this is not a problem currently.
    ! Perhaps one should consider this in tbtrans

    bk(1) = real(iA - 1, dp) / real(this%B(1), dp)
    bk(2) = real(iB - 1, dp) / real(this%B(2), dp)
    bk(3) = real(iC - 1, dp) / real(this%B(3), dp)
    
  end function bloch_unfold_get_k

  pure function unfold_k2pi(bk, Nk, ik) result(k)
    use units, only: Pi
    real(dp), intent(in) :: bk
    integer, intent(in) :: Nk, ik
    real(dp) :: k
    ! bk should already be in *local* coordinates
    k = 2._dp * Pi * (bk + real(ik - 1, dp) / real(Nk, dp))
  end function unfold_k2pi

end module bloch_unfold_m
