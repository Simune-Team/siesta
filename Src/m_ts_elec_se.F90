! This module will control what is used of the electrodes in the transiesta SCF
!
! Hence we here collect the routines for reading and expanding the self-energies
! in the GF-files.

module m_ts_elec_se

  use precision, only : dp

  use m_ts_electype
  use m_ts_cctype

  implicit none

  private

  public :: UC_expansion
  !public :: UC_expansion_Sigma_Bulk
  !public :: UC_expansion_Sigma
  !public :: UC_expansion_Sigma_GammaT
  public :: update_UC_expansion_A

contains

  subroutine UC_expansion(cE,El,nwork,work, non_Eq)
! ********************
! * INPUT variables  *
! ********************
    type(ts_c_idx), intent(in) :: cE
    type(Elec), intent(in out) :: El

! ********************
! * WORK variables   *
! ********************
    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)

    logical,  intent(in), optional :: non_Eq

    integer :: nou, no, nq
    logical :: lnon_Eq

    if ( cE%fake ) return

    call timer('ts_expand',1)

    nou = El%no_used
    no  = TotUsedOrbs(El)
    nq  = Rep(El)
    if ( El%pre_expand .and. nq > 1 ) then
       nou = no
       nq = 1
    end if

    lnon_Eq = .false.
    if ( present(non_Eq) ) lnon_Eq = non_Eq

    if ( lnon_Eq ) then
       call UC_expansion_Sigma_GammaT(cE%e, &
            nou,no,El, nq, &
            El%HA,El%SA,El%GA,El%Sigma,El%Gamma,nwork,work)
    else
       if ( El%Bulk ) then
          call UC_expansion_Sigma_Bulk(nou,no,El, nq, &
               El%HA,El%SA,El%GA,El%Sigma,nwork,work)
       else
          call UC_expansion_Sigma(cE%e,nou,no,El, nq, &
               El%HA,El%SA,El%GA,El%Sigma,nwork,work)
       end if
    end if

    call timer('ts_expand',2)

  end subroutine UC_expansion

  subroutine UC_expansion_Sigma_Bulk(no_u,no_s,El, &
       nq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(inout) :: Sigma(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s)
! ********************
! * LOCAL variables  *
! ********************
    integer :: ierr
    integer :: ipvt(no_s)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2 ) call die('Size of work-array is &
         &too small')

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)

    else

       call update_UC_expansion_A(no_u,no_s,El, &
            El%na_used,El%lasto_used,nq,GS,nwork,work(1,1))
       
       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1),no_s,ipvt,Sigma,no_s,ierr)
       if ( ierr /= 0 ) &
            write(*,*) 'Inversion of surface Greens function failed'
       
    end if

  end subroutine UC_expansion_Sigma_Bulk


  subroutine UC_expansion_Sigma(ZEnergy,no_u,no_s,El, &
       nq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(inout) :: Sigma(no_s,no_s)

    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    integer :: ierr
    integer :: io, jo
    integer :: ipvt(no_s)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')
       
    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         El%na_used,El%lasto_used,nq,H,S,GS,nwork,work(1,1,1))

    if ( nq == 1 ) then

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)

    else
       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)

       if ( ierr /= 0 ) &
            write(*,*) 'Inversion of surface Greens function failed'

    end if

    ! Do:
    ! \Sigma = Z*S - H - \Sigma_bulk
    do jo = 1 , no_s
       do io = 1 , no_s
          Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
       end do
    end do

  end subroutine UC_expansion_Sigma

  subroutine UC_expansion_Sigma_GammaT(ZEnergy,no_u,no_s,El, &
       nq,H,S,GS,Sigma,GammaT,nwork,work)
    use intrinsic_missing, only: EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(inout) :: Sigma(no_s,no_s)
    complex(dp), intent(inout) :: GammaT(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    complex(dp), parameter :: zi = dcmplx(0._dp,1._dp)
    integer :: ierr
    integer :: io,jo
    integer :: ipvt(no_s)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         El%na_used,El%lasto_used,nq,H,S,GS,nwork,work(1,1,1))

    if ( nq == 1 ) then

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)
       
    else

       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'Inversion of surface Greens function failed'
       end if

    end if

    if ( El%Bulk ) then

       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       do jo = 1 , no_s
          do io = 1 , no_s
             work(io,jo,1) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = i ( \Sigma - \Sigma^\dagger)
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = zi * ( &
                  work(io,jo,1)-dconjg(work(jo,io,1)) )
             GammaT(io,jo) = zi * ( &
                  work(jo,io,1)-dconjg(work(io,jo,1)) )
          end do
          GammaT(jo,jo) = zi * ( &
               work(jo,jo,1)-dconjg(work(jo,jo,1)) )
       end do

    else
       
       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       do jo = 1 , no_s
          do io = 1 , no_s
             Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = i ( \Sigma - \Sigma^\dagger)
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = zi * ( &
                  Sigma(io,jo)-dconjg(Sigma(jo,io)) )
             GammaT(io,jo) = zi * ( &
                  Sigma(jo,io)-dconjg(Sigma(io,jo)) )
          end do
          GammaT(jo,jo) = zi * ( &
               Sigma(jo,jo)-dconjg(Sigma(jo,jo)) )
       end do

    end if

  end subroutine UC_expansion_Sigma_GammaT

  subroutine update_UC_expansion(ZEnergy,no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,nwork,work)
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u), nq
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S, GS
! ********************
! * OUTPUT variables *
! ********************
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq
    integer :: iow,iau,ia3,ia2,ia1,iuo
    integer :: jow,jau,ja3,ja2,ja1,juo
    complex(dp) :: ph
    real(dp) :: qPi(3,nq), wq

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! We do not need to copy over GS, as it is
       ! used correctly
       !work(:,:,1) = GS(:,:,1)
       work(:,:,2) = ZEnergy * S(:,:,1) - H(:,:,1)

    else

       do iq = 1 , nq 
          qPi(1:3,iq) = 2._dp * Pi * (q_exp(El,iq) + El%bkpt_cur)
       end do
       wq = 1._dp / real(nq,dp)

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%RepA3
        do ia2 = 1 , El%RepA2
        do ia1 = 1 , El%RepA1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%RepA3
            do ja2 = 1 , El%RepA2
            do ja1 = 1 , El%RepA1
              ph = wq * cdexp(dcmplx(0._dp, &
                   (ja1-ia1)*qPi(1,iq) + (ja2-ia2)*qPi(2,iq) + (ja3-ia3)*qPi(3,iq) ) )
              do juo = 1 + lasto(jau-1) , lasto(jau)
                jow = jow + 1
                
                work(jow,iow,1) = ph * GS(juo,iuo,iq)
                
                work(jow,iow,2) = ph * (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
                
              end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%RepA3
         do ia2 = 1 , El%RepA2
         do ia1 = 1 , El%RepA1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%RepA3
             do ja2 = 1 , El%RepA2
             do ja1 = 1 , El%RepA1
               ph = wq * cdexp(dcmplx(0._dp, &
                    (ja1-ia1)*qPi(1,iq) + (ja2-ia2)*qPi(2,iq) + (ja3-ia3)*qPi(3,iq) ) )
               do juo = 1 + lasto(jau-1) , lasto(jau)
                  jow = jow + 1
                  
                  work(jow,iow,1) = work(jow,iow,1) + ph * GS(juo,iuo,iq)
   
                  work(jow,iow,2) = work(jow,iow,2) + ph * &
                       (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
   
               end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
        end do !iau
       end do !q-points

    end if

  end subroutine update_UC_expansion

  subroutine update_UC_expansion_A(no_u,no_s,El, &
       na_u,lasto,nq,A,nwork,work)
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u), nq
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: A
! ********************
! * OUTPUT variables *
! ********************
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq
    integer :: iow,iau,ia3,ia2,ia1,iuo
    integer :: jow,jau,ja3,ja2,ja1,juo
    complex(dp) :: ph
    real(dp) :: qPi(3,nq), wq

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')
       
       call die('This should probably not be called')

    else

       do iq = 1 , nq 
          qPi(1:3,iq) = 2._dp * Pi * (q_exp(El,iq) + El%bkpt_cur)
       end do
       wq = 1._dp / real(nq,dp)

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%RepA3
        do ia2 = 1 , El%RepA2
        do ia1 = 1 , El%RepA1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%RepA3
            do ja2 = 1 , El%RepA2
            do ja1 = 1 , El%RepA1
               ph = wq * cdexp(dcmplx(0._dp, &
                    (ja1-ia1)*qPi(1,iq) + (ja2-ia2)*qPi(2,iq) + (ja3-ia3)*qPi(3,iq) ) )
              do juo = 1 + lasto(jau-1) , lasto(jau)
                jow = jow + 1
                
                work(jow,iow) = ph * A(juo,iuo,iq)
                
              end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%RepA3
         do ia2 = 1 , El%RepA2
         do ia1 = 1 , El%RepA1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%RepA3
             do ja2 = 1 , El%RepA2
             do ja1 = 1 , El%RepA1
               ph = wq * cdexp(dcmplx(0._dp, &
                    (ja1-ia1)*qPi(1,iq) + (ja2-ia2)*qPi(2,iq) + (ja3-ia3)*qPi(3,iq) ) )
               do juo = 1 + lasto(jau-1) , lasto(jau)
                  jow = jow + 1
                  
                  work(jow,iow) = work(jow,iow) + ph * A(juo,iuo,iq)
   
               end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
        end do !iau
       end do !q-points

    end if

  end subroutine update_UC_expansion_A

end module m_ts_elec_se
