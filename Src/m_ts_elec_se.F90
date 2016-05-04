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

    complex(dp) :: E
    integer :: nou, no, nq
    logical :: lnon_Eq

    if ( cE%fake ) return

    call timer('ts_expand',1)

    nou = El%no_used
    no  = TotUsedOrbs(El)
    nq  = product(El%Bloch)
    if ( El%pre_expand > 0 .and. nq > 1 ) then
       nou = no
       nq = 1
    end if

    ! Save energy
    E = cE%e

    lnon_Eq = .false.
    if ( present(non_Eq) ) lnon_Eq = non_Eq

    if ( lnon_Eq ) then
#ifdef TBT_PHONON
       E = dcmplx(real(cE%e,dp)**2,El%Eta)
#else
       E = dcmplx(real(cE%e,dp),El%Eta)
#endif
       call UC_expansion_Sigma_GammaT(E, &
            nou,no,El, nq, &
            El%GA,El%Sigma,El%Gamma,nwork,work)
    else
       if ( El%Bulk ) then
          call UC_expansion_Sigma_Bulk(nou,no,El, nq, &
               El%GA,El%Sigma,nwork,work)
       else
          if ( cE%idx(1) /= 1 ) then ! .not. CONTOUR_EQ
#ifdef TBT_PHONON
             E = dcmplx(real(cE%e,dp)**2,El%Eta)
#else
             E = dcmplx(real(cE%e,dp),El%Eta)
#endif
          end if
          call UC_expansion_Sigma(E,nou,no,El, nq, &
               El%GA,El%Sigma,nwork,work)
       end if
    end if

    call timer('ts_expand',2)

  end subroutine UC_expansion

  subroutine UC_expansion_Sigma_Bulk(no_u,no_s,El, &
       nq,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
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

    if ( nq == 1 ) then
#ifndef TS_NOCHECKS
       if ( no_u /= no_s ) call die('no_E/=no_s')
#endif

       ! When no repetition we save it "as is"
       call zcopy(no_s*no_s,GS(1,1,1),1,Sigma(1,1),1)

    else

#ifndef TS_NOCHECKS
       if ( nwork < no_s ** 2 ) &
            call die('elec_se-Sig-Bulk: worksize too small. Error')
#endif

       if ( El%no_u /= El%no_used ) then
       
          call update_UC_expansion_A(no_u,no_s,El, &
               El%na_used,El%lasto_used,nq,GS,nwork,work(1,1))

          call EYE(no_s,Sigma)
          
          ! We have the matrix to invert in the first no_s**2 values.
          call zgesv(no_s,no_s,work(1,1),no_s,ipvt,Sigma,no_s,ierr)
          if ( ierr /= 0 ) &
               write(*,'(a,i0)') &
               'Inversion of surface Green function failed: ',ierr

       else

          call update_UC_expansion_A(no_u,no_s,El, &
               El%na_used,El%lasto_used,nq,GS,no_s*no_s,Sigma)
          
       end if
                 
    end if

  end subroutine UC_expansion_Sigma_Bulk


  subroutine UC_expansion_Sigma(ZEnergy,no_u,no_s,El, &
       nq,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
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

#ifndef TS_NOCHECKS
    if ( nwork < no_s ** 2 * 2 ) &
         call die('elec_se-Sigma: worksize too small. Error')
#endif

    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         El%na_used,El%lasto_used,nq,El%HA,El%SA,GS,nwork, &
         Sigma(1,1),work(1,1,2))
    
    if ( nq == 1 ) then
#ifndef TS_NOCHECKS
       if ( no_u /= no_s ) call die('no_E/=no_s')
#endif

       ! When no repetition we save it "as is"
       call zcopy(no_s*no_s,GS(1,1,1),1,Sigma(1,1),1)

    else if ( El%no_u /= El%no_used ) then

       ! Invert Sigma only when the electrode size is
       ! reduced, and not pre-expanded

       call zgetrf(no_s, no_s, Sigma, no_s, ipvt, ierr )
       if ( ierr /= 0 ) &
            write(*,'(a,i0)') &
            'Inversion of surface Green (A) function failed: ',ierr
       call zgetri(no_s, Sigma, no_s, ipvt, work(1,1,1), no_s**2, ierr)
       if ( ierr /= 0 ) &
            write(*,'(a,i0)') &
            'Inversion of surface Green (B) function failed: ',ierr
       
    end if

    ! Do:
    ! \Sigma = Z*S - H - \Sigma_bulk
!$OMP parallel do default(shared), private(io,jo), collapse(2)
    do jo = 1 , no_s
       do io = 1 , no_s
          Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
       end do
    end do
!$OMP end parallel do

  end subroutine UC_expansion_Sigma

  subroutine UC_expansion_Sigma_GammaT(ZEnergy,no_u,no_s,El, &
       nq,GS,Sigma,GammaT,nwork,work)
    use intrinsic_missing, only: EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: nq
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
    integer, pointer :: p_G(:)

#ifndef TS_NOCHECKS
    if ( nwork < no_s ** 2 * 2 ) &
         call die('elec_se-GT: worksize too small. Error')
#endif

#ifdef TBTRANS
    call die('elec_se: GT: This routine should never be called in &
         &TBtrans. Will produce erroneous results.')
#endif

    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         El%na_used,El%lasto_used,nq,El%HA,El%SA,GS,nwork, &
         Sigma(1,1), work(1,1,2))

    if ( nq == 1 ) then
#ifndef TS_NOCHECKS
       if ( no_u /= no_s ) call die('no_E/=no_s')
#endif

       ! When no repetition we save it "as is"
       call zcopy(no_s*no_s,GS(1,1,1),1,Sigma(1,1),1)

    else if ( El%no_u /= El%no_used ) then
       
       ! Invert Sigma only when the electrode size is
       ! reduced, and not pre-expanded
       
       call zgetrf(no_s, no_s, Sigma, no_s, ipvt, ierr )
       if ( ierr /= 0 ) &
            write(*,'(a,i0)') &
            'Inversion of surface Green (A) function failed (G): ',ierr
       call zgetri(no_s, Sigma, no_s, ipvt, work(1,1,1), no_s**2, ierr)
       if ( ierr /= 0 ) &
            write(*,'(a,i0)') &
            'Inversion of surface Green (B) function failed (G): ',ierr
       
    end if

    ! Get pivoting table for the scattering matrix
    ! Note that we here pivot directly into the
    ! the same order of the Green function
    ! to not do it "twice"
    p_G => El%inDpvt%r

!$OMP parallel default(shared), private(io,jo)

    if ( El%Bulk ) then

       ! Do:
       ! work = Z*S - H - (Z*S - H - \Sigma_bulk)
!$OMP do collapse(2)
       do jo = 1 , no_s
          do io = 1 , no_s
             work(io,jo,2) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do
!$OMP end do

       ! Do
       ! \Gamma ^ T = \Sigma - \Sigma^\dagger
       if ( associated(p_G) ) then
!$OMP do
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = work(p_G(io),p_G(jo),2) &
                  - dconjg(work(p_G(jo),p_G(io),2))
             GammaT(io,jo) = work(p_G(jo),p_G(io),2) &
                  - dconjg(work(p_G(io),p_G(jo),2))
          end do
          io = p_G(jo)
          GammaT(jo,jo) = work(io,io,2)-dconjg(work(io,io,2))
       end do
!$OMP end do nowait
       else
!$OMP do
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = work(io,jo,2) &
                  - dconjg(work(jo,io,2))
             GammaT(io,jo) = work(jo,io,2) &
                  - dconjg(work(io,jo,2))
          end do
          GammaT(jo,jo) = work(jo,jo,2)-dconjg(work(jo,jo,2))
       end do
!$OMP end do nowait
       end if
    else
       
       ! Do:
       ! \Sigma = Z*S - H - (Z*S - H - \Sigma_bulk)
!$OMP do collapse(2)
       do jo = 1 , no_s
          do io = 1 , no_s
             Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do
!$OMP end do 

       ! Do
       ! \Gamma ^ T = \Sigma - \Sigma^\dagger
       if ( associated(p_G) ) then
!$OMP do 
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = Sigma(p_G(io),p_G(jo)) &
                  - dconjg(Sigma(p_G(jo),p_G(io)))
             GammaT(io,jo) = Sigma(p_G(jo),p_G(io)) &
                  - dconjg(Sigma(p_G(io),p_G(jo)))
          end do
          io = p_G(jo)
          GammaT(jo,jo) = Sigma(io,io)-dconjg(Sigma(io,io))
       end do
!$OMP end do nowait
       else
!$OMP do 
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = Sigma(io,jo) &
                  - dconjg(Sigma(jo,io))
             GammaT(io,jo) = Sigma(jo,io) &
                  - dconjg(Sigma(io,jo))
          end do
          GammaT(jo,jo) = Sigma(jo,jo)-dconjg(Sigma(jo,jo))
       end do
!$OMP end do nowait
       end if
    end if

!$OMP end parallel 

  end subroutine UC_expansion_Sigma_GammaT

  subroutine update_UC_expansion(ZEnergy,no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,nwork,work1,work2)
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
    complex(dp), intent(inout) :: work1(no_s,no_s), work2(no_s,no_s)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq
    integer :: iow,iau,ia1,ia2,ia3,iuo
    integer :: jow,jau,ja1,ja2,ja3,juo
    complex(dp) :: p(3), pZ, qPi
    real(dp) :: rPi(3), wq

    if ( nq == 1 ) then
#ifndef TS_NOCHECKS
       if ( no_u /= no_s ) call die('no_E/=no_s')
#endif

       ! In case the pre-expansion is not done on H, S
       if ( El%pre_expand == 1 .and. product(El%Bloch) > 1 ) then

          iq = product(El%Bloch)

          ! Note that this is because the interface for H and S
          iow = El%no_used
!$OMP parallel do default(shared), private(iuo,juo), collapse(2)
          do juo = 1 , iow
             do iuo = 1 , no_s
                work1(iuo,juo) = ZEnergy * S(iuo,juo,1) - H(iuo,juo,1)
             end do
          end do
!$OMP end parallel do

          call update_UC_expansion_A(iow,no_s,El,na_u,lasto,&
               iq,work1(1,1),nwork,work2(1,1))
          
       else

          ! We do not need to copy over GS, as it is
          ! used correctly
       
!$OMP parallel do default(shared), private(iuo,juo), collapse(2)
          do juo = 1 , no_s
             do iuo = 1 , no_s
               !work1(iuo,juo,1) = GS(iuo,juo,1)
                work2(iuo,juo) = ZEnergy * S(iuo,juo,1) - H(iuo,juo,1)
             end do
          end do
!$OMP end parallel do

       end if

    else

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green function
!$OMP parallel default(shared), private(iq,wq,rPi,qPi,p,pZ), &
!$OMP&private(jow,jau,ja1,ja2,ja3,juo)

       ! Save some multiplications
       wq = log(1._dp / real(nq,dp))

       rPi = 2._dp * Pi * (q_exp(El,1) + El%bkpt_cur)
       qPi = cdexp(dcmplx(0._dp,rPi(1)))

! Create the single execution for creating tasks
!$OMP single
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%Bloch(3)
        do ia2 = 1 , El%Bloch(2)
        do ia1 = 1 , El%Bloch(1)
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
!$OMP task firstprivate(iow,iuo,ia1,ia2,ia3)
           p(3) = cdexp(dcmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3)))
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%Bloch(3)
            p(2) = p(3)*cdexp(dcmplx(0._dp,ja3*rPi(3)))
            do ja2 = 1 , El%Bloch(2)
            p(1) = p(2)*cdexp(dcmplx(0._dp,ja2*rPi(2)))
            do ja1 = 1 , El%Bloch(1)
              ! This takes one additional phase per iteration
              p(1) = p(1)*qPi
              pZ = p(1) * ZEnergy
              do juo = 1 + lasto(jau-1) , lasto(jau)
               jow = jow + 1
                
               work1(jow,iow) = p(1) * GS(juo,iuo,1)
                
               work2(jow,iow) = pZ * S(juo,iuo,1) - p(1) * H(juo,iuo,1)
                
              end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
!$OMP end task
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau

!$OMP end single ! keep implicit barrier (rPi,qPi)

       do iq = 2 , nq

        rPi = 2._dp * Pi * (q_exp(El,iq) + El%bkpt_cur)
        qPi = cdexp(dcmplx(0._dp,rPi(1)))

!$OMP single

        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%Bloch(3)
         do ia2 = 1 , El%Bloch(2)
         do ia1 = 1 , El%Bloch(1)
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
!$OMP task firstprivate(iow,iuo,ia1,ia2,ia3)
            p(3) = cdexp(dcmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3)))
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%Bloch(3)
             p(2) = p(3)*cdexp(dcmplx(0._dp,ja3*rPi(3)))
             do ja2 = 1 , El%Bloch(2)
             p(1) = p(2)*cdexp(dcmplx(0._dp,ja2*rPi(2)))
             do ja1 = 1 , El%Bloch(1)
               p(1) = p(1)*qPi
               pZ = p(1) * ZEnergy
               do juo = 1 + lasto(jau-1) , lasto(jau)
                  jow = jow + 1
                  
                  work1(jow,iow) = work1(jow,iow) + p(1) * GS(juo,iuo,iq)
   
                  work2(jow,iow) = work2(jow,iow) + &
                       pZ * S(juo,iuo,iq) - p(1) * H(juo,iuo,iq)
   
               end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
!$OMP end task
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
        end do !iau

!$OMP end single ! keep implicit barrier (rPi,qPi)

       end do !q-points

!$OMP end parallel

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
    integer :: iau,iow,ia1,ia2,ia3,iuo
    integer :: jau,jow,ja1,ja2,ja3,juo
    complex(dp) :: p(3), qPi
    real(dp) :: rPi(3), wq

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')
       
       call die('This should probably not be called')

    else

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green function
!$OMP parallel default(shared), private(iq,wq,rPi,qPi,p), &
!$OMP&private(jow,jau,ja1,ja2,ja3,juo)

       ! Save some multiplications
       wq = log(1._dp / real(nq,dp))

       rPi = 2._dp * Pi * (q_exp(El,1) + El%bkpt_cur)
       qPi = cdexp(dcmplx(0._dp,rPi(1)))

! Create the single execution for creating tasks
!$OMP single
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%Bloch(3)
        do ia2 = 1 , El%Bloch(2)
        do ia1 = 1 , El%Bloch(1)
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
!$OMP task firstprivate(iuo,iow,ia1,ia2,ia3)
           p(3) = cdexp(dcmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3)))
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%Bloch(3)
            p(2) = p(3)*cdexp(dcmplx(0._dp,ja3*rPi(3)))
            do ja2 = 1 , El%Bloch(2)
            p(1) = p(2)*cdexp(dcmplx(0._dp,ja2*rPi(2)))
            do ja1 = 1 , El%Bloch(1)
             p(1) = p(1)*qPi
             do juo = 1 + lasto(jau-1) , lasto(jau)
               jow = jow + 1
               
               work(jow,iow) = p(1) * A(juo,iuo,1)
               
             end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
!$OMP end task
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau

!$OMP end single ! keep implicit barrier (rPi,qPi)

       do iq = 2 , nq

        rPi = 2._dp * Pi * (q_exp(El,iq) + El%bkpt_cur)
        qPi = cdexp(dcmplx(0._dp,rPi(1)))

!$OMP single

        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%Bloch(3)
         do ia2 = 1 , El%Bloch(2)
         do ia1 = 1 , El%Bloch(1)
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
!$OMP task firstprivate(iuo,iow,ia1,ia2,ia3)
            p(3) = cdexp(dcmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3)))
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%Bloch(3)
             p(2) = p(3)*cdexp(dcmplx(0._dp,ja3*rPi(3)))
             do ja2 = 1 , El%Bloch(2)
             p(1) = p(2)*cdexp(dcmplx(0._dp,ja2*rPi(2)))
             do ja1 = 1 , El%Bloch(1)
              p(1) = p(1)*qPi
              do juo = 1 + lasto(jau-1) , lasto(jau)
                jow = jow + 1
                
                work(jow,iow) = work(jow,iow) + p(1) * A(juo,iuo,iq)
                 
              end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
!$OMP end task
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
        end do !iau

!$OMP end single ! keep implicit barrier (rPi,qPi)

       end do !q-points

!$OMP end parallel

    end if

  end subroutine update_UC_expansion_A

end module m_ts_elec_se
