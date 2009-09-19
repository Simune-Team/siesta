C...............................................................
C
      subroutine write_axsf1(ii2,io1,nstep,mdfirst,mdlast,mdstep,
     .                       nat,nz,mdmod,varcel,cc_ang,coord,veloc,
     .                       obox,rbox,rinv)
C
C     writes output info into the AXSF file for the case
C     with output box.
C     The AXSF is as for periodic structure, with fixed cell
C     originating at obox and spanned by three vectors rbox.
C     Read again nstep MD steps from ii2,
C     write those from mdfirst till mdlast, step mdstep.
C
      implicit none
      integer ii2,io1,is1,mdfirst,mdlast,mdstep,nat,nz(nat),mdmod,
     .        ii,jj,istep,nstep,idum,nbox,ibox,iat
      double precision cc_ang(3,3),cc_bohr(3,3),cc_velo(3,3),
     .                 coord(3,nat),veloc(3,nat),coort(3),
     .                 obox(3),rbox(3,3),rinv(3,3),b2ang
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      logical varcel
      character symbol*2
      external fillbox

      rewind (ii2) 
      rewind (io1) 
C --- write selected atoms first into a scratch file (is1), for the case
C     there are zero. Then the label 'ATOMS' with no  atoms following
C     will crash XCrySDen.
      open (is1,form='formatted',status='scratch')
C
      write (io1,"('ANIMSTEPS',i6)") (mdlast-mdfirst)/mdstep+1
      write (io1,"('CRYSTAL')")
      write (io1,"('PRIMVEC')")
      do ii=1,3
        write (io1,'(3x,3f18.9)') (rbox(jj,ii),jj=1,3)
      enddo
      if (mdmod.eq.1) then  ! (read data from MD, all in Bohr)
        do istep = 1,nstep
          read (ii2) idum,coord,veloc
          if (varcel) read (ii2) cc_bohr,cc_velo
          if ((istep.ge.mdfirst).and.
     .        (istep.le.mdlast).and.
     .        (mod(istep-mdfirst,mdstep).eq.0)) then
            if (varcel) then
C         ... with cc_bohr read from actual MD step and converted to Ang:
              call fillbox(is1,obox,rbox,rinv,cc_bohr*b2ang,nat,
     .                     coord*b2ang,nbox)
            else
C         ... with cc_ang as read from the XV file
              call fillbox(is1,obox,rbox,rinv,cc_ang,nat,
     .                     coord*b2ang,nbox)
            endif  ! if (varcel)
            if (nbox.gt.0) then
              write (io1,"('PRIMCOORD',i6)") (istep-mdfirst)/mdstep+1
              write (io1,"(i6,'  1')") nbox
              rewind (is1)
              do ibox = 1,nbox
                read  (is1,201) iat,     (coort(jj),jj=1,3)
                write (io1,201) nz(iat), (coort(jj)-obox(jj),jj=1,3)
              enddo
            else
              write (6,*) ' Empty box for MD step ',istep,' !'
            endif
          endif  !  if write this MD step
        enddo  !  do istep
      elseif (mdmod.eq.2) then  ! (read data from ANI, coord. in Ang)
        do istep = 1,nstep
          read (ii2,301) nat
          do iat=1,nat
            read (ii2,'(a2,2x,3f12.6)') symbol,(coord(ii,iat),ii=1,3)
          enddo
          if ((istep.ge.mdfirst).and.
     .        (istep.le.mdlast).and.
     .        (mod(istep-mdfirst,mdstep).eq.0)) then
C           coordinates in ANI are in Angstroem, and so passed to fillbox
            call fillbox(is1,obox,rbox,rinv,cc_ang,nat,coord,nbox)
            if (nbox.gt.0) then
              write (io1,"('PRIMCOORD',i6)") (istep-mdfirst)/mdstep+1
              write (io1,"(i6,'  1')") nbox
              rewind (is1)
              do ibox = 1,nbox
                read  (is1,201) iat,     (coort(jj),jj=1,3)
                write (io1,201) nz(iat), (coort(jj)-obox(jj),jj=1,3)
              enddo
            else
              write (6,*) ' Empty box for MD step ',istep,' !'
            endif
          endif  !  if write this MD step
        enddo  !  do istep
      endif
      return
  201 format (i4,3f20.8)
  301 format(i5,/)

      end
