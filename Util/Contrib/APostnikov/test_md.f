C...............................................................
C
      subroutine test_md(ii2,nat,varcel,mdstep)        
C
C     reads from ii2 (MD file), makes some consistency tests,
C     finds out whether the records are for fixed cell or variable cell,
C     returns the number of MD steps.
C
      implicit none
      integer ii2,nat,mdstep,istep,iprev,ii
      logical varcel
      double precision dummy
C     double precision cc_bohr(3,3),cc_velo(3,3),
C    .                 coord(3,nat),veloc(3,nat)

      varcel = .false.
       write (6,*) ' Assume fixed cell...'
  101 continue
      rewind (ii2)
C --- read through input file, to get the total number of MD steps.
C     It must be specified at the beginning of AXSF file.
C       read from unformatted MD, as written in subr. iomd (coord. in Bohr):
        mdstep = 1
  104  continue
C ...   variable format ( formt=F in subr. iomd):
C       read (ii2,err=807,end=107) istep,coord,veloc
        read (ii2,err=807,end=107) istep,(dummy,ii=1,3*nat),
     .                                   (dummy,ii=1,3*nat)
        write (6,*) ' MD record No. ',mdstep,',  istep=',istep
        if (mdstep.gt.1.and.istep.ne.iprev+1) then
C         we check that the sequential number of MD step is increased by 1,
C         hence the structure of record is correct. If this is not the case,
C         try to switch from fixed cell to variable one.
          write (6,*) ' Oops... Try variable cell'
          varcel = .true.
          goto 101
        endif
C       if (varcel) read (ii2,err=807,end=108) cc_bohr,cc_velo
        if (varcel) read (ii2,err=807,end=108) (dummy,ii=1,9),
     .                                         (dummy,ii=1,9)
        mdstep = mdstep + 1
        iprev = istep
        goto 104
  807 continue
C --- read error, presumably due to a wrong guess of varcel.
      if (.not.varcel) then
        write (6,202) mdstep
        varcel = .true.
        goto 101
      else
        write (6,203) mdstep
        stop
      endif
  108 continue
C --- irregular end of records in MD file: 
      write (6,*) ' Uncomplete record in MD step No.',mdstep
      mdstep = mdstep - 1  !  No. of full records
      if (mdstep.gt.0) then
        write (6,*) ' Keep ',mdstep,' records.'
        return
      else
        write (6,*) ' Check the MD file; bye'
        stop
      endif
  107 continue
C --- regular end of records in MD or ANI file: 
      mdstep = mdstep - 1  !  Attempt to read record Nr. mdstep failed
      if (mdstep.gt.0) then
        write (6,*) '  Cleanly read in ',mdstep,'  MD steps'
        return
      else
        write (6,*) ' End of record in the firt step: MD file empty?'
        stop
      endif
      return
      
  201 format(' MD record No. ',i6,'  istep =',i6)
  202 format(' ...results in read error for step=',i6,
     .       ' This is not yet so bad; I check now whether MD was ',
     .       ' for variable cell...')
  203 format(' ...and still get read error for step=',i6,
     .       ' It seems that the MD file is either corrupted,',
     .       ' or not compatible with XV.')
      end 
C
C...............................................................
C
      subroutine test_ani(ii2,nat,mdstep)        
C
C     reads from ii2 (ANI file), makes some consistency tests,
C     returns the number of MD steps.
C
      implicit none
      integer ii2,nat,mdstep,na,iat,ii
      logical varcel
      double precision dummy
C     double precision coord(3,nat)
      character symbol*2

      rewind (ii2)
C --- read through input file, to get the total number of MD steps.
C     It must be specified at the beginning of AXSF file.
C       read from unformatted MD, as written in subr. pixmol (coord. in Ang):
      mdstep = 1
  101 continue
      read (ii2,301,err=106,end=107) na
      if (na.ne.nat) then
        write (6,*) ' Error reading ANI file, step=',mdstep,
     .              ' : na=',na,' not matching Nr of atoms in XV'
        stop
      endif
      do iat=1,na
C       read (ii2,302,err=108,end=109) symbol,(coord(ii,iat),ii=1,3)
        read (ii2,302,err=108,end=109) symbol,(dummy,ii=1,3)
      enddo
      mdstep = mdstep + 1
      goto 101
  106 continue
C --- Error reading ANI file: 
      write (6,*) ' Error reading header of record ',mdstep,
     .            ' in the ANI file.'
      mdstep = mdstep - 1  !  No. of full records
      if (mdstep.gt.0) then
        write (6,*) ' Keep ',mdstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif
  107 continue
C --- regular end of records ANI file: 
      mdstep = mdstep - 1  !  Attempt to read record Nr. mdstep failed
      write (6,*) ' Cleanly read in ',mdstep,'  MD steps'
      return
  108 continue
C --- Error reading ANI file: 
      write (6,*) ' Error reading coordinates block in record ',
     .              mdstep,', for atom ',iat,' in the ANI file.'
      mdstep = mdstep - 1  !  No. of full records
      if (mdstep.gt.0) then
        write (6,*) ' Keep ',mdstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif
  109 continue
C --- irregular end of records in ANI file: 
      write (6,*) ' Uncomplete record in MD step No.',
     .              mdstep,', for atom ',iat,' in the ANI file.'
      mdstep = mdstep - 1  !  No. of full records
      if (mdstep.gt.0) then
        write (6,*) ' Keep ',mdstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif

  301 format(i5,/)
  302 format(a2,2x,3f12.6)

      end 
