c $Id: iomd.f,v 1.4 1999/11/26 18:28:17 wdpgaara Exp $

      subroutine iomd( na, isa, iza, xa, va, cell, vcell, varcel,
     .                 istep, istep0, istepf, temp, eks, getot)

c *******************************************************************
c Saves positions, cell, and energies in a MD run (accumulative)
c J.Kohanoff August 1998, slightly modified by E. Artacho, Feb. 1999
c *******************************************************************
c Two possibilities: |*everything into one single unformatted file
c                    | (for postprocessing programs, space saving)
c (formtt parameter) |*separated ascii files
c ********** INPUT **************************************************
c real*8  cell(3,3)  : Unit cell vectors
c real*8  vcell(3,3) : Velocities thereof
c integer na         : Number of atoms
c integer isa(na)    : Atomic species index
c integer iza(na)    : Atomic numbers
c real*8  xa(3,na)   : Atomic positions
c real*8  va(3,na)   : Velocities thereof
c logical varcel     : .true. when variable cell
c real*8  temp       : Temperature of ions
c real*8  eks        : Kohn-Sham energy
c real*8  getot      : Total energy
c integer istep      : Present time step
c integer istep0     : First time step
c integer istepf     : Last time step
c *******************************************************************

      use fdf

      implicit          none
      character         paste*33
      integer           na, isa(na), iza(na)
      integer           istep, istep0, istepf
      logical           varcel
      double precision  cell(3,3), xa(3,na), va(3,na), vcell(3,3),
     .                  temp, eks, getot
      external          io_assign, io_close, paste

c Internal variables and arrays
      logical    formtt
      parameter  ( formtt = .false. )

      character  sname*30, fnpos*33, fncel*33, fnene*33
      integer    ia, iupos, iuene, iucel, iv, ix
      logical    frstme, formt
      save       frstme, formt, fnpos, fnene, 
     .           fncel, iupos, iuene, iucel

      data frstme /.true./


c Find name of file
      if (frstme) then
        formt = formtt
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fnene = paste( sname, '.MDE' )
        if (formt) then
          fnpos = paste( sname, '.MDX' )
          if (varcel) fncel = paste( sname, '.MDC' )
        else
          fnpos = paste( sname, '.MD' )
        endif
        frstme = .false.
      endif

c Open file (only the first MD step)

      if(istep . eq . istep0) then
        call io_assign( iuene )
        open(iuene, file=fnene, form='formatted', position='append', 
     .    status='unknown')
        write(iuene,"(4a,/)") 'Step','   Temperature',
     .                      '     KS energy','  Total energy'
        if ( formt ) then
          call io_assign( iupos )
          open(iupos, file=fnpos, form='formatted', position='append',
     .      status='unknown')
          if ( varcel ) then 
            call io_assign( iucel )
            open(iucel,file=fncel,form='formatted',position='append',
     .        status='unknown' )
          endif
        else
          call io_assign( iupos )
          open(iupos,file=fnpos,form='unformatted',status='unknown')
          call windu(iupos)
        endif
      endif

c Write data on files

      write(iuene,'(i4,5(2x,f12.6))') istep, temp, eks, getot
      if ( formt ) then
        write(iupos,*) istep
        do ia = 1,na
          write(iupos,'(i3,i6,3f13.6,3x,3f13.6)')
     .      isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
        enddo
        if ( varcel ) then
          write(iucel,*) istep
          write(iucel,'(2(3x,3f13.6))')
     .      ((cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3),iv=1,3)
        endif
      else
        write(iupos) istep, xa, va
        if ( varcel ) write(iupos) cell, vcell
      endif

c Close file (only at the end of the MD run)

      if(istep. eq. istepf) then
        call io_close( iuene )
        call io_close( iupos )
        if ( formt .and. varcel ) call io_close( iucel )
      endif

      return
      end


c *******************************************************************
c Wind to end of file of an unformatted file
      subroutine windu(iu)
      integer iu

 1    read(iu,end=2)
      goto 1
 2    continue
      end
