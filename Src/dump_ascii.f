      subroutine dump_ascii

      use types
      use atmfuncs, only:species, nspecies, spp

      implicit none

      integer aux(maxnorbs), lun
      character*20 filename

      integer is, j, i, l

      write(6,'(/,a,a,/)')
     $     '*** Dumping basis and potential info to ASCII files...',
     $     '*** Not quite ready yet...'


      do is = 1, nspecies
         spp => species(is)
         write(filename,'(a,a)') trim(spp%label), ".PAOs"
         call io_assign(lun)
         open(lun,file=filename,status='replace',form='formatted')
!....
         call io_close(lun)

         write(filename,'(a,a)') trim(spp%label), ".KBs"
         call io_assign(lun)
         open(lun,file=filename,status='replace',form='formatted')
!....
         call io_close(lun)

      enddo

      end subroutine dump_ascii








