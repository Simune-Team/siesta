      subroutine dump_atom

      use atmfuncs

      implicit none

      integer is, j, i, lun

      call io_assign(lun)
      open(lun,file='ATOM.dump',form='formatted',status='unknown')
      rewind(lun)
      write(lun,*) 'Number of species (nsmax,ismax): ', nsmax, ismax
      write(lun,*) 'izsave: ', izsave
      write(lun,*) 'nomax: ', nomax
      write(lun,*) 'nkbmax: ', nkbmax
      write(lun,*) 'lmxosave: ', lmxosave
      write(lun,*) 'lmxkbsave: ', lmxkbsave
      write(lun,*) 'semicsave: ', semicsave
      write(lun,*) 'izvaltb: ', izvaltb
      write(lun,*) 'smasstb: ', smasstb
      write(lun,*) 'chargesave: ', chargesave
      write(lun,*) 'slfe: ', slfe
      write(lun,*) 'label_save: ', label_save(1:nsmax)
      write(lun,*) 'basistype_save: ', basistype_save(1:nsmax)

      write(lun,*) '* nkblsave:', nkblsave
      write(lun,*) '* lsemicsave:', lsemicsave
      write(lun,*) '* npolorbsave:', npolorbsave
      write(lun,*) '* nzetasave:', nzetasave
      write(lun,'(a)') '* cnfigtb'
      write(lun,999) cnfigtb
      write(lun,'(a)') '* lambdatb'
      write(lun,999) lambdatb

      write(lun,'(a)') '* qtb'
      write(lun,999) qtb
      write(lun,'(a)') '* qltb'
      write(lun,999) qltb
      write(lun,'(a)') '* table'
      write(lun,999)  table
      write(lun,'(a)') '* tabpol'
      write(lun,999) tabpol
      write(lun,'(a)') '* tab2'
      do is = 1, nsmax
         do j=-nkbmx*(lmaxd+1),nzetmx*nsemx*(lmaxd+1)
            do i = 1, ntbmax
               write(lun,'(3i5,a,g15.5)')
     $                i, j, is, '::', tab2(i,j,is)
            enddo
         enddo
      enddo
ccccc      write(lun,999) tab2
      write(lun,'(a)') '* tab2pol'
      write(lun,999) tab2pol

      write(lun,'(a)') '* coretab'
      write(lun,999) coretab
      write(lun,'(a)') '* chloctab'
      write(lun,999) chloctab
      write(lun,'(a)') '* corrtab'
      write(lun,999) corrtab
      write(lun,'(a)') '* rctb'
      write(lun,999) rctb
      write(lun,'(a)') '* rcotb'
      write(lun,999) rcotb
      write(lun,'(a)') '* rcpoltb'
      write(lun,999) rcpoltb

      call io_close(lun)

 999  format(5g15.5)
      end subroutine dump_atom
