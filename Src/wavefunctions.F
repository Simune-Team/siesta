 
       module wavefunctions

       use precision

       implicit none 
     
       type  occ_wavef
          integer, dimension (:,:), pointer :: nocck !Number of occupied wavef. 
                                                     !per kpoint and spin
          integer, dimension (:,:,:), pointer :: iwf !wavefunction index
          double precision, dimension (:,:,:), pointer :: phi !wavefunctions
       end type  occ_wavef

       type(occ_wavef), save  :: wavef, init
 
       CONTAINS

       subroutine alloc_wavef(func,nuotot,nocc,nk,nspin)

       type (occ_wavef), intent(inout)  :: func
       integer, intent(in)  :: nuotot
       integer, intent(in)  :: nk
       integer, intent(in)  :: nspin
       integer, intent(in)  :: nocc(:,:)

       integer   nwf, ik, ispin, mxnwf

       allocate(func%nocck(nk,nspin))
             
       nwf=0
       mxnwf=0
       do ispin=1,nspin
          do ik=1,nk
            nwf=nwf+nocc(ik,ispin)
            mxnwf=max(mxnwf,nocc(ik,ispin))
          enddo 
       enddo 
       
       allocate(func%iwf(mxnwf,nk,nspin))
       allocate(func%phi(2,nuotot,nwf))
       
       end subroutine alloc_wavef
   
       subroutine iowavef(task,func,nuotot,nk,nspin)

       use parallel
       use fdf
#ifdef MPI
       use mpi
#endif
       use elec_dyn_options

       type (occ_wavef), intent(inout)  :: func
       integer, intent(in)  :: nuotot
       integer, intent(in)  :: nk
       integer, intent(in)  :: nspin
       character*(*), intent(in):: task

C Internal variables and arrays
      character fname*33, sname*30, paste*33
      logical   exist1, frstme
c     integer   Node, Nodes, unit1, ie, 
      integer   unit1, ie, 
     .          nuototread,nkread,nspinread
      integer   nwf, ik, ispin, mxnwf, io, ix
      external  chkdim, io_assign, io_close, paste, timer,
     .          memory

      save      frstme, fname
      data      frstme /.true./


c Get the Node number

#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Find file name
      if (frstme) then
        if (Node.eq.0) then
          sname = fdf_string('SystemLabel','siesta')
        endif
#ifdef MPI
        call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World,
     .    MPIerror)
#endif
        fname = paste(sname,'.TDWF')
        frstme = .false.
      endif

#ifdef MPI
       stop 'iowavef: This subroutine not yet prepared for parallel run'
#endif

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file=fname, exist=exist1)
          if(.not.exist1) then 
           write(6,'(/,a)')
     .      'iowavef: file containing the KS orbitals not found'
           write(6,'(a)')
     .      'iowavef: simulation can not be restarted'
             stop 
          endif
          call io_assign(unit1)
          open( unit1, file=fname, form='unformatted',
     .      status='unknown' )
          read(unit1) restarttime
          read(unit1) nuototread, nkread, nspinread
          if(nkread.ne.nk) stop 'iowavef: nk inconsistent '
          if(nuototread.ne.nuotot) stop 'iowavef: nuotot inconsistent ' 
          if(nspinread.ne.nspin) stop 'iowavef: nspin inconsistent '

          allocate(func%nocck(nk,nspin))

          do ispin=1,nspin
              read(unit1) (func%nocck(ik,ispin),ik=1,nk)
          enddo

          nwf=0   
          mxnwf=0
          do ispin=1,nspin
            do ik=1,nk
              nwf=nwf+func%nocck(ik,ispin)
              mxnwf=max(mxnwf,func%nocck(ik,ispin))
            enddo
          enddo

          allocate(func%iwf(mxnwf,nk,nspin))
          allocate(func%phi(2,nuotot,nwf))

           
          do ispin=1,nspin
           do ik=1,nk
             do ie=1,func%nocck(ik,ispin)
                read(unit1)  nwf
                func%iwf(ie,ik,ispin)=nwf
                read(unit1)
     .     ((func%phi(ix,io,nwf),ix=1,2), io=1,nuotot)
             enddo
           enddo
         enddo
         call io_close(unit1)
        endif




       elseif(task.eq.'write'.or.task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file=fname, form='unformatted',
     .      status='unknown' )
          rewind(unit1)
          write(unit1) tottime(itimestep,initstep,timestep,restarttime)
          write(unit1) nuotot,nk,nspin
          do ispin=1,nspin
              write(unit1) (func%nocck(ik,ispin), ik=1,nk)
          enddo   
          do ispin=1,nspin
           do ik=1,nk
             do ie=1,func%nocck(ik,ispin)
                nwf=func%iwf(ie,ik,ispin)
                write(unit1)  nwf
                write(unit1)  
     .     ((func%phi(ix,io,nwf),ix=1,2), io=1,nuotot) 
             enddo
           enddo 
         enddo 
         call io_close(unit1)
        endif
          
       endif  

       end subroutine iowavef
 

       end module wavefunctions       

