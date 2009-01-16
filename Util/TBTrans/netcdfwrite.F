#ifndef CDF

      subroutine netcdf_write_teig(closeme,neig,Teig1,En1,ispin,nspin)
c INPUT
      logical closeme
      integer neig,ispin,nspin
      real*8 En1,Teig1(*)      
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#else ! Netcdf defined 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine netcdf_write_teig(closeme,neig,Teig1,En1,ispin,nspin)
      
c      implicit none

c      include 'netcdf.inc'
      use netcdf

c INPUT
      logical closeme
      integer neig,ispin,nspin
      real*8 En1,Teig1(neig)

c LOCAL  
      logical firsttime
      data firsttime /.true./
      save firsttime
      integer icall
      save icall

cccccccccccccccccccccccc BEGIN NETCDF STUFF ccccccccccccccccccccccccccc      
* error status return
      integer  iret
* netCDF id
      integer  ncid
* dimension ids
      integer  NEn_dim
      integer  neigch_dim
* dimension lengths
c      integer  NEn_len
c      integer  neigch_len
c      parameter (NEn_len = NF90_UNLIMITED)
c      parameter (neigch_len = neig)
* variable ids
      integer  Teig_id
      integer  En_id
* rank (number of dimensions) for each variable
      integer  Teig_rank
      integer  En_rank
      parameter (Teig_rank = 2)
      parameter (En_rank = 1)
* variable shapes
      integer  Teig_dims(Teig_rank)
      integer  En_dims(En_rank)
      save ncid,Teig_id, En_id
* data variables
* attribute vectors
* enter define mode
cccccccccccccccccccccccc END NETCDF STUFF ccccccccccccccccccccccccccc


c      neigch_len = neig
      if(closeme) then
         write(6,*) "Done writing netcdf file: Closing"
         iret = nf90_close(ncid)
         call check_err(iret)
         firsttime=.true.
         return
      end if
         

      
      if(firsttime) then        ! open nc file and make it ready
         icall=0
         if(nspin.eq.1) 
     .        iret = nf90_create('teig.nc', NF90_CLOBBER, ncid)
         if(nspin.eq.2) then
            if(ispin.eq.1) 
     .        iret = nf90_create('teig.UP.nc', NF90_CLOBBER, ncid)
            if(ispin.eq.2) 
     .        iret = nf90_create('teig.DOWN.nc', NF90_CLOBBER, ncid)
         endif
         call check_err(iret)
*     define dimensions
         iret = nf90_def_dim(ncid, 'NEn', NF90_UNLIMITED, NEn_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'neigch', neig, neigch_dim)
         call check_err(iret)
*     define variables
         Teig_dims(2) = NEn_dim
         Teig_dims(1) = neigch_dim
         iret = nf90_def_var(ncid, 'Teig', NF90_REAL, 
     .        Teig_dims, Teig_id)
         call check_err(iret)
         En_dims(1) = NEn_dim
         iret = nf90_def_var(ncid, 'En', NF90_REAL, 
     .        En_dims, En_id)
         call check_err(iret)
*     assign attributes
         iret = nf90_put_att(ncid, Teig_id, 
     .        'long_name', 'Eigenchannel transmissions')
         call check_err(iret)
         iret = nf90_put_att(ncid, En_id, 'long_name','Energy')
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 
     .        'source', 'TBTrans')
         call check_err(iret)
*     leave define mode
         iret = nf90_enddef(ncid)
         call check_err(iret)
         firsttime=.false.
      end if                    !first time

cccccccccccccccccccccccccccccccccccccccccccccccc
       

*     Write record variables
      icall=icall+1
c      write(6,*) icall
      iret = nf90_put_var(ncid,Teig_id, Teig1(1:),start=(/1,icall/))
      call check_err(iret)
      
      iret = nf90_put_var(ncid,En_id,En1,start=(/icall/))
      call check_err(iret)
      

      CONTAINS
      
      
      subroutine check_err(status)
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'netcdf Error: Stopped'
      endif
      end subroutine check_err

      
      end subroutine netcdf_write_teig             


#endif
