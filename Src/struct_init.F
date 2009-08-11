      module m_struct_init
      private
      public :: struct_init
      
      CONTAINS

      subroutine struct_init()
      USE siesta_options
      use m_ioxv
      USE alloc,          only: re_alloc, de_alloc

      use siesta_geom
      use atmfuncs,       only : izofis
      use atomlist,       only : elem, iza
      use fdf
      use periodic_table, only : symbol
      use iopipes,        only : coordsfrompipe
      use m_iostruct,     only : read_struct
      use siesta_cml
      use files,          only : slabel
      use zmatrix,        only : lUseZmatrix
      use parallel,       only : IOnode

      implicit none

      real(dp), external :: volcel
      external :: automatic_cell, bonds, iozm, shaper
      
      character(len=80) :: dummy_str
      logical           :: Server_calculation, foundxv, foundzm
      real(dp)          :: local_charnet

      real(dp):: bcell(3,3), tmp_cell(3,3)
      integer :: i, ia, ix
      integer:: nbcell  ! Auxiliary to call shaper routine
      real(dp), pointer :: xang(:,:)

!------------------------------------------------------------------------- BEGIN

      !! Read number of atoms and coordinates, and unit cell
      use_struct_file = fdf_get("UseStructFile",.false.)
      use_struct_file = fdf_get("MD.UseStructFile",
     $     use_struct_file)    ! For legacy symbol
      if (use_struct_file) then
         call read_struct( na_u, tmp_cell) ! Sets na_u, xa, and isa
      else
         call coor(na_u,tmp_cell)  ! Sets na_u, xa, and isa
      endif
      ucell = tmp_cell  ! initialize module variable

!     Prepare iza here: it might be needed by ioxv      
      nullify( iza, va )
      call re_alloc( iza, 1, na_u, 'iza', 'struct_init' )
      call re_alloc( va, 1, 3, 1, na_u, 'va', 'struct_init' )

      ! Options readed here instead of siesta_options
      default   = fdf_get('UseSaveData', .false.)
      usesavexv = fdf_get('MD.UseSaveXV', default)
      usesavezm = fdf_get('MD.UseSaveZM', default)
      writic    = fdf_get('WriteCoorInitial', .true.)

      ! Read Z-matrix coordinates and forces from file
      if (lUseZmatrix) then
        foundzm = .false.
        if (usesavezm) then
            call iozm('read',ucell,vcell,xa,foundzm)
            if (IOnode) then
                if (.not.foundzm) then
                   write(6,'(/,a)') 'siesta: WARNING: ZM file not found'
                else
                    write(6,'(/,a)') 
     .        "! Info in XV file prevails over previous structure input"
                endif
            endif
        endif
      endif

! Read cell shape and atomic positions from a former run ..............
      foundxv = .false.
      if (usesavexv) then
        call ioxv('read', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
        if (IOnode) then
           if (.not.foundxv) then
              write(6,'(/,a)') 'siesta: WARNING: XV file not found'
           else
              write(6,"(a)")
     $       "! Info in XV file prevails over previous structure input"
           endif
         endif
      endif
! ..................

! Read cell shape and atomic positions from driver program through pipe
      dummy_str = fdf_get('MD.TypeOfRun',"none")
      Server_calculation = (leqi(dummy_str,'forces'))

      if (Server_calculation) then
        call coordsFromPipe( na_u, xa, ucell )
      end if
! .....................

! Dump initial coordinates to output ..................................
      if ( writic.and.(IOnode) ) then
        write(6,'(/a)') 'siesta: Atomic coordinates (Bohr) and species'
        write(6,"('siesta: ',2x,3f10.5,i3,3x,i6)")
     .           ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na_u)
      endif
! ..................

! Automatic cell generation ...........................................
      if (volcel(ucell) .lt. 1.0d-8) then
         local_charnet = fdf_get('NetCharge',0.0_dp)
         call automatic_cell(ucell,scell,na_u,xa,isa,local_charnet)
      endif

! Initialize atomic velocities to zero ................................
      if (.not. foundxv) then
        ! AG ** What happens with iozm call?
        va(1:3,1:na_u) = 0.0_dp
        vcell(1:3,1:3) = 0.0_dp
      endif
! ..................


! Find system shape ...................................................
      call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
      if (IOnode) then
        write(6,'(/,2a)') 'siesta: System type = ', shape
      endif

! Find interatomic distances (output in file BONDS)
      if (IOnode) then
        rmax_bonds = fdf_physical("MaxBondDistance", 6.0_dp, "Bohr")
        call bonds( ucell, na_u, isa, xa,
     $       rmax_bonds, trim(slabel)// ".BONDS" )
      endif

! Output of initial system details:

      if (cml_p) then
! We need the names of the elements on node 0
        nullify(elem)
        allocate( elem(na_u) )
!        call re_alloc( elem, 1, na_u, 'elem', 'struct_init' )
        do i = 1, na_u
           elem(i) = symbol(izofis(isa(i)))
        enddo

        call cmlStartModule(xf=mainXML, title='Initial System')
        nullify(xang)
        call re_alloc( xang, 1, 3, 1, na_u, 'xang', 'struct_init' )
        xang = xa(1:3,1:na_u)/Ang
        call cmlAddMolecule(xf=mainXML, natoms=na_u, 
     .       coords=xang, elements=elem, refs=cisa,
     .       style='x3', fmt='(f12.6)')
        call de_alloc( xang, 'xang', 'struct_init' )

        call cmlAddLattice(xf=mainXML, cell=ucell, 
     .       units='siestaUnits:angstrom', dictref='siesta:ucell')
        call cmlAddProperty(xf=mainXML, property=trim(shape), 
     .       dictref='siesta:shape')
        call cmlEndModule(xf=mainXML)
      endif


      END subroutine struct_init
      end module m_struct_init
