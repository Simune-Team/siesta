      subroutine setatomnodes(na,lasto,Node,Nodes)
C
C  Sets up the distribution of atoms over nodes for a 
C  spatial decomposition in parallel.
C
C  On entry :
C
C  na            = total number of atoms
C  lasto         = pointer to atom orbitals
C  Node          = local node number
C  Nodes         = total number of nodes
C
C  On exit (in modules) :
C
C  ncellpernode               = number of cells on the current node
C  ncellnodeptr(ncellpernode) = pointer to cell indices for each local cell
C  lbuffercell(ncellpernode)  = if .true. this is a buffer region cell
C  nspmin(3)                  = lower bound to real cells (not buffers)
C  nspmax(3)                  = upper bound to real cells (not buffers)
C  natomsG2L(na)              = pointer from global atom no. to local one
C  natomsNode(na)             = pointer to node responsible for each atom
C
C  Julian Gale, NRI, Curtin University, March 2004
C
      use precision
      use spatial
      use alloc
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use mpi_siesta, only: mpi_integer, mpi_max
#endif
      implicit none
C
C  Passed variables
C
      integer, intent(in)                      :: na
      integer, intent(in)                      :: lasto(0:na)
      integer, intent(in)                      :: Node
      integer, intent(in)                      :: Nodes
C
C  Local variables
C
      integer                                  :: i
      integer                                  :: ia
      integer                                  :: icx
      integer                                  :: icy
      integer                                  :: icz
      integer                                  :: ii
      integer                                  :: ind
      integer                                  :: ix
      integer                                  :: iy
      integer                                  :: iz
      integer                                  :: j
      integer                                  :: maxorb
      integer                                  :: maxxy
      integer                                  :: maxx
      integer                                  :: n
      integer                                  :: n2
      integer                                  :: n3
      integer                                  :: n5
      integer                                  :: n235(3)
      integer                                  :: nn235(3)
      integer                                  :: nadd
      integer                                  :: natomnow
      integer                                  :: natomold
      integer, dimension(:), allocatable, save :: npgridxptr
      integer, dimension(:), allocatable, save :: npgridyptr
      integer, dimension(:), allocatable, save :: npgridzptr
      integer                                  :: nspdiff
      integer                                  :: nspnobuff
      logical                                  :: ldebug
      logical                                  :: lfixx
      logical                                  :: lfixy
      logical                                  :: lfixz
      logical                                  :: lnon235
      logical                                  :: lok
      logical                                  :: lokx
      logical                                  :: loky
      logical                                  :: lokz
      real(dp)                                 :: diffxy
      real(dp)                                 :: diffxz
      real(dp)                                 :: diffyz
      real(dp)                                 :: fdiff
      real(dp)                                 :: pratio
      real(dp)                                 :: ratiox
      real(dp)                                 :: ratioy
      real(dp)                                 :: ratioz
      real(dp)                                 :: ratioxy
      real(dp)                                 :: ratioxz
      real(dp)                                 :: ratioyz
      real(dp)                                 :: rnow
      real(dp)                                 :: spcelltot
      real(dp)                                 :: targetx
      real(dp)                                 :: targety
      real(dp)                                 :: targetz
      real(dp)                                 :: targetratioxy
      real(dp)                                 :: targetratioxz
      real(dp)                                 :: targetratioyz
#ifdef MPI
      integer                                  :: MPIerror
      integer, dimension(:), allocatable, save :: ntmp
#endif
C
      if (Nodes.eq.1) then
C**********************
C  Non-parallel case  *
C**********************
        natompernode = na
        if (natompernode.gt.maxatompernode) then
          maxatompernode = natompernode
          call re_alloc(natomsL2G,1,maxatompernode,
     .      name='natomsL2G')
        endif
        do i = 1,natompernode
           natomsL2G(i) = i
        enddo
C
C  For spatial algorithm set pointers to cells
C
        ncellpernode = nspcell(1)*nspcell(2)*nspcell(3)
        if (ncellpernode.gt.maxcellpernode) then
          maxcellpernode = ncellpernode
          call re_alloc(lbuffercell,1,maxcellpernode,
     .      name='lbuffercell')
          call re_alloc(ncellnodeptr,1,maxcellpernode,
     .      name='ncellnodeptr')
        endif
        nspmax(1) = nspcell(1) - nbufferx
        nspmin(1) = nbufferx
        nspmax(2) = nspcell(2) - nbuffery
        nspmin(2) = nbuffery
        nspmax(3) = nspcell(3) - nbufferz
        nspmin(3) = nbufferz
        maxxy = nspcell(1)*nspcell(2)
        maxx  = nspcell(1)
        n = 0
        do iz = 1,nspcell(3)
          lokz = (iz.gt.nbufferz.and.iz.lt.(nspcell(3)-nbufferz+1))
          do iy = 1,nspcell(2)
            loky = (iy.gt.nbuffery.and.iy.lt.(nspcell(2)-nbuffery+1))
            do ix = 1,nspcell(1)
              lokx = (ix.gt.nbufferx.and.ix.lt.(nspcell(1)-nbufferx+1))
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              n = n + 1
              ncellnodeptr(n) = ind
              if (lokx.and.loky.and.lokz) then
                lbuffercell(n) = .false.
              else
                lbuffercell(n) = .true.
              endif
            enddo
          enddo
        enddo
C
C  Build lists linking atoms with parallel structure
C
        call re_alloc(natomsG2L,1,na,name='natomsG2L')
        call re_alloc(natomsNode,1,na,name='natomsNode')
        natomsNode(1:na) = 0
        do ii = 1,na
          natomsG2L(ii) = ii
        enddo
C
C  Build dummy lists linking orbitals with parallel structure 
C
        maxorb = lasto(na)
        call re_alloc(nL2G,1,maxorb,1,Nodes,name='nL2G')
        call re_alloc(nG2L,1,maxorb,name='nG2L')
        call re_alloc(nNode,1,maxorb,name='nNode')
        call re_alloc(nOrbPerNode,1,Nodes,name='nOrbPerNode')
C
        nOrbPerNode(1) = maxorb
        do ii = 1,maxorb
          nL2G(ii,1) = ii
          nG2L(ii) = ii
          nNode(ii) = 0
        enddo
      else
C**************************
C  Spatial parallel case  *
C**************************
C
C  Set debugging flag
C
c        ldebug = .false.
        ldebug = (Node.eq.0)
C
C  Find factors of number of processors (2,3,5)
C
        n2 = 0
        n3 = 0
        n5 = 0
        n = Nodes
        lok = .true.
        do while (lok.and.n.gt.1)
          if (mod(n,2).eq.0) then
            n2 = n2 + 1
            n = n/2
          elseif (mod(n,3).eq.0) then
            n3 = n3 + 1
            n = n/3
          elseif (mod(n,5).eq.0) then
            n5 = n5 + 1
            n = n/5
          else
            lok = .false.
          endif
        enddo
C
C  Set flag to indicate whether there is a factor other than 2,3 or 5
C
        lnon235 = (n.ne.1)
        if (ldebug) then
          write(6,'(/,''  Spatial decomposition of processors :'',/)')
          write(6,'(''  No. of multiples of 2 = '',i4)') n2
          write(6,'(''  No. of multiples of 3 = '',i4)') n3
          write(6,'(''  No. of multiples of 5 = '',i4)') n5
          if (lnon235) then
            write(6,'(''  Other factor          = '',i4)') n
          endif
        endif
C
C  If non 2,3,5, decide which direction to fix
C
        lfixx = .false.
        lfixy = .false.
        lfixz = .false.
        if (lnon235) then
          spcelltot = dble(nspcell(1)-2*nbufferx)*
     .                dble(nspcell(2)-2*nbuffery)*
     .                dble(nspcell(3)-2*nbufferz)
          pratio = dble(n)/dble(Nodes)
          ratiox = pratio - dble(nspcell(1)-2*nbufferx)/spcelltot 
          ratioy = pratio - dble(nspcell(2)-2*nbuffery)/spcelltot
          ratioz = pratio - dble(nspcell(3)-2*nbufferz)/spcelltot
          ratiox = abs(ratiox)
          ratioy = abs(ratioy)
          ratioz = abs(ratioz)
          if (ratiox.lt.ratioy.and.ratiox.lt.ratioz) then
            lfixx = .true.
            npgrid(1) = n
          elseif (ratioy.lt.ratioz) then
            lfixy = .true.
            npgrid(2) = n
          else
            lfixz = .true.
            npgrid(3) = n
          endif
        endif
C
C  Calculate optimal target ratios for x-y and x-z
C
        targetratioxy = dble(nspcell(2)-2*nbuffery)/
     .    dble(nspcell(1)-2*nbufferx)
        targetratioxz = dble(nspcell(3)-2*nbufferz)/
     .    dble(nspcell(1)-2*nbufferx)
        targetratioyz = dble(nspcell(3)-2*nbufferz)/
     .    dble(nspcell(2)-2*nbuffery)
C
C  Group factors in proportion to numbers of cells
C
        npgrid(1:3) = 1
        n235(1) = n5
        n235(2) = n3
        n235(3) = n2
        nn235(1) = 5
        nn235(2) = 3
        nn235(3) = 2
        do ind = 1,3
          ii = n235(ind)
          n = nn235(ind)
          do i = 1,ii
            if (lnon235) then
              if (lfixx) then
                ratioyz = dble(npgrid(3))/dble(npgrid(2))
                diffyz = targetratioyz - ratioyz
                if (diffyz.lt.0.0_dp) then
                  npgrid(2) = n*npgrid(2)
                else
                  npgrid(3) = n*npgrid(3)
                endif
              elseif (lfixy) then
                ratioxz = dble(npgrid(3))/dble(npgrid(1))
                diffxz = targetratioxz - ratioxz
                if (diffxz.lt.0.0_dp) then
                  npgrid(1) = n*npgrid(1)
                else
                  npgrid(3) = n*npgrid(3)
                endif
              else
                ratioxy = dble(npgrid(2))/dble(npgrid(1))
                diffxy = targetratioxy - ratioxy
                if (diffxy.lt.0.0_dp) then
                  npgrid(1) = n*npgrid(1)
                else
                  npgrid(2) = n*npgrid(2)
                endif
              endif
            else
              ratioxy = dble(npgrid(2))/dble(npgrid(1))
              ratioxz = dble(npgrid(3))/dble(npgrid(1))
              diffxy = targetratioxy - ratioxy
              diffxz = targetratioxz - ratioxz
              if (diffxy.lt.0.0_dp.and.diffxz.lt.0.0_dp) then
                npgrid(1) = n*npgrid(1)
              elseif (diffxy.gt.diffxz) then
                npgrid(2) = n*npgrid(2)
              else
                npgrid(3) = n*npgrid(3)
              endif
            endif
          enddo
        enddo
        if (ldebug) then
          write(6,'(/,''  Processor grid = '',3(i4,1x))') 
     .      (npgrid(n),n=1,3)
        endif
C
C  Calculate target number of atoms per node in each direction
C
        targetx = dble(na)/dble(npgrid(1))
        targety = dble(na)/dble(npgrid(2))
        targetz = dble(na)/dble(npgrid(3))
C
C  Allocate local memory
C
        allocate(npgridxptr(0:npgrid(1)))
        allocate(npgridyptr(0:npgrid(2)))
        allocate(npgridzptr(0:npgrid(3)))
C
        npgridxptr(0) = nbufferx
        npgridyptr(0) = nbuffery
        npgridzptr(0) = nbufferz
C
        maxxy = nspcell(1)*nspcell(2)
        maxx  = nspcell(1)
C
C  Assign cells to nodes based on cumulative numbers of atoms : X
C
        natomnow = 0
        natomold = 0
        npgridxptr(1:npgrid(1)) = 0
        ii = 0
        do ix = nbufferx+1,nspcell(1) - nbufferx
          do iz = nbuffery+1,nspcell(3) - nbuffery
            do iy = nbufferz+1,nspcell(2) - nbufferz
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          rnow = dble(natomnow - natomold)
          fdiff = abs(rnow-targetx)/targetx
          if (rnow.ge.targetx.or.fdiff.lt.0.10d0) then
            ii = ii + 1
            npgridxptr(ii) = ix
            natomold = natomnow
          endif
        enddo
        npgridxptr(npgrid(1)) = nspcell(1) - nbufferx
C
C  Assign cells to nodes based on cumulative numbers of atoms : Y
C
        natomnow = 0
        natomold = 0
        npgridyptr(1:npgrid(2)) = 0
        ii = 0
        do iy = nbuffery+1,nspcell(2) - nbuffery
          do iz = nbufferz+1,nspcell(3) - nbufferz
            do ix = nbufferx+1,nspcell(1) - nbufferx
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          rnow = dble(natomnow - natomold)
          fdiff = abs(rnow-targety)/targety
          if (rnow.ge.targety.or.fdiff.lt.0.10d0) then
            ii = ii + 1
            npgridyptr(ii) = iy
            natomold = natomnow
          endif
        enddo
        npgridyptr(npgrid(2)) = nspcell(2) - nbuffery
C
C  Assign cells to nodes based on cumulative numbers of atoms : Z
C
        natomnow = 0
        natomold = 0
        npgridzptr(0:npgrid(3)) = nbufferz
        ii = 0
        do iz = nbufferz+1,nspcell(3) - nbufferz
          do iy = nbuffery+1,nspcell(2) - nbuffery
            do ix = nbufferx+1,nspcell(1) - nbufferx
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          rnow = dble(natomnow - natomold)
          fdiff = abs(rnow-targetz)/targetz
          if (rnow.ge.targetz.or.fdiff.lt.0.10d0) then
            ii = ii + 1
            npgridzptr(ii) = iz
            natomold = natomnow
          endif
        enddo
        npgridzptr(npgrid(3)) = nspcell(3) - nbufferz
C
C  Build lists linking cells to nodes
C
        n = 0
        ncellpernode = 0
        do iz = 1,npgrid(3)
          do iy = 1,npgrid(2)
            do ix = 1,npgrid(1)
              n = n + 1
              if (Node.eq.n-1) then
                nspmax(1) = npgridxptr(ix) 
                nspmin(1) = npgridxptr(ix-1)
                nspmax(2) = npgridyptr(iy)
                nspmin(2) = npgridyptr(iy-1)
                nspmax(3) = npgridzptr(iz)
                nspmin(3) = npgridzptr(iz-1)
C
                nspnobuff = nspcell(1) - 2*nbufferx
                nspdiff = min(nspnobuff-nspmax(1)+nspmin(1),nbufferx)
                if (nspdiff.gt.0) nspmin(1) = nspmin(1) - nspdiff
                nspdiff = min(nspnobuff-nspmax(1)+nspmin(1),nbufferx)
                if (nspdiff.gt.0) nspmax(1) = nspmax(1) + nspdiff
                nspnobuff = nspcell(2) - 2*nbuffery
                nspdiff = min(nspnobuff-nspmax(2)+nspmin(2),nbuffery)
                if (nspdiff.gt.0) nspmin(2) = nspmin(2) - nspdiff
                nspdiff = min(nspnobuff-nspmax(2)+nspmin(2),nbuffery)
                if (nspdiff.gt.0) nspmax(2) = nspmax(2) + nspdiff
                nspnobuff = nspcell(3) - 2*nbufferz
                nspdiff = min(nspnobuff-nspmax(3)+nspmin(3),nbufferz)
                if (nspdiff.gt.0) nspmin(3) = nspmin(3) - nspdiff
                nspdiff = min(nspnobuff-nspmax(3)+nspmin(3),nbufferz)
                if (nspdiff.gt.0) nspmax(3) = nspmax(3) + nspdiff
C
                nspmax(1) = min(nspmax(1),nspcell(1)-1)
                nspmin(1) = max(nspmin(1),1)
                nspmax(2) = min(nspmax(2),nspcell(2)-1)
                nspmin(2) = max(nspmin(2),1)
                nspmax(3) = min(nspmax(3),nspcell(3)-1)
                nspmin(3) = max(nspmin(3),1)
C
                nadd = (npgridzptr(iz) - npgridzptr(iz-1) + 2*nbufferz)*
     .                 (npgridyptr(iy) - npgridyptr(iy-1) + 2*nbuffery)*
     .                 (npgridxptr(ix) - npgridxptr(ix-1) + 2*nbufferx)
                if (ncellpernode+nadd.gt.maxcellpernode) then
                  maxcellpernode = ncellpernode + nadd
                  call re_alloc(lbuffercell,1,maxcellpernode,
     .              name='lbuffercell')
                  call re_alloc(ncellnodeptr,1,maxcellpernode,
     .              name='ncellnodeptr')
                endif
C
                do icz = npgridzptr(iz-1)+1-nbufferz,
     .              npgridzptr(iz)+nbufferz
                  lokz = (icz.gt.npgridzptr(iz-1).and.
     .              icz.le.npgridzptr(iz))
                  do icy = npgridyptr(iy-1)+1-nbuffery,
     .                npgridyptr(iy)+nbuffery
                    loky = (icy.gt.npgridyptr(iy-1).and.
     .                icy.le.npgridyptr(iy))
                    do icx = npgridxptr(ix-1)+1-nbufferx,
     .                  npgridxptr(ix)+nbufferx
                      lokx = (icx.gt.npgridxptr(ix-1).and.
     .                  icx.le.npgridxptr(ix))
                      ind = (icz-1)*maxxy + (icy-1)*maxx + icx
                      ncellpernode = ncellpernode + 1
                      ncellnodeptr(ncellpernode) = ind
                      if (lokx.and.loky.and.lokz) then
                        lbuffercell(ncellpernode) = .false.
                      else
                        lbuffercell(ncellpernode) = .true.
                      endif
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
        if (ldebug) then
          write(6,'(''  Cells per Processor = '',2(i8,1x))') 
     .      Node,ncellpernode
        endif
C
C  Find size for atom link array and then allocate
C
        natompernode = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            natompernode = natompernode + nspcellat(ind)
          endif
        enddo
        if (natompernode.gt.maxatompernode) then
          maxatompernode = natompernode
          call re_alloc(natomsL2G,1,maxatompernode,
     .      name='natomsL2G')
        endif
C
C  Build lists linking atoms to nodes
C
        natompernode = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              natompernode = natompernode + 1
              natomsL2G(natompernode) = nspcellatptr(n,ind)
            enddo
          endif
        enddo
        if (ldebug) then
          write(6,'(''  Atoms per Processor = '',2(i8,1x))') 
     .      Node,natompernode
        endif
C
C  Build lists linking atoms with parallel structure
C
        call re_alloc(natomsG2L,1,na,name='natomsG2L')
        call re_alloc(natomsNode,1,na,name='natomsNode')
        natomsNode(1:na) = 0
        natomsG2L(1:na) = 0
        j = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              ia = nspcellatptr(n,ind)
              j = j + 1
              natomsG2L(ia) = j
              natomsNode(ia) = Node
            enddo
          endif
        enddo
C
C  Globalise data
C
        maxorb = lasto(na)
#ifdef MPI
        allocate(ntmp(max(maxorb,Nodes)))
        ntmp(1:na) = natomsNode(1:na)
        call MPI_AllReduce(ntmp,natomsNode,na,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        ntmp(1:na) = natomsG2L(1:na)
        call MPI_AllReduce(ntmp,natomsG2L,na,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
#endif
C
C  Build lists linking orbitals with parallel structure 
C
        call re_alloc(nL2G,1,maxorb,1,Nodes,name='nL2G')
        call re_alloc(nG2L,1,maxorb,name='nG2L')
        call re_alloc(nNode,1,maxorb,name='nNode')
        call re_alloc(nOrbPerNode,1,Nodes,name='nOrbPerNode')
C
        nOrbPerNode(1:Nodes) = 0
        nNode(1:maxorb) = 0
        nL2G(1:maxorb,1:Nodes) = 0
        nG2L(1:maxorb) = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              ia = nspcellatptr(n,ind)
              do j = lasto(ia-1)+1,lasto(ia)
                nOrbPerNode(Node+1) = nOrbPerNode(Node+1) + 1
                nL2G(nOrbPerNode(Node+1),Node+1) = j
                nG2L(j) = nOrbPerNode(Node+1)
                nNode(j) = Node
              enddo         
            enddo         
          endif
        enddo
C
C  Globalise data
C
#ifdef MPI
        ntmp(1:maxorb) = nNode(1:maxorb)
        call MPI_AllReduce(ntmp,nNode,maxorb,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        do ii = 1,Nodes
          ntmp(1:maxorb) = nL2G(1:maxorb,ii)
          call MPI_AllReduce(ntmp,nL2G(1,ii),maxorb,MPI_integer,
     .      MPI_max,MPI_Comm_World,MPIerror)
        enddo
        ntmp(1:maxorb) = nG2L(1:maxorb)
        call MPI_AllReduce(ntmp,nG2L,maxorb,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        ntmp(1:Nodes) = nOrbPerNode(1:Nodes)
        call MPI_AllReduce(ntmp,nOrbPerNode,Nodes,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        deallocate(ntmp)
#endif
C
C  Free local memory
C
        deallocate(npgridxptr)
        deallocate(npgridyptr)
        deallocate(npgridzptr)
      endif
C
      return
      end
