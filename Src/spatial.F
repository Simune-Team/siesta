      module spatial

C
C  This module holds most of the data relating to the spatial decomposition of the system
C
C  Julian Gale, NRI, March 2004
C
        use precision
        implicit none

        logical                                 :: lspatial

        integer,                           save :: maxatompernode = 0
        integer,                           save :: maxcellpernode = 0
        integer,                           save :: maxspcell = 0
        integer,                           save :: maxspcellat = 0
        integer,                           save :: natompernode
        integer,                           save :: nbufferx = 0
        integer,                           save :: nbuffery = 0
        integer,                           save :: nbufferz = 0
        integer,                           save :: ncellpernode
        real(dp),                          save :: rspatial
C
C  Atom pointers
C
C  natomsNode -> pointer from atom to Node number
C  natomsL2G  -> pointer from local atom number to global atom number
C  natomsG2L  -> pointer from global atom number to local atom number on relevant node
C
        integer,  dimension(:),   pointer, save :: natomsNode
        integer,  dimension(:),   pointer, save :: natomsL2G
        integer,  dimension(:),   pointer, save :: natomsG2L
C
        integer,  dimension(:),   pointer, save :: ncellnodeptr
        integer,  dimension(:),   pointer, save :: nG2L
        integer,  dimension(:,:), pointer, save :: nL2G
        integer,  dimension(:),   pointer, save :: nNode
        integer,  dimension(:),   pointer, save :: nOrbPerNode
        integer,                           save :: npgrid(3)
        integer,                           save :: nspcell(3)
        integer,  dimension(:),   pointer, save :: nspcellat
        integer,  dimension(:,:), pointer, save :: nspcellatptr
        integer,  dimension(:,:), pointer, save :: nspcellatptrcell
        integer,                           save :: nspmax(3)
        integer,                           save :: nspmin(3)
        logical,  dimension(:),   pointer, save :: lbuffercell
        real(dp),                          save :: spcell(3)
        real(dp),                          save :: spmin(3)
        real(dp), dimension(:),   pointer, save :: xinbox
        real(dp), dimension(:),   pointer, save :: yinbox
        real(dp), dimension(:),   pointer, save :: zinbox

      end module spatial
