      subroutine read_basis(nsp,iz,lmxkb,nkbl,erefkb,
     .         lmxo, nzeta, rco, lambda,
     .         atm_label,polorb,semic,lsemic,charge,smass,
     .         basistype)
!
!
C Read information which will used by subroutine ATOM to construct the 
C basis set and the Kleinman-Bylander projectors
C Written by D.Sanchez-Portal. Aug. 1998
C Modified by DSP, July 1999
C  WARNING: the information about the basis set has been split in a
C  somewhat arbitrary fashion: 
C  Arrays lmxo, nzeta and lambda refer to PAO like orbitals and doble-z,
C  triple-z, etc... obtained from the PAO function using any of the 
C  augmentation options.
C  Array polorb contains the number of polarization orbitals desired 
C  for each angular momentum.
C************************INPUT******************************************
C   INTEGER NSP         : Expected total number of different chemical 
C                           species 
C
C************************OUTPUT*****************************************
C   INTEGER IZ(NSP)   : Atomic number of the different species.
C   INTEGER LMXKB(NSP): Maximum angular momentum for the Kleinman-
C                         Bylander projector.
C   INTEGER NKBL(0:LMAXD,NSP): Number of KB projectors for each 
C                         angular momentum
C   REAL*8  EREFKB(NKBMX,0:LMAXD,NSP): Reference energies for the 
C                         calculation of the KB projectors.
C   INTEGER LMXO(NSAMX) : Maximum angular momentum for the PAO orbitals
C                          (WARNING: this does not include polarization
C                            orbitals)
C   INTEGER NZETA(0:LMAXD,MSEMX,NSP): Number of basis functions for each 
C                                  angular momentum (WARNING: this does 
C                                 not include polarization orbitals).
C   REAL*8 RCO(NZETMX,0:LMAXD,NSEMX,NSP):  Cut-off radius for each orbital
C   REAL*8 LAMBDA(NZETMX,0:LMAXD,NSEMX,NSP): Compresion factor fot each
C                                                   orbital
C   CHARACTER*20 ATM_LABEL(NSP):    Species label
C   POLORB(0:LMAXD,NSEMX,NSP):  Number of polarization orbitals 
C                                 perturbatively  calculated from
C                                 the PAO orbitals with angular momentum
C                                 L (This orbitals will have L'=L+1) 
C   LOGICAL SEMIC(NSP)       :  Presence or not of semicore orbitals
C   INTEGER LSEMIC(0:LMAXD,NSP):  Number of semicore shells for each angular
C                                       momentum
C   REAL*8  CHARGE(NSP)  :  Charge state of the ion, only for PAO 
C                             calculation purposes.
C   REAL*8  SMASS(NSP)  :   Atomic mass for each species.
C   CHARACTER*10 BASISTYPE(NSP): Augmentation procedure for each atomic
C                               species
C    Distances in Bohr.
C    Energies in Rydbergs.
C------------------------------------------------------------------C
C  Modules
C
      use precision
      use atmfuncs, only: lmaxd, nsemx, nzetmx, nkbmx
      use fdf

      implicit none

C
C***********************************************************************
C
      integer, intent(in) :: nsp
      integer, intent(out) :: iz(nsp),lmxkb(nsp),lmxo(nsp),
     .     polorb(0:lmaxd,nsemx,nsp),nzeta(0:lmaxd,nsemx,nsp),
     .     lsemic(0:lmaxd,nsp), nkbl(0:lmaxd,nsp)

      double precision, intent(out) ::
     .     rco(nzetmx,0:lmaxd,nsemx,nsp),
     .     lambda(nzetmx,0:lmaxd,nsemx,nsp),
     .     charge(nsp), smass(nsp), erefkb(nkbmx,0:lmaxd,nsp) 
       
      character, intent(out) ::
     .     atm_label(nsp)*20, basistype(nsp)*10 


      logical, intent(out)  ::  semic(nsp)


C      Internal variables

       integer
     .     i, is, l, latm, lmxval,nsm,
     .     config(0:4)

       double precision
     .   ql(0:3), Zval
!
!      Automatic array...
!
       character(len=15), dimension(nsp) :: basis_sizes

       character  basis_size*15,  basistype_generic*10
        
       logical old_block, new_block, chem_block, error
        
       character(len=15), parameter  :: basis_size_default='standard'
       character(len=10), parameter  :: basistype_default='split'
       character(len=1)   :: sym(0:4) = (/ 's','p','d','f','g' /)
C
!------------------------------

C Checking if is there any block containing information
C about the basis set
           old_block=fdf_defined('PAO_basis_and_PS_lmax')
           new_block=fdf_defined('PAO.Basis') 

           if((old_block).and.(new_block)) then 
             write(6,'(2a)')
     .  'read_basis: ERROR both information blocks: ', 
     .      'PAO.Basis and PAO_basis_and_PS_lmax'
             write(6,'(2a)')
     .  'read_basis: ERROR cannot be simultaneously specified ', 
     .      'in the input file'
             call die
           endif 

C  Check that Chemical_species_labels block exists (it must).

          chem_block=fdf_defined('Chemical_species_label')

          if (.not.chem_block) then 
          write(6,100)
          write(6,101) 
          write(6,'(a)')
     .   'read_basis: ERROR: Chemical_species_label block not found'
          write(6,'(a)')
     .   'read_basis: ERROR: The chemical species must be specified.' 
          call die
          endif 

C**********Read species and labels
C
          call rechemsp(atm_label,iz,nsp,nsp)

C**Standard valence charge for each atom and checking for semicore orbitals**
C
        do is=1,nsp
          semic(is)=.false.
          do l=0,lmaxd
            lsemic(l,is)=0
          enddo
        enddo 
        do is=1,nsp  
           do l=0,3
              ql(l)=0.0d0
           enddo   
           Zval=0.0d0 

           if(iz(is).ne.-100) then  
            call lmxofz(abs(iz(is)),lmxval,latm) 
            call qvlofz(abs(iz(is)),ql) 
            do l=0,lmxval
               Zval=Zval+ql(l)
            enddo 
       
            call semicore(Zval,iz(is),atm_label(is),
     .                      semic(is),lsemic(0,is)) 

            if (semic(is)) then 
               do l=0,lmaxd
                 if(lsemic(l,is).ne.0) 
     .               lmxval=max(lmxval,l)
               enddo 
            endif

            if(lmxval.gt.lmaxd) then 
              write(6,'(a)') 
     .  'read_basis: ERROR bad dimensions in atmfuncs'
              call chkdim('read_basis','lmaxd',
     .        lmaxd,lmxval,1) 
            endif 
            if(semic(is).and.old_block) then 
              write(6,'(a)')
     .  'read_basis: ERROR semicore states are incompatible with'
              write(6,'(a)')
     .  'read_basis: ERROR old block PAO_basis_and_PS_lmax'
            endif 

           else

             semic(is)=.false. 

           endif
        
           enddo
C
C************************************************************************



C********  Reading atomic masses **********************************
C********* Default atomic masses given by the function atmass****** 
C
          do is=1,nsp
            if(iz(is).gt.0) then
              smass(is) = atmass(iz(is))
            else
              smass(is) = 1.0d30
            endif  
          enddo  
C
C********Read atomic masses from fdf if a change is wanted***************
C
         call remass(smass,nsp,nsp)
C
C****Augmentation procedure for the construction of the basis set********
C
         basistype_generic=
     .        fdf_string('PAO.BasisType',basistype_default)

         call type_name(basistype_generic) 
         do is=1, nsp
             basistype(is)=basistype_generic
         enddo 

           if(.not.old_block) then 
C******** DESIRED SIZE FOR THE BASIS SET********************************** 
C
           basis_size=fdf_string('PAO.BasisSize',basis_size_default)
           call size_name(basis_size)
           do is=1,nsp
                basis_sizes(is)=basis_size  
           enddo  
           call resizes(atm_label,basis_sizes,nsp,nsp)  
          
C********Maximum angular momentum for the Kleinman-Bylander projectors**** 
C
           do is=1,nsp
               lmxkb(is)=-1
           enddo
           call relmxkb(atm_label,lmxkb,nsp,nsp)  
           endif 
C

C*****READING BASIS SET INFORMATION ***************************************
C*********Initializing lmxo to -1****************************************** 
C********This allows us to check if for a given species some***************
C*********information will be read in the next block***********************
C
           do is=1,nsp           
               lmxo(is)=-1
           enddo    
C*******Initializing arrays nzeta and polorb, which give respectively**** 
Cthe number of PAO and polarization basis functions for a given species**
C***********and angular momentum*****************************************  
           do is=1,nsp
            do nsm=1,nsemx
              do l=0,lmaxd
                  nzeta(l,nsm,is)=0
                  polorb(l,nsm,is)=0
              enddo 
            enddo
           enddo

C****Reading from new block****************************************
C
           if((.not.old_block).and.(new_block)) then 

           call rePAOBasis(atm_label,iz,lmxo,nzeta,rco,lambda,
     .         polorb,charge,basistype,semic,lsemic,nsp)

           endif 

C*******************Reading from old block******************************
C           
          if(old_block) then   
              call reOldBlock(iz,lmxo,lmxkb,nzeta,rco,lambda,
     .        polorb,charge,nsp)
          endif 
C
C*********Reading automatic information about atomic basis sets***
C*********Checking the basis set for each atom******************** 
C         
          error=.false.
          do is=1,nsp
             if((lmxo(is).eq.-1).and.(old_block)) then 
               write(6,'(2a)') 
     .   'read_basis: ERROR not basis set provided for species ',
     .       atm_label(is) 
               write(6,'(a)') 
     .   'read_basis: ERROR in block PAO_basis_and_PS_lmax' 
               call die
             endif 

             if((lmxo(is).eq.-1).and.(iz(is).eq.-100)) then 
                 write(6,'(2a)')
     .      'read_basis: ERROR details for the construction of ',
     .      'the floating Bessel functions' 
               write(6,'(2a)')
     .      'read_basis: ERROR basis set for the species ',
     .       atm_label(is)
               write(6,'(2a)')
     .      'read_basis: ERROR must be specified using ',
     .      'either PAO.Basis or'
               write(6,'(a)')
     .      'read_basis: ERROR PAO_basis_and_PS_lmax blocks' 

               call die
             endif   

             if (lmxo(is).eq.-1) then 
               call autobasis(basis_sizes(is),
     .           abs(iz(is)),
     .           semic(is), lsemic(0,is),
     .           lmxo(is),nzeta(0,1,is),rco(1,0,1,is),
     .           lambda(1,0,1,is),
     .           polorb(0,1,is)) 
                 charge(is)=0.0d0
             endif  

C************Checking the size of the basis set********************
C      
          if(iz(is).ne.-100) then 
             do l=0,3
                ql(l)=0.0d0 
             enddo
             call lmxofz(abs(iz(is)),lmxval,latm)
             call qvlofz(abs(iz(is)),ql)   
             call cnfig(abs(iz(is)),config)

             if (semic(is)) then
               do l=0,lmaxd
                 if(lsemic(l,is).ne.0)
     .               lmxval=max(lmxval,l)
               enddo
             endif
             if (lmxval.gt.lmxo(is)) then 
               write(6,'(/,2a)')
     .     'read_basis: ERROR: For species ', atm_label(is)
               write(6,'(a)')
     .     'read_basis: ERROR: Maximum angular momentum of the basis' 
               write(6,'(a,i4)')
     .     'read_basis: ERROR: orbitals must be at least: ',lmxval 
               error=.true.                
             endif  
            
            do nsm=1,lsemic(0,is)+1           
              
             if(nzeta(0,nsm,is).eq.0) then 
                write(6,'(a)')
     .     'read_basis: ERROR: At least one shell of basis set orbitals' 
                write(6,'(a,i1,3a)')
     .     'read_basis: ERROR: is need for the ',
     .      config(0)-(lsemic(0,is)+1)+nsm,sym(0),
     .     ' in ',atm_label(is)
                error=.true.    
             endif 
            enddo 

             do l=1,lmxval
                do nsm=1,lsemic(l,is)+1
                   if(
     .           (ql(l).gt.1.0d-4).or.
     .             ((nsm.lt.lsemic(l,is)+1).and.(semic(is))) ) then 

                    if(nzeta(l,nsm,is).eq.0) then 
                     write(6,'(2a)')
     .    'read_basis: ERROR: For species ', atm_label(is)
                     write(6,'(a)')
     .    'read_basis: ERROR: There must be at least one shell of ' 
                     write(6,'(a,i1,a)')
     .    'read_basis: ERROR: basis set orbitals for the state',
     .     config(l)-(lsemic(l,is)+1)+nsm,sym(l)  
                      error=.true.    
                  endif
                 endif 
               enddo 
             enddo  
            endif 

         enddo 
C
C         Reading information for multiples KB projectors 
C             
          call reKBblock(atm_label,lmxkb,nkbl,erefkb,
     .               lmxo, lsemic, polorb, nsp, nsp)

          do is=1,nsp              
             if(iz(is).lt.0) then 
                lmxkb(is)=0
                 do l=0,lmaxd
                    nkbl(l,is)=0
                 enddo
             endif
          enddo 
C
          if (error) call die

100   format(/,'redbasis: ',71(1h*))
101   format('redbasis:                  INPUT ERROR')
102   format('rebasis: ',71(1h*))
!
!
!     OLD REDBASIS_SUBS follows
!
!
      CONTAINS

      subroutine size_name(basis_size)

      character basis_size*(*)


      if(basis_size.eq.'STANDARD') basis_size='dzp'
      if(basis_size.eq.'standard') basis_size='dzp'
      if(basis_size.eq.'DZP')  basis_size='dzp'
 
      if(basis_size.eq.'DZ') basis_size='dz'

      if(basis_size.eq.'MINIMAL') basis_size='sz'
      if(basis_size.eq.'minimal')  basis_size='sz'
      if(basis_size.eq.'SZ')  basis_size='sz'
            
      if(basis_size.eq.'SZP') basis_size='szp' 

      if( (basis_size.ne.'szp').and.(basis_size.ne.'sz').and.
     .    (basis_size.ne.'dz') .and.(basis_size.ne.'dzp') ) then 

         write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .   'size_name: Incorrect basis-size option specified,',
     .   ' active options are:',
     .   'SZ','SZP','DZ','and','DZP or STANDARD'

         call die
      endif

      end subroutine size_name


      subroutine resizes(atm_label,basis_sizes,maxs,ns)

c Reading atomic basis sizes for different species.
c
c Reads fdf block. Not necessarily all species have to be given. The 
c ones not given at input will be assumed to have the basis sizes 
c given by the general input PAO.BasisSize, or its default value. 
c      
c ********* INPUT ***************************************************
c integer maxs                   : Maximum number of species
c integer ns                     : Number of species
c character*20 atm_label         : Labels of the different species
c character*15 basis_sizes(maxs) : Basis sizes as given by the input 
c                                  PAO.BasisSize, or its default value 
c ********* OUTPUT **************************************************
c character*15 basis_sizes(maxs) : Basis sizes specified for each species
c                                  in the block PAO.BasisSizes
c *******************************************************************

      integer           maxs, ns
      character*15      basis_sizes(maxs) 
      character*20      atm_label(maxs)
      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80, label_read*20
      integer           ni, nn, nr, nv, ist, nst, iu, ns_read, ilabel
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)


c check for block and read

      if ( fdf_block('PAO.BasisSizes',iu) ) then

         nst = 0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line 
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (nn .ge. 2)  then 
              label_read=names(lc(0)+1:lc(1))
              if(label_read.ne.'%endblock') then 

c determine the species index from the label 
              ns_read=0 
              do ilabel=1,ns
                 if(atm_label(ilabel).eq.label_read) then 
                  ns_read=ilabel
                 endif  
              enddo 
              if(ns_read.eq.0) then 
                 write(6,"(/,2a)") 
     .      'resizes: WARNING: Species not defined. Ignored label= ',
     .        label_read
              else
                 basis_sizes(ns_read) = names(lc(1)+1:lc(2))
                 call size_name(basis_sizes( ns_read ) )
                 write(6,'(4a)')
     .            'resizes: Read basis size for species ',
     .             names(lc(0)+1:lc(1)),' = ',basis_sizes( ns_read ) 

                   nst = nst + 1
                   if (nst .gt. maxs) then
                    call die('resizes: BAD DIMENSIONS. maxs too small')
                   endif
 
              endif  


               
              else
                 return
              endif 
          
 
            else
              return
            endif
         enddo

         call die('resizes: Too many entries in PAO.BasisSizes')

      endif

 50   continue

      end subroutine resizes

      subroutine rechemsp(atm_label,iz,maxs,ns)

c *******************************************************************
c Reading the labels and atomic numbers for different species. 
c The labels are neccesary to 
c identify the files which contain the information about the 
c corresponding pseudopotentials and also for the
c names of some output files. 
c
c Reads fdf block. 
c
c  All the species must be specified!!!!
c      
c Written by D. Sanchez-Portal. Aug. 1998. 
c ********* INPUT ***************************************************
c integer maxs                   : Maximum number of species
c integer ns                     : Number of species
c ********* OUTPUT **************************************************
c character*20 atm_label(maxs)   : Labels for the species
c integer      iz(maxs)          : Atomic number of the 
c                                    species
c *******************************************************************

C
      implicit          none
      integer           maxs, ns, iz(maxs)
      character*20      atm_label(maxs) 
      external          parse

c ----------------------------------------------------------------------------

      character         line*130, names*80
      integer           ni, nn, nr, nv, ist, nst, iu, nsread
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)


c check for block and read

      if ( fdf_block('Chemical_species_label',iu) ) then

         nst = 0
         nsread=0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (ni .ge. 2)  then
              nst=nst+1
              if (nst .gt. maxs) then
                  call die('rechemsp: BAD DIMENSIONS. maxs too small')
              endif
              if ((integs(1).gt.ns) .or. (integs(1).lt.1) ) then
                write(6,"(/,a,i4)") 
     .   'rechemsp: WARNING: Species out of range. Ignored is =',
     .           integs(1)
              else
                 nsread=nsread+1
                 atm_label(integs(1)) = names(lc(0)+1:lc(1))
                 iz(integs(1))=integs(2)
                 write(6,"(a, i4, a, i4, 2a)")
     .   'rechemsp: Read atomic number and label for species', 
     .    integs(1), ' as ' , integs(2), '  ',names(lc(0)+1:lc(1))
              endif 


            else
              if (nsread.lt.ns) then 
                write(6,"(/,a,i4,a)")
     .        'rechemsp: ERROR: Atomic number and label read only for ',
     .         nsread,' species'
                write(6,"(a)") 
     .        'rechemsp: ERROR: There are some undefined species'
                call die
              endif  

              return
            endif 


         enddo

        call die(
     $       'rechemsp: Too many entries in Chemical_species_label')

      endif

 50   continue

      end subroutine rechemsp

      subroutine reKBblock(atm_label,lmxkb,nkbl,erefkb,
     .               lmxo,lsemic,polorb,maxs,ns)
c *******************************************************************
c Reading the number of KB projectors for each angular momentum
c and the reference energies for their construction.
c
c Reads fdf block. 
c The information should be compatible with that provided in 
c the block 'PS.lmax'.
c Not necessarily all species have to be given.
c The values not specify in this block are taken as the
c default value:
c *Maximum angular momentum of the KB projectors,
c lmxkb=lmxo+1, where lmxo is the maximum angular momentum
c of the orbitals in the basis set.
c *Only one KB projector for each angular momentum,
c nkbl(l,is)=1
c *Reference energy:
c erefkb(1,l,is)>=1.0d3 Ry, with this default value or bigger, 
C the program 
c will use the eigenstates as the reference states to build
c the projectors
c with erefkb(1,l,is)=< -1.0d3 Ry, with this default value
c or lower, the program will use the energy derivative of 
c the previous state to build the projector
c Written by D. Sanchez-Portal, July 1999.
c ********* INPUT ***************************************************
c integer maxs                   : Maximum number of species
c integer ns                     : Number of species
c character*20  atm_label(maxs)  : Label for the different species
c integer lmxo(maxs)             : Maximum angular momentum for the 
c                                  basis orbitals
c lsemic(0:lmaxd,ns)          : Number of semicore+valence shells
c                                  for each angular momentum 
c integer polorb(0:lmaxd,nsemx,ns): Numebr of polarisation orbitals
c                                  perturbatively calculated from
c                                  the PAO orbitals with angular
c                                  momentum L
c ********* OUTPUT **************************************************
c integer lmxkb(maxs)          : Maximum angular momentum for KB projectors
c                                as given in the PS.lmax block
c integer nkbl(0:lmaxd,maxs)   : Number of KB projector for each l.
c real*8  erefkb(nkbmx,0:lmaxd,maxs) : Reference energies to build 
c                                      the KB projectors.
c *******************************************************************

      implicit          none
      integer           maxs, ns
      integer           lmxkb(maxs),nkbl(0:lmaxd,maxs),lmxo(ns)
      integer           polorb(0:lmaxd,nsemx,ns)
      real*8            erefkb(nkbmx,0:lmaxd,maxs)
      character*20      atm_label(maxs)
      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
!
!     Automatic array
!
      integer           ns_control(ns)
!
      character         line*130, names*80, label_read*20
      character*10      defunit, unitstr
      integer           ni, nn, nr, nv, ist, iu, ns_read, ilabel
      integer           integs(4), lastc, lc(0:3), nzt, izt, l
      integer           nsh,l_control(0:lmaxd), l_read, ish, ish2
      integer           nzt_overf, l_overf, lmax, l_contr, lpol
      integer           lsemic(0:lmaxd,ns)  
      double precision  reals(4), values(4)
      logical           overf_zet, overf_l
!      external          fdf_convfac 
!      double precision  fdf_convfac

c Default unit for the energy

         parameter(defunit='Ry')

c check for block and read

      do ist=1,ns
         ns_control(ist)=0
      enddo
      do ist=1,ns
        do l=0,lmaxd
          nkbl(l,ist)=0
        enddo
      enddo

c initialise nzt to keep compiler happy
      nzt = 0
      lpol = 0

      if ( fdf_block('PS.KBprojectors',iu) ) then

         overf_zet=.false.
         overf_l=.false.
         nzt_overf=0
         l_overf=0
 
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )


C**Reading reference energies (if presents) for the last***********
C**projectors of the previous chemical specie**************************
             if((ist.ne.1).and.(nr.ge.nzt).and.(ni.eq.0)) then
               if(nn.ge.1) then 
                  unitstr=names(lc(0)+1:lc(1))
                else
                  unitstr=defunit
                endif
                do izt=1,nzt
                  erefKB(izt,l_read,ns_read)=
     .            reals(izt)*fdf_convfac(unitstr,defunit)
                enddo
                read(iu,'(a)') line
                lastc = index(line,'#') - 1
                if (lastc .le. 0) lastc = len(line)
                call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
              endif
C*********************************************************************


            if ((nn.ge.1).and.(ni.ge.1)) then

              label_read=names(lc(0)+1:lc(1))
              if(label_read.ne.'%endblock') then

c determine the species index from the label
              ns_read=0
              do ilabel=1,ns
                 if(atm_label(ilabel).eq.label_read) then
                  ns_read=ilabel
                 endif
              enddo
              if(ns_read.eq.0) then

                 write(6,"(/,2a)")
     .      'reKBblock: WARNING: Species not defined. Ignored label= ',
     .        label_read

               else
c reading the number of shells
                nsh=integs(1)
                
                lmax=0
                do ish=1,nsh
                  read(iu,'(a)') line
                  lastc = index(line,'#') - 1
                  if (lastc .le. 0) lastc = len(line)
                  call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
                  if((ish.ne.1).and.(nr.ge.nzt)) then
                      if(nn.ge.1) then 
                         unitstr=names(lc(0)+1:lc(1))
                      else
                         unitstr=defunit
                      endif
                      do izt=1,nzt
                        erefKB(izt,l_read,ns_read)=
     .                         reals(izt)*fdf_convfac(unitstr,defunit)
                      enddo
                      read(iu,'(a)') line
                      lastc = index(line,'#') - 1
                      if (lastc .le. 0) lastc = len(line)
                call parse( line(1:lastc), nn, lc, names, nv, values,
     .                   ni, integs, nr, reals )
                  endif
            
                  if(ni.lt.2) goto 40
                   l_read=integs(1)
                   do ish2=1,ish-1
                     if(l_control(ish2).eq.l_read) then
                      write(6,'(a)')
     .          'reKBblock: ERROR duplicate information for species ',
     .                atm_label(ns_read)
                      call die
                     endif
                   enddo
                   l_control(ish)=l_read

                   nzt=integs(2)
                   if(l_read.gt.lmax)  lmax=l_read
                   if(nzt.gt.nzetmx) overf_zet=.true.
                   if(nzt.gt.nzt_overf) nzt_overf=nzt
                   if(l_read.gt.lmaxd) overf_l=.true.
                   if(l_read.gt.l_overf) l_overf=l_read
                   if(.not.overf_l) nkbl(l_read,ns_read)=nzt
             
                   if((.not.overf_zet).and.(.not.overf_l)) then
                      do izt=1,nzt
                        erefKB(izt,l_read,ns_read)=huge(1.d0)
                      enddo
                   endif

30              continue
                enddo 
                if((lmxkb(ns_read).ne.-1).and.
     .             (lmxkb(ns_read).lt.lmax) ) then  
                     write(6,'(/,a)')
     .  'reKBblock: ERROR: Data read in PS.KBprojectors block'
                     write(6,'(a)')
     .  'reKBblock: ERROR are incompatible with data from PS.lmax block'
                     call die
                endif
                lmxkb(ns_read)=max(lmax,lmxkb(ns_read))
                ns_control(ns_read)=1
                do l=0,lmxkb(ns_read)
                    l_contr=0
                    do ish=1,nsh
                        if(l.eq.l_control(ish)) l_contr=1
                    enddo
                    if(l_contr.eq.0) then
                        nkbl(l,ns_read)=1
                        erefKB(1,l,ns_read)=huge(1.d0)
                    endif
                enddo

             endif




             else   
               goto 50
             endif

           else 
             goto 50
           endif

         enddo

         if((overf_zet).or.(overf_l)) then
            write(6,'(a)')
     .      'reKBblock: ERROR bad dimensions in atmfuncs'
            if(overf_zet) call chkdim('reKBblock','nkbmx',
     .      nzetmx,nzt_overf,1)
            if(overf_l) call chkdim('reKBblock','lmaxd',
     .      lmaxd,l_overf,1)

         endif

         call die('reKBblock: Too many entries in PS.KBprojectors')

 40   continue 

        if((overf_zet).or.(overf_l)) then
           write(6,'(a)')
     .     'reKBblock: ERROR bad dimensions in atmfuncs'
           if(overf_zet) call chkdim('reKBblock','nkbmx',
     .     nzetmx,nzt_overf,1)
          if(overf_l) call chkdim('reKBblock','lmaxd',
     .    lmaxd,l_overf,1)
        endif


       write(6,'(a)')
     .    'reKBblock: ERROR reading PS.KBprojectors block'
       write(6,'(a)')
     .    'reKBblock: ERROR check the block sintaxis'
       call die
       
 50   continue


       if((overf_zet).or.(overf_l)) then
               write(6,'(a)')
     .      'reKBblock: ERROR bad dimensions in atmfuncs'
         if(overf_zet) call chkdim('reKBblock','nkbmx',
     .          nzetmx,nzt_overf,1)
         if(overf_l) call chkdim('reKBblock','lmaxd',
     .          lmaxd,l_overf,1)

         call die
        endif


        endif

        do ist=1,ns
           if(ns_control(ist).eq.0) then 
             if(lmxkb(ist).eq.-1) lmxkb(ist)=lmxo(ist)+1
             do l=0,lmxkb(ist)
               nkbl(l,ist)=lsemic(l,ist)+1
               do izt=1, lsemic(l,ist)+1
                      erefKB(izt,l,ist)=huge(1.d0)
               enddo 
             enddo 
             do l=0,lmxo(ist)
               do izt=1, lsemic(l,ist)+1
                   if(polorb(l,izt,ist).gt.0) lpol=l+1
               enddo
             enddo
             call chkdim('reKBblock','lmaxd',
     .          lmaxd,lpol+1,1)
             if(lpol+1.gt.lmxkb(ist)) then
                  lmxkb(ist)=lpol+1
                  nkbl(lpol+1,ist)=1
                  erefKB(1,lpol+1,ist)=huge(1.d0)
             endif
           endif
        enddo
  
      end subroutine reKBblock


      subroutine relmxkb(atm_label,lmxkb,maxs,ns)

c *******************************************************************
c Reading the maximum angular momentum of the Kleinman-Bylander 
c projectors for the different species.
c
c Reads fdf block. Not necessarily all species have to be given. 
c The values not specify in this block are taken as the  
c default value, lmxo+1, where lmxo is the maximum angular momentum 
c of the orbitals in the basis set.
c      
c Written by D. Sanchez-Portal. Aug. 1998. 
c ********* INPUT ***************************************************
c integer maxs                   : Maximum number of species
c integer ns                     : Number of species
c character*20  atm_label(maxs)  : Label for the different species
c ********* OUTPUT **************************************************
c integer lmxkb(maxs)          : Maximum angular momentum for KB projectors
c                                as given in the PS.lmax block   
c *******************************************************************

      implicit          none
      integer           maxs, ns
      integer           lmxkb(maxs)
      character*20      atm_label(maxs)
      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80, label_read*20
      integer           ni, nn, nr, nv, ist, nst, iu, ns_read, ilabel
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)


c check for block and read

      if ( fdf_block('PS.lmax',iu) ) then

         nst = 0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if ((ni .ge. 1).and.(nn.ge.1))  then
              label_read=names(lc(0)+1:lc(1))


c determine the species index from the label
              ns_read=0
              do ilabel=1,ns
                 if(atm_label(ilabel).eq.label_read) then
                  ns_read=ilabel
                 endif
              enddo
              if(ns_read.eq.0) then
                 write(6,"(/,2a)")
     .      'relmxkb: WARNING: Species not defined. Ignored label= ',
     .        label_read
              else
                 lmxkb(ns_read)=integs(1)
                 write(6,"(a, i4, 2a)")
     .            'relmxkb: Read lmxkb= ',integs(1), 
     .            ' for species ', label_read 
                   nst = nst + 1
                   if (nst .gt. maxs) then
                     call die('relmxkb: BAD DIMENSIONS. maxs too small')
                   endif

              endif 

            else
              return
            endif
         enddo

         call die('relmxkb: Too many entries in AtomicMass')
      endif

 50   continue

      end subroutine relmxkb

        subroutine semicore(Zval,Z,atm_label,semic,lsemic)
c *******************************************************************
c Reading information for basis generation.
c
c Reading the valence charge provided by the pseudopotential file
c it can be deduced if any semicore internal states have to be included 
c into the valence shell.
c
c Written by D. Sanchez-Portal. Aug. 1998.
c Modified by DSP, July 1999
c ********* INPUT ***************************************************
c double precision Zval          : Expected valence charge of the atom
c character*20 atm_label         : Labels for the different species
c ********* OUTPUT **************************************************
c logical      semic             : if .true. there are semicore states 
c                                  if .false. there are not.
c integer      lsemic            : Angular momentum of the semicore shell
c *******************************************************************

       implicit none      
       double precision Zval, b, a, Zval_read, tiny, charge_loc
       character atm_label*20, paste*50, fname*50, namatm*2, 
     .   icorr*2, irel*3, nicore*4,
     .   method(6)*10, text*70
       integer lsemic(0:lmaxd), nr, npotd, npotu, i,ndiff, Z
       logical semic,found
       integer lun
       
       parameter(tiny=1.0d-5) 

 
C****Open the file with the pseudopotential information**************

           fname = paste(atm_label,'.vps')
           inquire(file=fname, exist=found)
           if (.not.found) then
             write(6,'(/,2a,a20)') 'semicore: WARNING: ',
     .         'Pseudopotential file not found: ', fname
             fname = paste(atm_label,'.psatom.data')
             write(6,'(2a)') 'ATOM: WARNING: Looking for ', fname
             inquire(file=fname, exist=found)
           endif
           if (.not.found) then
             write(6,'(/,2a,a20,/)') 'semicore: ERROR: ',
     .         'Pseudopotential file not found: ', fname
             call die
           endif

           call io_assign(lun)
           open(lun,file=fname,form='unformatted',status='unknown')

           read(lun) namatm, icorr, irel, nicore,
     .     (method(i),i=1,6), text,
     .     npotd, npotu, nr, b, a, Zval_read
           call io_close(lun)
C***************************************************************************
             if(abs(Zval-Zval_read).lt.tiny) then 
                   semic=.false.
                   return 
             else 

             ndiff=nint(abs(Zval-Zval_read))
             if(abs(ndiff-abs(Zval-Zval_read)).gt.tiny) then 
                write(6,'(2a)') 
     .   'semicore: ERROR expected valence charge for species ', 
     .    atm_label
                write(6,'(2a)')
     .  'semicore: ERROR and the value read from the file ', fname
                write(6,'(a,f6.3,a,f6.3)')
     .  'semicore: ERROR differ:  Zval(expected)= ', Zval,
     .  ' Zval(read)= ',Zval_read 
                call die
             endif 

             semic=.true.
             charge_loc=Zval_read-Zval
             write(6,'(3a)')
     .       'semicore: ',
     .       ' semicore shell included in the valence for species ',
     .        atm_label

             call find_semi(abs(z),charge_loc,lsemic) 
             do i=0,lmaxd
              if(lsemic(i)+1.gt.nsemx) then 
                write(6,'(a)')
     .                   'semicore: ERROR parameter nsemx in atmfuncs'
                write(6,'(a,i2)')
     .                'semicore: ERROR must be increased to at least',
     .                lsemic(i)+1
              endif
             enddo

 
             endif 
     

         end subroutine semicore
!
        subroutine autobasis(basis_size,iz,
     .       semic, lsemic,
     .       lmxo,nzeta,rco,lambda,
     .       polorb)
         
C Generates automatic basis set from the information given 
C in the block PAO.BasisSizes
C Written by D. Sanchez-Portal, Aug. 1998.
C Modified by DSP, July 1999

         implicit none
        
         character basis_size*15
         integer lmxo, nzeta(0:lmaxd,nsemx), polorb(0:lmaxd,nsemx),
     .     lsemic(0:lmaxd), iz,nsm
         logical semic 
         double precision rco(nzetmx,0:lmaxd,nsemx), 
     .       lambda(nzetmx,0:lmaxd,nsemx)
       
C*******************internal variables*********************************
         
         double precision ql(0:3)
         integer lmxval, latm, l, nzt, izt, lmax, nval, npol
         logical pol_contr, overflow
         
         overflow=.false.
         call lmxofz(iz,lmxval,latm)

C Initialise ql - JDG Nov 99
         ql = 0.0d0

         call qvlofz(iz,ql) 
         
         lmax=lmxval
         if(semic) then
           do l=0,lmaxd
               if(lsemic(l).gt.0) lmax=max(lmax,l)
           enddo 
         endif
         lmxval=max(lmxval,lmax)
         if(lmax.gt.lmaxd) overflow=.true.
         lmxo=lmxval
         call size_name(basis_size)  
         if((basis_size.eq.'dzp').or.(basis_size.eq.'dz') ) then  
               izt=nzetmx
               if(izt.lt.2)  write(6,'(/2a)') 
     .        'autobasis: ERROR: Parameter nzetmx in atmfuncs',
     .        'must be increased to 2'       
         endif

         pol_contr=.false. 
         do l=0,lmxval
             
             nval=lsemic(l)+1
             if(l.gt.0) npol=lsemic(l-1)+1

             if(abs(ql(l)).lt.1.0d-4) then 
                
                  if(l.eq.0) then 
                     if ((basis_size.eq.'sz').or.
     .                    (basis_size.eq.'szp')) nzt=1

                     if ((basis_size.eq.'dz').or.
     .                    (basis_size.eq.'dzp')) nzt=2   

                     if(.not. overflow) then 
                       nzeta(l,nval)=nzt
                       do izt=1,nzt
                           rco(izt,l,nval)=0.0d0
                           lambda(izt,l,nval)=0.0d0
                        enddo
                     endif 

                  elseif(  ((basis_size.eq.'szp').or.
     .                      (basis_size.eq.'dzp'))
     .                        .and. (.not.pol_contr)
     .                        .and. (l.ne.0)) then
                          pol_contr=.true.   
                          if(.not. overflow) then 
                             polorb(l-1,npol)=1
                             nzeta(l,nval)=0  
                          endif 
                  elseif( (pol_contr)
     .                        .and. (l.ne.0)) then 
                          if(.not. overflow) then
                             nzeta(l,nval)=0  
                          endif 
                  endif


 
             endif 
            
             if(ql(l).gt.1.0d-4) then   


                if ((basis_size.eq.'sz').or.
     .                    (basis_size.eq.'szp')) nzt=1 
    
                if ((basis_size.eq.'dz').or.
     .                    (basis_size.eq.'dzp')) nzt=2
                if(.not. overflow) then 
                   nzeta(l,nval)=nzt
                   do izt=1,nzt
                       rco(izt,l,nval)=0.0d0 
                       lambda(izt,l,nval)=0.0d0 
                   enddo 
                endif 

             endif 
 

             if((semic).and.(lsemic(l).gt.0)) then
                if(.not. overflow) then
                   do nsm=1,lsemic(l)
                      nzeta(l,nsm)=1
                      rco(1,l,nsm)=0.0d0
                      lambda(1,l,nsm)=0.0d0
                   enddo
                endif
             endif

                         
           enddo      

           if((basis_size.eq.'szp').or.(basis_size.eq.'dzp') ) then
            if(.not.pol_contr) then 
              lmxval=lmxval+1  
c             lmax=max(lmax, lmxval) 
c             if(lmax.gt.lmaxd) overflow=.true. 
              if(.not.overflow) then 
                 npol=lsemic(lmxval-1)+1
                 polorb(lmxval-1,npol)=1 
              endif 
            endif 
           endif

         
           if(overflow) then 
               write(6,'(/a)')
     .          'autobasis: ERROR: Parameter lmaxd in atmfuncs' 
               write(6,'(a,i4)')
     .          'autobasis: ERROR: must be increased to at least ',lmax 
               call die
           endif
    
         end  subroutine autobasis

        subroutine type_name(basistype)

C Written by D. Sanchez-Portal, Aug. 1998.

         character basistype*(*) 

         if(basistype.eq.'NODES') then
               basistype='nodes'
         elseif(basistype.eq.'nodes') then
         elseif(basistype.eq.'NONODES') then
              basistype='nonodes'
         elseif(basistype.eq.'nonodes') then
         elseif(basistype.eq.'SPLIT') then
              basistype='split'
         elseif(basistype.eq.'split') then
         elseif(basistype.eq.'SPLITGAUSS') then
              basistype='splitgauss'
         elseif(basistype.eq.'splitgauss') then
         elseif(basistype.eq.'USER') then
              basistype='user'
         elseif(basistype.eq.'user') then
         else

              write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .        'type_name: Incorrect basis-type option specified,',
     .        ' active options are:',
     .        'NODES','SPLIT','USER','SPLITGAUSS','NONODES'
              call die
         endif

        end subroutine type_name

         subroutine rePAOBasis(atm_label,iz,lmxo,nzeta,rco,lambda,
     .         polorb,charge,basistype,semic,lsemic,ns)

c *******************************************************************
c Reading information for basis generation of different species from
c block PAO.Basis
c
c Reads fdf block. Not necessarily all species have to be given. The 
c basis set for those species not given in this block will be generated
c automatically.
c      
c Written by D. Sanchez-Portal. Aug. 1998. 
c *******************************************************************

      implicit  none 
 
      integer   ns, lmxo(ns), polorb(0:lmaxd,nsemx,ns)
      integer   nzeta(0:lmaxd,nsemx,ns), iz(ns)
      integer   lsemic(0:lmaxd,ns)
      double precision
     .   rco(nzetmx,0:lmaxd,nsemx,ns), 
     .   lambda(nzetmx,0:lmaxd,nsemx,ns),
     .   charge(ns)
      character atm_label(ns)*20, basistype(ns)*10
      logical    semic(ns)
      external   parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80, label_read*20
      integer           ni, nn, nr, nv, ist, iu,  ns_read, ilabel
      integer           integs(4), lastc, lc(0:3), nzt, izt, l, ipol 
      integer           nsh,l_control((lmaxd+1)*nsemx)
      integer           l_read, ish, ish2, nsm_control((lmaxd+1)*nsemx)
      integer           nzt_overf, l_overf, lmax, l_contr, nsm_read
      integer           nsm_overf, nsm_max, nsm, nsm_contr
      integer           config(0:lmaxd)
      double precision  reals(4), values(4)
      logical           overf_zet, overf_l, overf_nsm

c initialise nzt to keep compiler happy
      nzt = 0

c check for block and read

      if ( fdf_block('PAO.Basis',iu) ) then
         
         overf_zet=.false.
         overf_l=.false.   
         overf_nsm=.false. 
         nzt_overf=0 
         l_overf=0
         nsm_overf=0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals ) 


C**Reading compression factors (if presents) for the last***********
C**orbital of the previous chemical specie************************** 
             if((ist.ne.1).and.(nr.ge.nzt).and.(ni.eq.0)) then 
                do izt=1,nzt 
                    lambda(izt,l_read,nsm_read,ns_read)=reals(izt)
                enddo
                 read(iu,'(a)') line
                 lastc = index(line,'#') - 1
                 if (lastc .le. 0) lastc = len(line)
                 call parse( line(1:lastc), nn, lc, names, nv, values,
     .                   ni, integs, nr, reals )
              endif
C*********************************************************************


            if ((nn.ge.1).and.(ni.ge.1)) then 

              label_read=names(lc(0)+1:lc(1))
              if(label_read.ne.'%endblock') then

c determine the species index from the label
              ns_read=0
              do ilabel=1,ns
                 if(atm_label(ilabel).eq.label_read) then
                  ns_read=ilabel
                 endif
              enddo
              if(ns_read.eq.0) then 

                 write(6,"(/,2a)")
     .      'rePAOBasis: WARNING: Species not defined. Ignored label= ',
     .        label_read    

               else
c reading the augmentation procedure for the basis set
                if(nn.ge.2) then 
                   basistype(ns_read)=names(lc(1)+1:lc(2))
                   call type_name(basistype(ns_read))
                endif 
c reading the number of shells                   
                nsh=integs(1) 
c reading charge state of the atom for the basis generation
                if (nr.gt.0) then 
                    charge(ns_read)=reals(1)   
                else
                    charge(ns_read)=0.0d0 
                endif
C calculating the usual valence configuration to determine
C the semicore-shell number from the the read principal 
C quantum number 
                do l=0,lmaxd
                    config(l)=l+1
                enddo                     
                call cnfig(abs(iz(ns_read)),config)
                lmax=0
                nsm_max=0
                do ish=1,nsh
                  read(iu,'(a)') line
                  lastc = index(line,'#') - 1
                  if (lastc .le. 0) lastc = len(line)
                  call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
                  if((ish.ne.1).and.(nr.ge.nzt)) then 
                      do izt=1,nzt 
                        lambda(izt,l_read,nsm_read,ns_read)=reals(izt)
                      enddo
                      read(iu,'(a)') line
                      lastc = index(line,'#') - 1
                      if (lastc .le. 0) lastc = len(line)
                call parse( line(1:lastc), nn, lc, names, nv, values,
     .                   ni, integs, nr, reals )  
                  endif 

                  if(ni.lt.2) goto 40

                  l_read=integs(1)

                  if((nn.gt.0).and.(iz(ns_read).ne.-100)) then 
                   if(names(lc(0)+1:lc(1)).ne.'P') then  
                     if(lc(1)-lc(0).ne.3) then 
                        write(6,'(a)')
     .          'rePAOBasis: ERROR sintaxis error in the specification'
                        write(6,'(a)')
     .          'rePAOBasis: ERROR of the principal quantum number'

                        call die
                     endif
                     nsm_read=qnumber(names(lc(0)+1:lc(1)))
                     nsm_read=lsemic(l_read,ns_read)+1
     .                       -(config(l_read)-nsm_read)
                   else
                     if(semic(ns_read)) then
                       write(6,'(a)')
     .          'rePAOBasis: ERROR with species with semicore states'
                       write(6,'(2a)')
     .          'rePAOBasis: ERROR you must specified the main',
     .           ' quantum number of the shell'
                       call die
                     endif
                     nsm_read=1
                   endif
                  else
                     if(semic(ns_read)) then
                       write(6,'(a)')
     .          'rePAOBasis: ERROR with species with semicore states'
                       write(6,'(2a)')
     .          'rePAOBasis: ERROR you must specified the main',
     .           ' quantum number of the shell'  
                       call die
                     endif
                     nsm_read=1
                  endif
                  if(iz(ns_read).eq.-100) nsm_read=1

                   do ish2=1,ish-1
                     if((l_control(ish2).eq.l_read).and.
     .                   (nsm_control(ish2).eq.nsm_read)) then  
                      write(6,'(a)')
     .          'rePAOBasis: ERROR duplicate information for species ',
     .                atm_label(ns_read)    

                      call die

                     endif 
                   enddo 
                   l_control(ish)=l_read
                   nsm_control(ish)=nsm_read
                   nzt=integs(2)
                   if(l_read.gt.lmax)  lmax=l_read 
                   if(nsm_read.gt.nsm_max) nsm_max=nsm_read
                   if(nsm_read.gt.nsemx) overf_nsm=.true.
                   if(nsm_read.gt.nsm_overf) nsm_overf=nsm_read 
                   if(nzt.gt.nzetmx) overf_zet=.true.
                   if(nzt.gt.nzt_overf) nzt_overf=nzt
                   if(l_read.gt.lmaxd) overf_l=.true. 
                   if(l_read.gt.l_overf) l_overf=l_read
                   if((.not.overf_l).and.
     .             (.not.overf_nsm)) nzeta(l_read,nsm_read,ns_read)=nzt           


                   if(((nn.ge.1).and.(names(1:1).eq.'P')).or.
     .          ((nn.ge.2).and.(names(lc(1)+1:lc(1)+1).eq.'P'))) then 
 
C****No polarization option if floating Bessel functions (Z=-100)
                    if(iz(ns_read).ne.-100) then

c                     if(l_read+1.gt.lmaxd) overf_l=.true. 
c                     if(l_read+1.gt.l_overf) l_overf=l_read+1
                      
                      if(ni.ge.3) then 
                        ipol=integs(3) 
                      else
                        ipol=1
                      endif   
                      if((.not.overf_l).and.
     .                  (.not.overf_nsm)) 
     .                    polorb(l_read,nsm_read,ns_read)=ipol
                      if(ipol.gt.nzetmx) 
     .                  overf_zet=.true.
                      if(ipol.gt.nzt_overf) 
     .                  nzt_overf=ipol
                   else
                    write(6,'(a)')
     .        'rePAOBasis: WARNING polarization option not active with ' 
                    write(6,'(a)')
     .        'rePAOBasis: WARNING  Z=-100 (Floating Bessel functions)'  
                   endif 

                  endif
 
               if((.not.overf_zet).and.(.not.overf_l)
     .                              .and.(.not.overf_nsm)) then  
                   read(iu,'(a)') line 
                   lastc = index(line,'#') - 1
                   if (lastc .le. 0) lastc = len(line) 
           call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )

                   if(nr.lt.nzt) goto 40 
                   do izt=1,nzt 
                      rco(izt,l_read,nsm_read,ns_read)=reals(izt) 
                      lambda(izt,l_read,nsm_read,ns_read)=0.0d0 
                   enddo 
                endif

30              continue        
                enddo  
                lmxo(ns_read)=lmax 

                do l=0,lmax
                  do nsm=1,nsm_max
                    l_contr=0
                    nsm_contr=0
                    do ish=1,nsh
                        if(l.eq.l_control(ish)) l_contr=1
                        if(nsm.eq.nsm_control(ish)) nsm_contr=1
                    enddo 
                    if((l_contr.eq.0).and.(nsm_contr.eq.0)) 
     .                   nzeta(l,nsm_read,ns_read)=0
                  enddo
                enddo

             endif 

             


             else    
               goto 50 
             endif 
           
           else  
             goto 50
           endif 

         enddo 

         if((overf_zet).or.(overf_l).or.(overf_nsm)) then    
            write(6,'(a)')
     .      'rePAOBasis: ERROR bad dimensions in atmfuncs'
            if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .      nzetmx,nzt_overf,1)
            if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .      lmaxd,l_overf,1)
            if(overf_nsm)  call chkdim('rePAOBasis','nsmx',
     .      nsemx-1,nsm_overf-1,1)
         endif 
         
         call die('rePAOBasis: Too many entries in PAO.Basis')
      endif     
 40   continue  


        if((overf_zet).or.(overf_l).or.(overf_nsm)) then  
           write(6,'(a)') 
     .     'rePAOBasis: ERROR bad dimensions in atmfuncs'
           if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .     nzetmx,nzt_overf,1)
           if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .     lmaxd,l_overf,1)
            if(overf_nsm)  call chkdim('rePAOBasis','nsmx',
     .      nsemx-1,nsm_overf-1,1)
        endif 


       write(6,'(a)')
     .    'rePAOBasis: ERROR reading PAO.Basis block'
       write(6,'(a)')
     .    'rePAOBasis: ERROR check the block sintaxis'
       call die

 50   continue 


       if((overf_zet).or.(overf_l).or.(overf_nsm)) then   
               write(6,'(a)') 
     .      'rePAOBasis: ERROR bad dimensions in atmfuncs'
         if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .          nzetmx,nzt_overf,1)
         if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .          lmaxd,l_overf,1)
         if(overf_nsm)  call chkdim('rePAOBasis','nsmx',
     .      nsemx-1,nsm_overf-1,1)

         call die
        endif 

      end subroutine rePAObasis

         subroutine reOldBlock(iz,lmxo,lmxkb,nzeta,rco,lambda,
     .        polorb,charge,ns)

c *******************************************************************
c Reading information for basis generation of different species
c from obsolete block PAO_basis_and_PS_lmax.
C This is a highly unrecommended input format due to its lack of 
C flexibility. This routine has only 
c been included for compatibility with old input files.
c      
c Written by D. Sanchez-Portal. Aug. 1998. 
C Modified by DSP, July. 1999
c *******************************************************************

      implicit  none 
 
      integer   ns, lmxo(ns), polorb(0:lmaxd,nsemx,ns)
      integer   nzeta(0:lmaxd,nsemx,ns), iz(ns), lmxkb(ns)
      double precision
     .   rco(nzetmx,0:lmaxd,nsemx,ns), 
     .   lambda(nzetmx,0:lmaxd,nsemx,ns),
     .   charge(ns)
      external   parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      integer           ist, iu
      integer           nzt, izt, l
      integer           nzt_overf, l_overf, izsr 
      integer           numb_pol, numb_pol_default, il, is

      logical           overf_zet, overf_l, polarization, pol_default
    
      parameter (pol_default=.false.)
      parameter (numb_pol_default=1) 
       
c check for block and read


       overf_zet=.false.
       overf_l=.false.   
       nzt_overf=0
       l_overf=0  


       polarization=fdf_boolean('PAO.PolarizationOrbitals',
     .       pol_default)
         if(polarization) then
            numb_pol=fdf_integer('PAO.SplitPolarizationOrbitals',
     .         numb_pol_default) 
            nzt_overf=numb_pol
            if(numb_pol.gt.nzetmx) overf_zet=.true.
          endif
          
         if ( fdf_block('PAO_basis_and_PS_lmax',iu) ) then


         do ist = 1, ns  
            read(iu,*) is, izsr, lmxo(is), lmxkb(is) 
            charge(is)=0.0d0
            if(iz(is).ne.izsr) then
               write(6,'(2a)') 'reOldBlock: ERROR:',
     .  ' Atomic numbers specified in Chemical_species_label block'
               write(6,'(3a)') 'reOldBlock: ERROR:',
     .  ' are NOT consistent with those specified in',
     .  ' PAO_basis_and_PS_lmax block'
               call die
            endif
            
                if(lmxo(is).gt.lmaxd) overf_l=.true. 
                if(lmxo(is).gt.l_overf) l_overf=lmxo(is) 
                if(.not.overf_l) then 
                   do il=0, lmxo(is)
                      polorb(il,1,is)=0
                   enddo 
                endif
                if(polarization) then  
C****No polarization orbitals if floating Bessel functions***********
                  if(iz(is).ne.-100) then 
c                  if(lmxo(is)+1.gt.lmaxd) overf_l=.true.
c                  if(lmxo(is)+1.gt.l_overf) l_overf=lmxo(is)+1 
                   if(.not.overf_l) polorb(lmxo(is),1,is)=numb_pol 
                  else
                   write(6,'(a)')
     .       'reOldBlock: WARNING polarization option not active with ' 
                   write(6,'(a)')
     .       'reOldBlock: WARNING Z=-100 (Floating Bessel functions)' 
                  endif 
                endif 
                 
                do il=0, lmxo(is)  
                   read(iu,*) l, nzt
                   if(nzt.gt.nzetmx) overf_zet=.true.
                   if(nzt.gt.nzt_overf) nzt_overf=nzt
                      
                   if((.not.overf_zet).and.(.not.overf_l)) then   
                     nzeta(l,1,is)=nzt
                     read(iu,*) (rco(izt,l,1,is),izt=1,nzt)
                     read(iu,*) (lambda(izt,l,1,is),izt=1,nzt)
                   else
                     read(iu,*)
                     read(iu,*)
                   endif 
                enddo  
               
           enddo 
 
           if((overf_zet).or.(overf_l)) then  
              write(6,'(a)') 
     .         'rePAOBasis: ERROR bad dimensions in atmfuncs'
               if(overf_zet) call chkdim('reOldBlock','nzetmx',
     .           nzetmx,nzt_overf,1)
               if(overf_l) call chkdim('reOldBlock','lmaxd',
     .           lmaxd,l_overf,1)
           endif 

          endif        

      end subroutine reOLDblock

      subroutine find_semi(jz,chg,lsemic )
          
           implicit none
           double precision qval(0:4), sum, charge_local, chg
           integer config(0:4), conf_last(0:4), niter,
     .      conf_val(0:4), l, lval, latm, iter,
     .      izlast, nsem, nsemcap, ll, jz, iz_loc, i,
     .      lsemic(0:lmaxd)      
           character sym(0:4)
           
           data sym / 's','p','d','f','g' /


           niter=20
           charge_local=chg
           nsemcap=0              
C Ground state configuration for the atom
              do l=0,lmaxd
                 lsemic(l)=0
              enddo 
              do l=0,3
                 qval(l)=0.0d0
                 config(l)=0
              enddo
              call qvlofz(jz,qval)
              call cnfig(jz,config)
              call lmxofz(jz,lval,latm)
              izlast=jz
              sum=0.0d0 
              do l=0,3
                conf_val(l)=config(l)
                conf_last(l)=config(l)
                sum=sum+qval(l)
              enddo 
              do iter=1,niter

C Restamos la carga de valencia 
                sum=0.0d0 
                do l=0,3
                    sum=sum+qval(l)
                enddo 
                iz_loc=izlast-nint(sum)
                do l=0,3
                   qval(l)=0.0d0
                   config(l)=0
                enddo
                lval=0
                latm=0 
                call qvlofz(iz_loc,qval)
                call lmxofz(iz_loc,lval,latm)
                call cnfig(iz_loc,config)

              nsem=0
              sum=0.0d0
              do l=0,3
                if(conf_last(l)-config(l).ne.0) then 
c                 write(6,*) 'Posible semicore:', config(l),'(',l,')'
                  nsem=nsem+1
                  sum=sum+2*(2*l+1)
                endif
              enddo 

              if((sum.gt.charge_local).and.(nsem.eq.1))then
                  write(6,*) 'ERROR: Expected semicore states:'
                  do l=0,3
                   write(6,*) (i,sym(l),i=config(l),conf_val(l)-1) 
                  enddo

                  call die

              elseif((sum.gt.charge_local).and.(nsem.eq.2)) then 
                  do l=0,3
                     do i=1,conf_last(l)-config(l)
                        if(charge_local.eq.2*(2*l+1) ) then 
                          ll=nint((sum-charge_local-2)/4.0d0)
                          config(ll)=config(ll)+1
                          goto 900
                        endif
                     enddo
                  enddo
                  write(6,*) 'ERROR: Expected semicore states:'
                  do l=0,3
                    write(6,*) (i,sym(l), i=config(l),conf_val(l)-1)
                  enddo
                  call die
               elseif(charge_local.eq.sum) then
                   goto 900
               endif
              nsemcap=max(nsemcap,nsem)
              if (nsem.eq.0) goto 900
              do l=0,3
                conf_last(l)=config(l)
              enddo 
              izlast=iz_loc
              charge_local=charge_local-sum
              enddo 
900           continue
              write(6,'(a,x,15(i1,a1,x))') 
     .        'find_semi: Semicore states:',
     .        ((i,sym(l), i=config(l),conf_val(l)-1),l=0,3)
              do l=0,3
                 lsemic(l)=conf_val(l)-config(l)
c                if(config(l).ne.conf_val(l)) then
c                  write(6,*) (i,sym(l), i=config(l),conf_val(l)-1)
c                endif
              enddo

           end subroutine find_semi
             

           integer function qnumber(a)              
               character(len=*) a
               integer n
               
               if (len(a).ne.3) then
                  write(6,*) 'qnumber: len(a)', len(a), ' ', a
               endif
               if(a(3:3).eq.'1') n=1
               if(a(3:3).eq.'2') n=2
               if(a(3:3).eq.'3') n=3
               if(a(3:3).eq.'4') n=4
               if(a(3:3).eq.'5') n=5
               if(a(3:3).eq.'6') n=6
               if(a(3:3).eq.'7') n=7
               if(a(3:3).eq.'8') n=8
               if(a(3:3).eq.'9') n=9

               qnumber=n

      end function qnumber

      subroutine remass(smass,maxs,ns)

c Reading atomic masses of different species.
c
c Reads fdf block. Not necessarily all species have to be given. The 
c ones not given at input will be assumed to have their natural mass
c (according to atmass subroutine).
c      
c Written by E. Artacho. March 1998.97. 
c ********* INPUT ***************************************************
c integer maxs             : Maximum number of species
c integer ns               : Number of species
c double smass(maxs)       : Atomic masses of each species (amu)
c                            as given by default by atmass
c ********* OUTPUT **************************************************
c double smass(maxs)       : Atomic masses of each species (amu)
c                            modified by explicit numbers at input
c *******************************************************************

      integer           maxs, ns
      double precision  smass(maxs)

      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80
      integer           ni, nn, nr, nv, ist, nst, iu
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)

c check for block and read
!
!     Should be wrapped!!
!
      if ( fdf_block('AtomicMass',iu) ) then

         nst = 0
         do ist = 1, ns+1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if ((ni .ge. 1) .and. (nr .ge. 1)) then
              nst = nst + 1
              if (nst .gt. maxs) then
                 call die('remass: BAD DIMENSIONS. maxs too small')
              endif
              if ((integs(1).gt.ns) .or. (integs(1).lt.1) ) then
                 write(6,"(/,a,i4)") 
     .           'remass: WARNING: Species out of range. Ignored is =',
     .           integs(1)
              else
                 smass(integs(1)) = reals(1)
                 write(6,"(a, i4, a, f12.5)")
     .            'remass: Read atomic mass for species ', integs(1),
     .            ' as ', reals(1)
              endif
            else
              return
            endif
         enddo
         call die('remass: ERROR: Too many entries in AtomicMass')
      endif

 50   continue

      end subroutine remass

      FUNCTION ATMASS(IZ)
      real*8               :: atmass
      integer, intent(in)  :: iz

C Returns the average atomic mass from the atomic number IZ.
C Dta taken from VCH periodic table.
C Written by J.M.Soler. April'97.

      integer, PARAMETER  :: NA=94
      character(len=100) message

      DOUBLE PRECISION AMASS(0:NA)
      DATA AMASS / 0.00,
     .     1.01,  4.00,  6.94,  9.01, 10.81, 12.01, 14.01, 16.00,
     .    19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07,
     .    35.45, 39.95, 39.10, 40.08, 44.96, 47.88, 50.94, 52.00,
     .    54.94, 55.85, 58.93, 58.69, 63.55, 65.39, 69.72, 72.61,
     .    74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22,
     .    92.91, 95.94, 98.91,101.07,102.91,106.42,107.87,112.41,
     .   114.82,118.71,121.75,127.60,126.90,131.29,132.91,137.33,
     .   138.91,140.12,140.91,144.24,146.92,150.36,151.97,157.25,
     .   158.93,162.50,164.93,167.26,168.93,173.04,174.97,178.49,
     .   180.95,183.85,186.21,190.2 ,192.22,195.08,196.97,200.59,
     .   204.38,207.2 ,208.98,208.98,209.99,222.02,223.02,226.03,
     .   227.03,232.04,231.04,238.03,237.05,244.06/


      IF (IZ.LT.0 .OR. IZ.GT.NA) THEN
         write(message,'(a,i4)') 'ATMASS: NO DATA FOR Z =',IZ
         call die(message)
      ELSE
         ATMASS=AMASS(IZ)
      ENDIF
      END function atmass

      end subroutine read_basis

