C $Id: redbasis.f,v 1.7 1999/03/01 17:21:50 wdpgaara Exp $

      subroutine read_basis(nsp,iz,lmxkb,lmxo,nzeta,rco,lambda,
     .         atm_label,polorb,semic,lsemic,charge,smass,
     .         basistype)

C**********************************************************************
C Read information which will used by subroutine ATOM to construct the 
C basis set and the Kleinman-Bylander projectors
C Written by D.Sanchez-Portal. Aug. 1998
C
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
C   INTEGER IZ(NSMAX)   : Atomic number of the different species.
C   INTEGER LMXKB(NSMAX): Maximum angular momentum for the Kleinman-
C                         Bylander projector.
C   INTEGER LMXO(NSAMX) : Maximum angular momentum for the PAO orbitals
C                          (WARNING: this does not include polarization
C                            orbitals)
C   INTEGER NZETA(0:LMAXD,NSAMX): Number of basis functions for each 
C                                  angular momentum (WARNING: this does 
C                                 not include polarization orbitals).
C   REAL*8 RCO(NZETMX,0:LMAXD,NSAMX): Cut-off radius for each orbital
C   REAL*8 LAMBDA(NZETMX,0:LMAXD,NSAMX): Compression factor for each 
C                                           orbital
C   CHARACTER*20 ATM_LABEL(NSMAX):    Species label
C   POLORB(LMAXD)          :  Number of polarization orbitals 
C                             with angular momentum L 
C                             (This orbitals are perturbatively 
C                              caluated from PAO witn L'=L-1)
C   LOGICAL SEMIC          :  Presence or not of semicore orbitals
C   INTEGER LSEMIC         :  Angular momentum for semicore orbitals
C   REAL*8  CHARGE(NSMAX)  :  Charge state of the ion, only for PAO 
C                             calculation purposes.
C   REAL*8  SMASS(NSMAX)  :   Atomic mass for each species.
C   CHARACTER*10 BASISTYPE(NSMAX): Augmentation procedure for each atomic
C                               species
C***********************UNITS*******************************************
C    Distances in Bohr.
C    Energies in Rydbergs.
C***********************************************************************


C***************************PRECISION***********************************

           implicit none

C***********************************************************************
C
C
C**************************PARAMETERS (in file atom.h)******************
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and
C                      projectors.
C INTEGER  NZETMX   : Maximum number of PAO orbitals or polarization 
C                      orbitals with the same angular momentum and for
C                      the same species.
C INTEGER  NRMAX    : Maximum number of points in the functions read
C                      from file '.vps' or '.psatom.data' (this number is
C                      determined by the parameter nrmax in the
C                      program atm, which generates the files with
C                      the pseudopotential information).
C**********************************************************************

        include 'atom.h'
        include 'fdf/fdfdefs.h'

C***********************************************************************
C
C
       integer
     .   nsp, iz(nsmax),lmxkb(nsmax),lmxo(nsmax),polorb(lmaxd,nsmax),
     .   nzeta(0:lmaxd,nsmax), lsemic(nsmax)  

       double precision
     .   rco(nzetmx,0:lmaxd,nsmax),lambda(nzetmx,0:lmaxd,nsmax),
     .   charge(nsmax), smass(nsmax) 
       
       character
     .   atm_label(nsmax)*20, basistype(nsmax)*10 


       logical semic(nsmax)

C
C
C***********************************************************************


C***********************Internal variables******************************
C
C
       integer
     .     numb_pol_default, ns_defect,
     .     ns_read, i, is, l, latm, lmxval
    
       double precision
     .   ql(0:3), atmass, Zval


       character
     .   basis_size*15,basis_size_default*15,
     .   basis_sizes(nsmax)*15, 
     .   basistype_default*10, 
     .   basistype_generic*10
        
       logical
     .   old_block, new_block, chem_block,
     .   old_polorb_default, error
        
       external 
     .   atmass
C
C
C***********************************************************************



C***********INTERNAL PARAMETERS***************************************
C****Default size of the automatic basis set**************************
         parameter(basis_size_default='standard')
C*********************************************************************

C****Default type of basis set*****************************************
         parameter(basistype_default='split')
C*********************************************************************

C****Default polarization orbitals (parameter for the old input)******
        parameter(old_polorb_default=.false.)
C*********************************************************************

C*Default number of polarization orbitals (parameter for the old input)*
        parameter(numb_pol_default=1)
C**********************************************************************

      
C******Checking the total number or species**********************
 
      ns_defect = 0
      ns_read = fdf_integer('NumberOfSpecies',ns_defect)
      if (ns_read.ne.nsp) then
        write(6,'(2a)')
     . 'read_basis: ERROR: Read number of species different from'
        write(6,'(a)')
     . 'read_basis: ERROR: that provided by the input of the routine'
        stop
       endif

       if (nsp.gt.nsmax) then 
          write(6,'(2a,i4)')
     . 'read_basis: ERROR: parameter nsmax in file atom.h ',
     . 'must be incremented to at least ', nsp
        stop 
       endif 

C**********Checking if is there any block containing information*****
C********************about the basis set******************************
           old_block=fdf_block('PAO_basis_and_PS_lmax',i)
           new_block=fdf_block('PAO.Basis',i) 

           if((old_block).and.(new_block)) then 
             write(6,'(2a)')
     .  'read_basis: ERROR both information blocks: ', 
     .      'PAO.Basis and PAO_basis_and_PS_lmax'
             write(6,'(2a)')
     .  'read_basis: ERROR cannot be simultaneously specified ', 
     .      'in the input file'
             stop
           endif 

C**********Exists the Chemical_species_labels block ?*****************
          chem_block=fdf_block('Chemical_species_label',i)

C********The block Chemical_species_label should exist always***********
          if (.not.chem_block) then 
          write(6,100)
          write(6,101) 
          write(6,'(a)')
     .   'read_basis: ERROR: Chemical_species_label block not found'
          write(6,'(a)')
     .   'read_basis: ERROR: The chemical species must be specified.' 
          stop 
          endif 

C**********Read species and labels********************************** 
C
          call rechemsp(atm_label,iz,nsmax,nsp)
C
C************************************************************************


C**Standard valence charge for each atom and checking for semicore orbitals**
C
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
            call semicore(Zval,atm_label(is),semic(is),lsemic(is))
            if (semic(is)) lmxval=max(lmxval,lsemic(is))
            if(lmxval.gt.lmaxd) then 
              write(6,'(a)') 
     .  'read_basis: ERROR bad dimensions in atom.h'
              call chkdim('read_basis','lmaxd',
     .        lmaxd,lmxval,1) 
            endif 

            if ( semic(is) ) then
               if ( ql(lsemic(is)).gt.1.0d-4 ) then 
                  write(6,'(a,i3)')
     .         'read_basis: ERROR a semicore shell with l=',lsemic(is)
                  write(6,'(3a)')
     .         'read_basis: ERROR cannot be included (by the moment) ',
     .         'in the valence for ', atm_label(is)
                  write(6,'(2a)')
     .         'read_basis: ERROR because there are also populated ',
     .         'valence states with the same angular momentum'  
                  stop 
               endif 
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
         call remass(smass,nsmax,nsp)
C
C************************************************************************


C****Augmentation procedure for the construction of the basis set********
C
         basistype_generic=
     .        fdf_string('PAO.BasisType',basistype_default)

         call type_name(basistype_generic) 
         do is=1, nsp
             basistype(is)=basistype_generic
         enddo 

C
C************************************************************************


           if(.not.old_block) then 
C******** DESIRED SIZE FOR THE BASIS SET********************************** 
C
           basis_size=fdf_string('PAO.BasisSize',basis_size_default)
           call size_name(basis_size)
           do is=1,nsp
                basis_sizes(is)=basis_size  
           enddo  
           call resizes(atm_label,basis_sizes,nsmax,nsp)  
C
C************************************************************************


          
C********Maximum angular momentum for the Kleinman-Bylander projectors**** 
C
           do is=1,nsp
               lmxkb(is)=-1
           enddo 
           call relmxkb(atm_label,lmxkb,nsmax,nsp)  

           endif 
C
C************************************************************************




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
           do is=1,nsmax
              do l=0,lmaxd
                  nzeta(l,is)=0
                  if(l.gt.0) polorb(l,is)=0
              enddo 
           enddo
C
C************************************************************************



C****Reading from new block****************************************
C
           if((.not.old_block).and.(new_block)) then 

           call rePAOBasis(atm_label,iz,lmxo,nzeta,rco,lambda,
     .         polorb,charge,basistype,nsp)

           endif 
C
C************************************************************************

C*******************Reading from old block******************************
C           
          if(old_block) then   
              call reOldBlock(iz,lmxo,lmxkb,nzeta,rco,lambda,
     .        polorb,charge,nsp)
          endif 
C
C************************************************************************


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
               stop
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
               stop
             endif   
            
             if (lmxo(is).eq.-1) then 
               call autobasis(atm_label(is),basis_sizes(is),
     .           abs(iz(is)),
     .           semic(is), lsemic(is),
     .           lmxo(is),nzeta(0,is),rco(1,0,is),lambda(1,0,is),
     .           polorb(1,is)) 
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
             if(semic(is)) lmxval=max(lmxval,lsemic(is))
             if (lmxval.gt.lmxo(is)) then 
               write(6,'(/,2a)')
     .     'read_basis: ERROR: For species ', atm_label(is)
               write(6,'(a)')
     .     'read_basis: ERROR: Maximum angular momentum of the basis' 
               write(6,'(a,i4)')
     .     'read_basis: ERROR: orbitals must be at least: ',lmxval 
               error=.true.                
             endif  

             if(nzeta(0,is).eq.0) then 
                write(6,'(a)')
     .     'read_basis: ERROR: All species should have at least one'
                write(6,'(a)')
     .     'read_basis: ERROR: s (l=0) shell in their basis set'
                error=.true.    
             endif 


             do l=1,lmxval
                if(
     .           (ql(l).gt.1.0d-4).or.
     .                   ((l.eq.lsemic(is)).and.(semic(is))) ) then 

                  if(nzeta(l,is).eq.0) then 
                     write(6,'(2a)')
     .    'read_basis: ERROR: For species ', atm_label(is)
                     write(6,'(a)')
     .    'read_basis: ERROR: There must be at least one shell of ' 
                     write(6,'(a,i4)')
     .    'read_basis: ERROR: basis set orbitals with l=',l
                      error=.true.    
                  endif
               endif 
             enddo  
            endif 
C
C************************************************************************          
C************************************************************************
C 
             if(lmxkb(is).eq.-1) lmxkb(is)=lmxo(is)+1   
             if(iz(is).lt.0) lmxkb(is)=0
C
C************************************************************************

          enddo  

          if(error) stop
C
C************************************************************************
 
            

100   format(/,'redbasis: ',71(1h*))
101   format('redbasis:                  INPUT ERROR')
102   format('rebasis: ',71(1h*))

      end 
