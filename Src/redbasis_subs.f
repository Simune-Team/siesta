C $Id: redbasis_subs.f,v 1.7 1999/02/26 21:19:21 daniel Exp $

C This file contains all the subroutines called by subroutine
C 'redbasis' to obtain all the input information for the generation 
C of the basis set and Kleinman-Bylander projector.
C
C  Written by D. Sanchez-Portal, Aug. 1998.
C***********************************************************************
c     subroutine size_name
c     subroutine resizes
c     subroutine rechemsp
c     subroutine relmxkb
c     subroutine semicore
c     subroutine autobasis
c     subroutine type_name
c     subroutine rePAOBasis
c     subroutine reOldBlock 
C***********************************************************************


      subroutine size_name(basis_size)

C  Written by D. Sanchez-Portal, Aug. 1998.

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
         stop
      endif

      return 
      end


      subroutine resizes(atm_label,basis_sizes,maxs,ns)

c *******************************************************************
c Reading atomic basis sizes for different species.
c
c Reads fdf block. Not necessarily all species have to be given. The 
c ones not given at input will be assumed to have the basis sizes 
c given by the general input PAO.BasisSize, or its default value. 
c      
c Written by D. Sanchez-Portal. Aug. 1998. 
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

      implicit          none
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


c enable FDF input/output

      include 'fdf/fdfdefs.h'

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
                     write(6,"(/,a)")
     .                'resizes: ERROR: BAD DIMENSIONS. Too small maxs'
                   stop 'resizes: ERROR: BAD DIMENSIONS. Too small maxs'
                   endif
 
              endif  


               
              else
                 return
              endif 
          

 
            else
              return
            endif
         enddo

         write(6,'(a)') 
     .        'resizes: ERROR: Too many entries in PAO.BasisSizes'
         stop 'resizes: ERROR: Too many entries in PAO.BasisSizes'

      endif

 50   continue

      return
      end

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

      implicit          none
      integer           maxs, ns, iz(maxs)
      character*20      atm_label(maxs) 
      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80
      integer           ni, nn, nr, nv, ist, nst, iu, nsread
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)


c enable FDF input/output

      include 'fdf/fdfdefs.h'

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
                 write(6,"(/,a)") 
     .                'rechemsp: ERROR: BAD DIMENSIONS. Too small maxs'
                 stop 'rechemsp: ERROR: BAD DIMENSIONS. Too small maxs'
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
               stop 'rechemsp: ERROR: There are some undefined species'
              endif  

              return
            endif 


         enddo

         write(6,'(a)') 
     .     'rechemsp: ERROR: Too many entries in Chemical_species_label'
      stop 'rechemsp: ERROR: Too many entries in Chemical_species_label'

      endif

 50   continue

      return
      end





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


c enable FDF input/output

      include 'fdf/fdfdefs.h'

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
                     write(6,"(/,a)")
     .                'relmxkb: ERROR: BAD DIMENSIONS. Too small maxs'
                  stop 'relmxkb: ERROR: BAD DIMENSIONS. Too small maxs'
                   endif

              endif 



            else
              return
            endif
         enddo

         write(6,'(a)') 
     .        'relmxkb: ERROR: Too many entries in AtomicMass'
         stop 'relmxkb: ERROR: Too many entries in AtomicMass'

      endif

 50   continue

      return
      end




        subroutine semicore(Zval,atm_label,semic,lsemic)
c *******************************************************************
c Reading information for basis generation.
c
c Reading the valence charge provided by the pseudopotential file
c it can be deduced if any semicore internal states have to be included 
c into the valence shell.
c
c Written by D. Sanchez-Portal. Aug. 1998.
c ********* INPUT ***************************************************
c double precision Zval          : Expected valence charge of the atom
c character*20 atm_label         : Labels for the different species
c ********* OUTPUT **************************************************
c logical      semic             : if .true. there are semicore states 
c                                  if .false. there are not.
c integer      lsemic            : Angular momentum of the semicore shell
c *******************************************************************

       implicit none      

       double precision Zval, b, a, Zval_read, tiny
       character atm_label*20, paste*50, fname*50, namatm*2, 
     .   icorr*2, irel*3, nicore*4,
     .   method(6)*10, text*70
       integer lsemic, nr, npotd, npotu, i,ndiff
       logical semic,found
       
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
             stop 'semicore: ERROR: Pseudopotential file not found'
           endif

           open(unit=1,file=fname,form='unformatted',status='unknown')

           read(1) namatm, icorr, irel, nicore,
     .     (method(i),i=1,6), text,
     .     npotd, npotu, nr, b, a, Zval_read
C*******************Close the file******************************************
             close(1)
C***************************************************************************
             if(dabs(Zval-Zval_read).lt.tiny) then 
                   semic=.false.
                   return 
             else 

             ndiff=nint(dabs(Zval-Zval_read))
             if(dabs(ndiff-dabs(Zval-Zval_read)).gt.tiny) then 
                write(6,'(2a)') 
     .   'semicore: ERROR expected valence charge for species ', 
     .    atm_label
                write(6,'(2a)')
     .  'semicore: ERROR and the value read from the file ', fname
                write(6,'(a,f6.3,a,f6.3)')
     .  'semicore: ERROR differ:  Zval(expected)= ', Zval,
     .  ' Zval(read)= ',Zval_read 
             stop
             endif 

             semic=.true.
             lsemic=(ndiff-2)/4 

             if ((4*lsemic+2).ne.ndiff) then 
                 write(6,'(a)')
     .   'semicore: ERROR at the moment the programme can only deal'
                 write(6,'(a)')
     .   'semicore: ERROR with one shell of semicore orbitals'
                 write(6,'(2a)') 
     .   'semicore: ERROR expected valence charge for species ', 
     .    atm_label
                 write(6,'(2a)')
     .  'semicore: ERROR and the value read from the file ', fname
                 write(6,'(a,f6.3,a,f6.3)')
     .  'semicore: ERROR differ:  Zval(expected)= ', Zval,
     .  ' Zval(read)= ',Zval_read  
             stop 
            endif 
              
            write(6,'(a,i3,2a)')
     .       'semicore: l=',lsemic,
     .       ' semicore shell included in the valence for species ',
     .        atm_label
       
 
             endif 
     

         return
         end
       
     




        subroutine autobasis(atm_label,basis_size,iz,
     .       semic, lsemic,
     .       lmxo,nzeta,rco,lambda,
     .       polorb)
         
C Generates automatic basis set from the information given 
C in the block PAO.BasisSizes
C Written by D. Sanchez-Portal, Aug. 1998.

         implicit none
         include 'atom.h'
        
         character atm_label*20, basis_size*15
         integer lmxo, nzeta(0:lmaxd), polorb(lmaxd),
     .     lsemic, iz
         logical semic 
         double precision rco(nzetmx,0:lmaxd), lambda(nzetmx,0:lmaxd)
       
C*******************internal variables*********************************
         
         double precision ql(0:3)
         integer lmxval, latm, l, nzt, izt, lmax
         logical pol_contr, overflow
         
         overflow=.false.
         call lmxofz(iz,lmxval,latm)
         call qvlofz(iz,ql) 
         
         lmax=lmxval
         if(semic) lmax=max(lsemic,lmxval) 
         if(lmax.gt.lmaxd) overflow=.true.
        
         call size_name(basis_size)  
         if((basis_size.eq.'dzp').or.(basis_size.eq.'dz') ) then  
               izt=nzetmx
               if(izt.lt.2)  write(6,'(/2a)') 
     .        'autobasis: ERROR: Parameter nzetmx in file atom.h ',
     .        'must be increased to 2'       
         endif

         pol_contr=.false. 
   
         do l=0,lmxval
             
             
             if(dabs(ql(l)).lt.1.0d-4) then  

                  if(l.eq.0) then 
                     if ((basis_size.eq.'sz').or.
     .                    (basis_size.eq.'szp')) nzt=1

                     if ((basis_size.eq.'dz').or.
     .                    (basis_size.eq.'dzp')) nzt=2   

                     if(.not. overflow) then 
                       nzeta(l)=nzt
                       do izt=1,nzt
                           rco(izt,l)=0.0d0
                           lambda(izt,l)=0.0d0
                        enddo
                     endif 

                  elseif(  ((basis_size.eq.'szp').or.
     .                      (basis_size.eq.'dzp'))
     .                        .and. (.not.pol_contr)
     .                        .and. (l.ne.0)) then 
                          pol_contr=.true.   
                          if(.not. overflow) then 
                             polorb(l)=1
                             nzeta(l)=0  
                          endif 
                  elseif( (pol_contr)
     .                        .and. (l.ne.0)) then 
                          if(.not. overflow) then
                             nzeta(l)=0  
                          endif 
                  endif
 
                  if((semic).and.(lsemic.eq.l)) then 
                      if(.not. overflow) then
                          nzeta(lsemic)=1
                          rco(1,lsemic)=0.0d0 
                          lambda(1,lsemic)=0.0d0 
                      endif
                  endif 


             endif 
            
             if(ql(l).gt.1.0d-4) then   


                if(semic) then
                      if(lsemic.eq. l) then 
                       write(6,'(a,i3)')
     .  'autobasis: ERROR a semicore shell with l=',lsemic
                                       write(6,'(3a)')
     .  'autobasis: ERROR cannot be included (by the moment) ',
     .     'in the valence for ', atm_label
                                       write(6,'(2a)')
     .  'autobasis: ERROR because there are also populated ',
     .    'valence states with the same simmetry'
                          stop
                      endif
                endif 


                if ((basis_size.eq.'sz').or.
     .                    (basis_size.eq.'szp')) nzt=1 
    
                if ((basis_size.eq.'dz').or.
     .                    (basis_size.eq.'dzp')) nzt=2
                if(.not. overflow) then 
                   nzeta(l)=nzt
                   do izt=1,nzt
                       rco(izt,l)=0.0d0 
                       lambda(izt,l)=0.0d0 
                   enddo 
                endif 

             endif 
                         
           enddo      

           if ((semic).and.(lsemic.gt.lmxval)) then    
                 if(.not. overflow) then
                     do l=lmxval+1,lsemic
                        nzeta(lsemic)=0
                     enddo  
                     nzeta(lsemic)=1
                     rco(1,lsemic)=0.0d0
                     lambda(1,lsemic)=0.0d0
                 endif  
           endif  

           if((basis_size.eq.'szp').or.(basis_size.eq.'dzp') ) then
            if(.not.pol_contr) then 
              lmxval=lmxval+1  
              lmax=max(lmax, lmxval) 
              if(lmax.gt.lmaxd) overflow=.true. 
              if(.not.overflow) then 
                 polorb(lmxval)=1 
              endif 
            endif 
           endif

           lmxo=lmax  
         
           if(overflow) then 
               write(6,'(/a)')
     .          'autobasis: ERROR: Parameter lmaxd in file atom.h ' 
               write(6,'(a,i4)')
     .          'autobasis: ERROR: must be increased to at least ',lmax 
               stop
           endif
    


         return
         end  



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
              stop
         endif



        return  
        end

         subroutine rePAOBasis(atm_label,iz,lmxo,nzeta,rco,lambda,
     .         polorb,charge,basistype,ns)

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
      include 'atom.h'
 
      integer   ns, lmxo(nsmax), polorb(lmaxd,nsmax)
      integer   nzeta(0:lmaxd,nsmax), iz(nsmax)
      double precision
     .   rco(nzetmx,0:lmaxd,nsmax), lambda(nzetmx,0:lmaxd,nsmax),
     .   charge(nsmax)
      character atm_label(nsmax)*20, basistype(nsmax)*10
      external   parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80, label_read*20
      integer           ni, nn, nr, nv, ist, iu,  ns_read, ilabel
      integer           integs(4), lastc, lc(0:3), nzt, izt, l, ipol 
      integer           nsh,l_control(0:lmaxd), l_read, ish, ish2
      integer           nzt_overf, l_overf, lmax, l_contr
      double precision  reals(4), values(4)
      logical           overf_zet, overf_l

c enable FDF input/output

      include 'fdf/fdfdefs.h'

c check for block and read

      if ( fdf_block('PAO.Basis',iu) ) then
         
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


C**Reading compression factors (if presents) for the last***********
C**orbital of the previous chemical specie************************** 
             if((ist.ne.1).and.(nr.ge.nzt).and.(ni.eq.0)) then 
                do izt=1,nzt 
                    lambda(izt,l_read,ns_read)=reals(izt)
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
                lmax=0
                do ish=1,nsh
                  read(iu,'(a)') line
                  lastc = index(line,'#') - 1
                  if (lastc .le. 0) lastc = len(line)
                  call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
                  if((ish.ne.1).and.(nr.ge.nzt)) then 
                      do izt=1,nzt 
                        lambda(izt,l_read,ns_read)=reals(izt)
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
     .          'rePAOBasis: ERROR duplicate information for species ',
     .                atm_label(ns_read)    
                      stop
                     endif 
                   enddo 
                   l_control(ish)=l_read
 
                   nzt=integs(2)
                   if(l_read.gt.lmax)  lmax=l_read  
                   if(nzt.gt.nzetmx) overf_zet=.true.
                   if(nzt.gt.nzt_overf) nzt_overf=nzt
                   if(l_read.gt.lmaxd) overf_l=.true. 
                   if(l_read.gt.l_overf) l_overf=l_read
                   if(.not.overf_l) nzeta(l_read,ns_read)=nzt           


                   if((nn.gt.0).and.(names(1:1).eq.'P') ) then 
 
C****No polarization option if floating Bessel functions (Z=-100)
                    if(iz(ns_read).ne.-100) then

                      if(l_read+1.gt.lmaxd) overf_l=.true. 
                      if(l_read+1.gt.l_overf) l_overf=l_read+1
                      
                      if(ni.ge.3) then 
                        ipol=integs(3) 
                      else
                        ipol=1
                      endif   
                      
                      if(.not.overf_l) polorb(l_read+1,ns_read)=ipol
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
 
               if((.not.overf_zet).and.(.not.overf_l)) then  
                   read(iu,'(a)') line 
                   lastc = index(line,'#') - 1
                   if (lastc .le. 0) lastc = len(line) 
           call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )

                   if(nr.lt.nzt) goto 40 
                   do izt=1,nzt 
                      rco(izt,l_read,ns_read)=reals(izt) 
                      lambda(izt,l_read,ns_read)=0.0d0 
                   enddo 
                endif

30              continue        
                enddo  
                lmxo(ns_read)=lmax 

                do l=0,lmax
                    l_contr=0
                    do ish=1,nsh
                        if(l.eq.l_control(ish)) l_contr=1
                    enddo 
                    if(l_contr.eq.0) nzeta(l,ns_read)=0
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
     .      'rePAOBasis: ERROR bad dimensions in atom.h'
            if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .      nzetmx,nzt_overf,1)
            if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .      lmaxd,l_overf,1)

         endif 
         
         write(6,'(a)') 
     .     'rePAOBasis: ERROR: Too many entries in PAO.Basis'
       stop  'rePAOBasis: ERROR: Too many entries in PAO.Basis'




      endif     
 40   continue  


        if((overf_zet).or.(overf_l)) then  
           write(6,'(a)') 
     .     'rePAOBasis: ERROR bad dimensions in atom.h'
           if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .     nzetmx,nzt_overf,1)
          if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .    lmaxd,l_overf,1)
        endif 


       write(6,'(a)')
     .    'rePAOBasis: ERROR reading PAO.Basis block'
       write(6,'(a)')
     .    'rePAOBasis: ERROR check the block sintaxis'
        stop


 50   continue 


       if((overf_zet).or.(overf_l)) then   
               write(6,'(a)') 
     .      'rePAOBasis: ERROR bad dimensions in atom.h'
         if(overf_zet) call chkdim('rePAOBasis','nzetmx',
     .          nzetmx,nzt_overf,1)
         if(overf_l) call chkdim('rePAOBasis','lmaxd',
     .          lmaxd,l_overf,1)

           stop 
        endif 


      return
      end

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
c *******************************************************************

      implicit  none 
      include 'atom.h'
 
      integer   ns, lmxo(nsmax), polorb(lmaxd,nsmax)
      integer   nzeta(0:lmaxd,nsmax), iz(nsmax), lmxkb(nsmax)
      double precision
     .   rco(nzetmx,0:lmaxd,nsmax), lambda(nzetmx,0:lmaxd,nsmax),
     .   charge(nsmax)
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
       
c enable FDF input/output

      include 'fdf/fdfdefs.h'



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
             stop
            endif
            
                if(lmxo(is).gt.lmaxd) overf_l=.true. 
                if(lmxo(is).gt.l_overf) l_overf=lmxo(is) 
                if(polarization) then  
C****No polarization orbitals if floating Bessel functions***********
                  if(iz(is).ne.-100) then 
                   if(lmxo(is)+1.gt.lmaxd) overf_l=.true.
                   if(lmxo(is)+1.gt.l_overf) l_overf=lmxo(is)+1 
                   if(.not.overf_l) polorb(lmxo(is)+1,is)=numb_pol 
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
                     nzeta(l,is)=nzt
                     read(iu,*) (rco(izt,l,is),izt=1,nzt)
                     read(iu,*) (lambda(izt,l,is),izt=1,nzt)
                   else
                     read(iu,*)
                     read(iu,*)
                   endif 
                   if(il.gt.0) polorb(il,is)=0
                enddo  
               
           enddo 
 
           if((overf_zet).or.(overf_l)) then  
              write(6,'(a)') 
     .         'rePAOBasis: ERROR bad dimensions in atom.h'
               if(overf_zet) call chkdim('reOldBlock','nzetmx',
     .           nzetmx,nzt_overf,1)
               if(overf_l) call chkdim('reOldBlock','lmaxd',
     .           lmaxd,l_overf,1)
           endif 

          endif        


      return
      end

