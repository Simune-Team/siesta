      module interfaces

      interface
      subroutine hsparse( negl, cell, nsc, na, isa, xa,
     .                    lasto, lastkb, iphorb, iphkb,
     .                    nlhmax, numh, 
     .                    listhptr, listh, Node, Nodes )
      integer, intent(in)          ::  na
      integer, intent(in)          ::  iphkb(:), iphorb(:)
      integer, intent(in)          ::  isa(na), lastkb(0:na),
     $                                 lasto(0:na)
      integer, intent(in)          ::  nsc(3), Node, Nodes
      double precision, intent(in) ::  cell(3,3), xa(3,na)
      logical, intent(in)          ::  negl

      integer, intent(inout)       ::  nlhmax
      integer, intent(out)         ::  listh(:), listhptr(:)
      integer, intent(out)         ::  numh(:)

      end subroutine hsparse

      subroutine ordern(usesavelwf,ioptlwf,natoms,nbasis,lasto,
     .                  isa,qa,rcoor,rh,cell,xa,iscf,istep,itmax,
     .                  ftol,eta,enum,maxuo,maxnh,numh,listhptr,
     .                  listh,h,s,chebef,noeta,rcoorcp,beta,ipcheb,
     .                  dm,edm,Ecorrec)

      integer, intent(INOUT)          :: istep   !!  ??????
      double precision, intent(inout) :: eta

      integer, intent(in)           ::
     .  ioptlwf, ipcheb, iscf, itmax, natoms, 
     .  nbasis, maxnh, maxuo

      integer, intent(in)           ::
     .  isa(natoms), lasto(0:natoms), listh(maxnh), 
     .  numh(:), listhptr(:)

      double precision, intent(out) :: dm(maxnh), Ecorrec, edm(maxnh)
      double precision, intent(in)  ::
     .  cell(3,3), 
     .  enum, ftol, h(maxnh), qa(natoms), rcoor, rh, 
     .  s(maxnh), xa(3,natoms), beta, rcoorcp

      logical, intent(in)           :: chebef, noeta, usesavelwf

      end subroutine ordern

      end interface

      end module interfaces

