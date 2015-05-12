#!/usr/bin/env python

# Should make it "backwards" compatible down to 2.6
from __future__ import print_function

# This small utility processes data from tbtrans and
# writes easy to manage files.

# It has been created so that it can easily be combined
# with other codes by a simple import statement.
# PLEASE
# !!!
# Do not edit this file without contributing 
# your edits to the community!
# What is helping you, could MOST likely also
# help others in the community!
# Respect the maintainer, the code and the other users.
# !!!

# When retrieving data from the file it will return
# so in units of the file.
# The standard units of SIESTA/TBtrans is Ry/Bohr
# To convert utilize the intrinsic conversion
# factors.
#  E = TBTFile(<fname>).E * TBTFile.Ry
# which brings the energy array from Rydberg to electron volt.

# Written by: Nick Papior Andersen, 2014

# Load needed modules
import itertools
import numpy as np
import netCDF4 as nc

# Get current date
import datetime
_TODAY = datetime.date.today().isoformat()

# A class to retain information given in a 
# regular <syslabel>.TBT.nc file.
class TBTFile(object):

    # constant fields for converting Ry/Bohr to eV/Ang
    Ry   = 13.60580
    Bohr = 0.529177

    """ TBT file to contain the NetCDF file """
    def __init__(self,file):
        """ Initialised a <syslabel>.TBT.nc file to handle data """
        self.file = file
        # Open the file
        self.nc = nc.Dataset(self.file,'r')
        # In case we want to do k-avg it will be easier to
        # have the k-weights locally.
        try:
            self.wk = np.array(self.nc.variables['wkpt'])
        except:
            # A Gamma-point calculation
            self.wk = np.array([1.])

    def date(self):
        """ Returns the creation date of the NetCDF file """
        return self.nc.date

    @property
    def cell(self):
        """ Returns the unit cell """
        return self._get_Data('cell')
    @property
    def xa(self):
        """ Returns the atomic coordinates """
        return self._get_Data('xa')

    @property
    def xa_dev(self):
        """ Returns the atomic coordinates only for the device """
        return self._get_Data('xa')[self.a_dev,:]

    @property
    def kpt(self):
        """ Returns the k-points """
        return self._get_Data('kpt')
    @property
    def wkpt(self):
        """ Returns the k-point weights """
        return self._get_Data('wkpt')

    @property
    def E(self):
        """ Returns the energies that resides in the file """
        return self._get_Data('E')

    @property
    def na_u(self):
        """ Returns number of atoms in the cell """
        return int(len(self.nc.dimensions['na_u']))

    @property
    def na(self):
        """ Returns number of atoms in the device region """
        return int(len(self.nc.dimensions['na_d']))

    @property
    def no_u(self):
        """ Returns number of orbitals in the unit-cell """
        return int(len(self.nc.dimensions['no_u']))
    @property
    def no(self):
        """ Returns the size of the projected region """
        return int(len(self.nc.dimensions['no_d']))

    @property
    def a_dev(self):
        """ Returns the atoms in the device region """
        return np.array(self.nc.variables['a_dev'][:],np.int)

    @property
    def lasto(self):
        """ Returns the last orbital of the equivalent atom """
        return self._get_Data('lasto')

    @property
    def pivot(self):
        """ Returns the pivot table """
        return self._get_Data('pivot')

    def o2a(self,orbs):
        """ Returns an array of of atoms for where the orbitals
        reside.
        """
        if isinstance(orbs,(int,np.int16,np.int32)):
            return np.where(orbs < self.lasto)[0][0]
        return np.array([np.where(o < self.lasto)[0][0] for o in orbs])

    def a2o(self,atoms):
        """ Returns an array of integers corresponding 
        to the orbitals of the input atoms """
        tatm = np.array(atoms)
        lasto = np.zeros((self.na_u+1,),np.int)
        lasto[1:] = self.lasto
        # Count the new size
        a_no = np.sum(lasto[tatm] - lasto[tatm-1])
        orbs = np.zeros((a_no,),np.int)
        i = 0
        for atm in atoms:
            n = lasto[atm]-lasto[atm-1]
            orbs[i:i+n] = np.arange(lasto[atm-1]+1,lasto[atm]+1)
            i += n
        return orbs

    def a2idx(self,atoms):
        """ Returns an array of integers corresponding 
        to the orbitals of the input atoms """
        return self.o2idx(self.a2o(atoms))

    def o2idx(self,orbs):
        """ Returns an array of integers corresponding 
        to the orbitals of the input orbitals """
        # Locate the equivalent orbitals in
        # the pivot table and create an index variable
        return np.where(np.in1d(self.pivot,orbs))[0]

    def _get_Data(self,var,tree=[],k_avg=False):
        """ Generic routine for retrieving the data in a tree """
        g = self.nc # just a reference copy
        if tree:
            for t in tree:
                if t: g = g.groups[t]
        return self._kavg(np.array(g.variables[var]),k_avg)

    @property
    def elecs(self):
        """ Returns a list of electrodes """
        elecs = self.nc.groups.keys()
        # in cases of not calculating all 
        # electrode transmissions we
        # find the last one
        tmp_elecs = self.nc.groups[elecs[0]].variables.keys()
        for tmp in tmp_elecs:
            if tmp.endswith('.T'):
                tmp = tmp.split('.')[0]
            else: continue
            if not tmp in elecs: elecs.append(tmp)
        return elecs

    def _kavg(self,DATA,k_avg):
        """ Returns the k-averaged part of the data """
        tmp = DATA.shape
        if isinstance(k_avg,bool):
            if k_avg:
                DATA.shape = (len(self.wk),-1)
                D = np.sum(DATA[:,:] * self.wk[:,None],axis=0)
                D.shape = tmp[1:]
                return D
        else:
            DATA.shape = (len(self.wk),-1)
            D = DATA[k_avg,:] * self.wk[k_avg,None]
            if len(D.shape) > 1:
                D = np.sum(D,axis=0)
            D.shape = tmp[1:]
            return D
        return DATA

    def has_Elec(self,El):
        """ Returns true if the Group 'El' exists in the top level """
        return El in self.nc.groups.keys()

    def T(self,E1,E2=None,k_avg=True):
        """ Returns the transmission function between electrodes
        E1 and E2 """
        if not E2: E2 = E1
        if not E1 in self.nc.groups:
            # Swap the electrodes, the user may have them confused
            E1, E2 = E2, E1

        if E1 == E2:
            # The user have requested the "reflection"
            return self._get_Data(E1+'.R',k_avg=k_avg)
        return self._get_Data(E2+'.T',[E1],k_avg=k_avg)

    def DOS(self,El=None,k_avg=True):
        """ Returns the DOS """
        if El: return self.ADOS(El,k_avg=k_avg)
        return self._get_Data('DOS',k_avg=k_avg)

    def ADOS(self,El,k_avg=True):
        """ Returns the DOS from the spectral function """
        return self._get_Data('ADOS',[El],k_avg=k_avg)

    def mu(self,El):
        """ Returns the chemical potential of the electrode """
        return self._get_Data('mu',[El])

    def kT(self,El):
        """ Returns the temperature of the electrode """
        return self._get_Data('kT',[El])

    def Jij(self,El,k_avg=True,E=None):
        """ Returns the orbital current of the electrode in a scipy
        csr matrix format """
        # Create csr sparse formats.
        # We import here as the user might not want to
        # rely on this feature.
        from scipy.sparse import csr_matrix
        
        # First read in the sparse matrix information
        col = np.array(self.nc.variables['list_col'][:],np.int) - 1
        tmp = np.cumsum(self.nc.variables['n_col'][:])
        s = len(tmp)
        ptr = np.zeros((s+1,),np.int)
        ptr[1:] = tmp[:]
        del tmp
        if E is None:
            J = self._get_Data('J',[El],k_avg=k_avg)
        else:
            # if the user requests a specific energy-point, then we
            # only return the csr matrix
            J = self._kavg(np.array(self.nc.groups[El].variables['J'][:,E,:]),k_avg)
            return csr_matrix((J[:],col,ptr),shape=(s,s))
        if len(J.shape) == 2:
            return J, csr_matrix((J[0,:],col,ptr),shape=(s,s))
        return J, csr_matrix((J[0,0,:],col,ptr),shape=(s,s))

    def Jij2Ja(self,Jij):
        """ Returns the total current flowing through each atom.
        
        It requires a sparse matrix return from ``self.Jij`` that
        contains the orbital current.

        It returns the atomic current on all atoms in the device
        region.
        Note that the returned atomic current is for the all atoms
        even those not in ``self.a_dev.``
        """
        # Convert to csr format (just ensure it)
        tmp = Jij.tocsr()
        Ja = np.zeros((self.na_u,),np.float)
        atoms = self.a_dev
        lasto = np.append([0],self.lasto)
        for ia in atoms:
            Ja[ia-1] = np.sum(np.sqrt((tmp[lasto[ia-1]:lasto[ia],:].multiply(
                            tmp[lasto[ia-1]:lasto[ia],:])).data))
        return Ja
    
    def Jij2Jab(self,Jij,symmetry=True):
        """ Returns the total current flowing through each atom.
        
        It requires a sparse matrix return from ``self.Jij`` that
        contains the orbital current.
        
        If ``symmetry`` is ``True`` it will only return one of the 
        bond currents, the one which is positive.

        It returns the atomic current on all atoms in the device
        region.
        Note that the returned atomic current is for the all atoms
        even those not in ``self.a_dev``, they are just zero outside
        the device region.
        """
        # Create csr sparse formats.
        # We import here as the user might not want to
        # rely on this feature.
        from scipy.sparse import lil_matrix

        # Convert to lil format (just ensure it)
        Jab = lil_matrix((self.na_u,self.na_u))

        # Faster to loop across data
        tmp = Jij.tocoo()
        if symmetry:
            for ja,ia,d in zip(self.o2a(tmp.row),self.o2a(tmp.col),tmp.data*.5):
                if ia == ja: continue
                if d > 0:
                    Jab[ja,ia] += d
                else:
                    Jab[ia,ja] -= d
        else:
            for ja,ia,d in zip(self.o2a(tmp.row),self.o2a(tmp.col),tmp.data):
                if ia == ja: continue # it is zero anyway
                Jab[ja,ia] += d
        return Jab.tocsr()

class TBTProjFile(TBTFile):
    """ We inherit the TBT.nc file object. Many of the quantities are similar """

    def mols(self):
        """ Returns a list of all projection molecules """
        groups = []
        for group in self.nc.groups.keys():
            if len(self.nc.groups[group].groups) > 0:
                # this is a group
                groups.append(group)
        return groups

    def projs(self,mol):
        """ Return all projections of this molecule """
        return list(self.nc.groups[mol].groups.keys())

    def lvl(self,mol,relative=True):
        """ Returns levels of projection molecule """
        if relative:
            return np.array(self.nc.groups[mol].variables['lvl'][:],np.int)
        # The user has requested the absolute projection states
        lvl = self.lvl(mol)
        # Correct for HOMO-index
        H_idx = int(self.nc.groups[mol].HOMO_index)
        return np.where(lvl < 0 , lvl, lvl - 1) + H_idx

    def eig(self,mol,lvl=None,k_avg=False):
        """ Returns the eigenvalues of the molecule """
        a = np.asarray(self._get_Data('eig',tree=[mol],k_avg=k_avg))
        if lvl is None:
            return a
        if isinstance(lvl,str):
            if lvl.lower() in ['used','lvl','level']:
                lvl = self.lvl(mol,relative=False)
            else:
                raise ValueError('lvl argument must be integer list or str in [used,lvl,level]')
        if len(a.shape) == 1:
            return a[lvl]
        return a[None,lvl]

    def elecs(self,mol,proj):
        """ Return all projections of this molecule """
        if proj is None: return [mol]
        return list(self.nc.groups[mol].groups[proj].groups.keys())

    def iter_T(self,mol,proj,El,k_avg=True):
        """ Returns an iterable of different transmissions """
        if proj is None:
            grp = self.nc.groups[El]
            tree = [El]
        else:
            grp = self.nc.groups[mol].groups[proj].groups[El]
            tree = [mol,proj,El]
        for var in grp.variables.keys():
            if var.endswith('.T') or var.endswith('.R'):
                yield var,self._get_Data(var,tree=tree,k_avg=k_avg)

    def ADOS(self,mol,proj,El,k_avg=True):
        """ Returns the spectral function DOS from the molecule, projection and electrode """
        return self._get_Data('ADOS',tree=[mol,proj,El],k_avg=k_avg)

# The user requests a specific set of atoms/orbitals
import re
# Compile the reg-exp we search for in the input
re_a   = re.compile('[,]?([0-9-]+\[[0-9,-]+\]|[0-9-]+)[,]?')
re_rng = re.compile('[,]?([0-9-]+)[,]?')

def parse_rng(s):
    """ Parses a string into a list of ranges 
    This range can be formatted like this:
      1,2,3-6
    in which case it will return:
      [1,2,3,4,5,6]
    """
    rng = []
    for la in re_rng.findall(s):
        # We accumulate a list of integers
        tmp = la.split('-')
        if len(tmp) > 2: 
            print('Error in parsing: "'+s+'".')
        elif len(tmp) > 1:
            bo, eo = tuple(map(int,tmp))
            rng.extend(range(bo,eo+1))
        else:
            rng.append(int(tmp[0]))

    return rng

def main():

    # We create all permutations of the available information
    # in the file
    import argparse as opt

    p = opt.ArgumentParser('Generate data files from TBT[|.Proj].nc',
                           formatter_class=opt.ArgumentDefaultsHelpFormatter)

    p.add_argument('-p','--projection',default=None,type=str,
                   help='Only process the molecule[.projection] provided')

    p.add_argument('-k','--k-index',default=None,type=int,
                   help='Only save the equivalent k-index')

    p.add_argument('--eigs',default=False,action='store_true',
                   help='Print out eigs of used levels for the projection')

    p.add_argument('-a','--atom',default=None,action='append',type=str,
                   help='Only save DOS for designated atoms[orbs], could be a list/range (several calls allowed) [1|1,2|1-3|1[2]|1-3[2]]')
    
    p.add_argument('--prefix',default='data',type=str,
                   help='The prefix for the data files when saving the data.')

    p.add_argument('--fmt',default='13.6e',
                   help='The output format for data precision')

    p.add_argument('netcdf', 
                   help='Output NetCDF file from TBtrans to handle')
    
    args = p.parse_args()

    if args.projection:
        Tf = TBTProjFile(args.netcdf)
        process = process_tbt_proj
    else:
        Tf = TBTFile(args.netcdf)
        process = process_tbt

    # Generate file-object
    print('Datafile '+args.netcdf+' created on: '+Tf.date())

    # Generate list of specified atoms
    orbs  = []
    atoms = []
    # we do not allow atoms out-side of the device region
    dev_a = Tf.a_dev
    if args.atom:

        # Re-create a single line of comma separated
        # requests
        args.atom = ','.join(args.atom)

        for lla in re_a.findall(args.atom):
            # We accumulate a list of atoms and orbitals
            
            # First we split in []
            if '[' in lla:
                t = lla.split('[')
                a = parse_rng(t[0])
                # parse the orbitals
                o = parse_rng(t[1].split(']')[0])
            else:
                a = parse_rng(lla)
                o = []
            
            # Update a to only contain atoms in the device region
            a = [ia for ia in a if ia in dev_a]
            # in case the atom list is now empty
            if not a: continue

            # Create list of atoms
            atoms.append(a)

            if o:
                # we have only a subset of the orbitals
                for ia in a:
                    ob = Tf.a2o([ia])
                    # we need to correct the python indices
                    orbs.append(ob[np.array(o,np.int)-1].tolist())
            else:
                orbs.append(Tf.a2o(a).tolist())

    # We also remove dublicates by passing it through a set
    atoms = list(set(itertools.chain.from_iterable(atoms)))
    orbs  = list(set(itertools.chain.from_iterable(orbs)))
    if orbs:
        # Sort the atoms (for clarity), this has no effect otherwise
        atoms.sort()
        # Print out which atoms we project on
        print('Will project DOS onto atoms: ',atoms)
        del atoms
        orbs.sort()
        orbs = Tf.o2idx(orbs)
    else:
        print('Will project DOS onto all device atoms.')
        orbs = np.arange(Tf.no)

    k_idx = args.k_index
    if isinstance(k_idx,int):
        # it must be an integer
        k_idx = np.array([k_idx],np.int) - 1
        if np.any(k_idx < 0) or np.any(len(Tf.wk) <= k_idx):
            raise ValueError('Specified k-index does not exist. Please correct.')
        args.suffix = 'TRANS'
        args.kpt = Tf.kpt[k_idx,:].flatten()
        print('Selected k-point [b]: [{0:.6f}, {1:.6f}, {2:.6f}]'.format(*args.kpt[:]))

    else:
        k_idx = True
        args.suffix = 'AVTRANS'
        args.kpt = None

    process(args,Tf,k_idx,orbs)

def process_tbt_proj(args,Tf,k_idx,orbs):
    """ Processes the TBT.Proj.nc file """

    # decipher molecule and projection
    mol = args.projection.split('.')[0]
    try:
        projs = [args.projection.split('.')[1]]
    except:
        projs = Tf.projs(mol)
    if len(projs) == 0: projs = [None]

    # We start by printing out the eigenvalues if requested
    if args.eigs:
        eig = Tf.eig(mol,lvl='used') * Tf.Ry
        print('Gamma eigenvalues [eV]:')
        print(eig)
        del eig

    # Get the energies (they are in Ry)
    E = Tf.E * Tf.Ry

    # Normalize to number of orbitals in sub-space (DOS is in 1/Ry)
    fac_DOS = 1. / len(orbs) / Tf.Ry
    
    for proj in projs:

        # Loop over all projections associated with this molecule
        for El in Tf.elecs(mol,proj):

            # LHS name
            LHS = El
            if proj: LHS += '.' + mol + '.' + proj

            # Read in DOS
            try:
                # Get ADOS and sum on designated orbitals
                ADOS = np.sum(Tf.ADOS(mol,proj,El,k_avg=k_idx)[:,orbs],axis=-1) * fac_DOS
            except: 
                ADOS = None

            fname = args.prefix+'.TBT.DOS.'+LHS
            save_txt(fname,E,ADOS=ADOS, fmt = args.fmt, kpt=args.kpt)

            for RHS,T in Tf.iter_T(mol,proj,El,k_avg=k_idx):
                # Save k-averaged data (remove .[RT] from variable)
                RHS = RHS[:-2]
                # As this will produce a lot of file names
                # we help the user by printing out the filename:
                fname = args.prefix+'.TBT.'+LHS+'_'+RHS+'.'+args.suffix
                print('Saving transmission data in: '+fname)
                save_txt(fname,E,T=T,ADOS=ADOS,fmt = args.fmt, kpt=args.kpt)

def process_tbt(args,Tf,k_idx,orbs):
    """ Processes the TBT.nc file """

    # Normalize to number of orbitals in sub-space
    fac_DOS = 1. / len(orbs) / Tf.Ry

    # Grab different electrodes in this file
    elecs = Tf.elecs

    # Get the energies (they are in Ry)
    E = Tf.E * Tf.Ry
    try:
        # Get DOS and sum on designated orbitals
        DOS = np.sum(Tf.DOS(k_avg=k_idx)[:,orbs],axis=-1) * fac_DOS
    except: 
        DOS = None

    # Save the DOS
    fname = args.prefix+'.TBT.DOS'
    save_txt(fname,E,DOS=DOS, fmt = args.fmt, kpt=args.kpt)

    for el1 in elecs:
        if not Tf.has_Elec(el1): continue
        try:
            # Get ADOS and sum on designated orbitals
            ADOS = np.sum(Tf.ADOS(el1,k_avg=k_idx)[:,orbs],axis=-1) * fac_DOS
        except: 
            ADOS = None

        fname = args.prefix+'.TBT.DOS.'+el1
        save_txt(fname,E,ADOS=ADOS, fmt = args.fmt, kpt=args.kpt)

        for el2 in elecs:
            try:
                T = Tf.T(el1,el2,k_avg=k_idx)
            except: continue

            # Save k-averaged data
            fname = args.prefix+'.TBT.'+el1+'_'+el2+'.'+args.suffix
            save_txt(fname,E,T,DOS,ADOS, fmt = args.fmt, kpt=args.kpt)


# General text file for saving averaged trans files
def save_txt(file,E,T=None,DOS=None,ADOS=None,fmt='13.5e',kpt=None):
    # Grab length between columns
    cl = int(fmt.split('.')[0])
    # create data stack
    l = [E]
    header  = 'E [eV]'.rjust(cl-2)
    if not T is None:
        l.append(T)
        header += ' ' + 'T [G0]'.rjust(cl)
    if not DOS is None: 
        l.append(DOS)
        header += ' ' + 'DOS'.rjust(cl)
    if not ADOS is None: 
        l.append(ADOS)
        header += ' ' + 'Spec.DOS'.rjust(cl)
    if len(l) == 1: return # this is only the energy grid
    # Append date of file-creation to the file.
    header += '\nDate: '+_TODAY
    if kpt is None:
        header += '\nk-averaged'
    else:
        header += '\nk-point [b]: [{0:.6f}, {1:.6f}, {2:.6f}]'.format(*kpt[:])
    data = np.vstack(l)
    np.savetxt(file,data.T,header=header,fmt='%'+fmt)
    
if __name__ == '__main__':
    main()


