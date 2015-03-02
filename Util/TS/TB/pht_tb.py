#!/usr/bin/env python

from __future__ import print_function
# Should make it "backwards" compatible down to 2.6

from tbt_tb import PeriodicTable, SIESTA_UNITS
from tbt_tb import TBT_Geom, TBT_Model, TBT_dH, OutputFile

# This utility function creates a phonon transport input
# to phtrans.

import numpy as np

# It is _exactly_ the same as tbt_tb with some small modifications
# to allow some differentations.

class PHT_Geom(TBT_Geom):
    """
    Geometry object for phonon geometries
    Very similar to the TBT_Geom, however with the change
    that all atoms default to 3 orbitals (x/y/z)

    See ``TBT_Geom`` for specifications.
    """
    def __init__(self,cell,xa,dR=None,n_orb=3,Z=1,update_sc=False):
        super(PHT_Geom, self).__init__(cell,xa,dR=dR,n_orb=n_orb,Z=Z,update_sc=update_sc)

class GULP(OutputFile):
    """ Class that extracts information from an GULP output """
    def __init__(self,filename):
        super(GULP,self).__init__(filename)

    def read_model(self,na_u,dtype=np.float):
        """ 
        Returns a new dynamical matrix in the following lines
        Always returns a coo matrix.
        """
        if not hasattr(self,'fh'):
            # The file-handle has not been opened
            with self:
                return self.read_model(na_u,dtype=dtype)

        from scipy.sparse import lil_matrix
        
        no = na_u * 3
        dyn = lil_matrix( (no,no) ,dtype = dtype)

        self.step_to('Real Dynamical matrix')

        # skip 1 line
        self.fh.readline()
        
        l = self.fh.readline()
        i = 0 ; j = 0
        while len(l.strip()) > 0:
            ls = [float(f) for f in l.split()]
            # add the values (12 values == 3*4)
            for k in range(4):
                dyn[i,j  ] = ls[k*3  ]
                dyn[i,j+1] = ls[k*3+1]
                dyn[i,j+2] = ls[k*3+2]
                j += 3
                if j >= no: 
                    i += 1
                    j = 0
            l = self.fh.readline()

        if dtype in [np.complex,np.complex64,np.complex128]:

            # Find the imaginary part
            self.step_to('Imaginary Dynamical matrix')

            # skip 1 line
            self.fh.readline()
            
            l = self.fh.readline()
            i = 0 ; j = 0
            while len(l.strip()) > 0:
                ls = np.array([float(f) for f in l.split()])
                # add the values
                for k in range(4):
                    dyn[i,j  ] += 1j*ls[k*3  ]
                    dyn[i,j+1] += 1j*ls[k*3+1]
                    dyn[i,j+2] += 1j*ls[k*3+2]
                    j += 3
                    if j >= no: 
                        i += 1
                        j = 0
                l = self.fh.readline()

        # Convert the GULP data to standard units
        dyn = dyn.tocoo()
        dyn.data[:] *= ( 521.469 * 1.23981e-4 ) ** 2

        return dyn

    def read_models(self,na_u):
        """
        Returns a list of dynamical matrices with associated k-points.
        """
        if not hasattr(self,'fh'):
            # The file-handle has not been opened
            with self:
                return self.read_models(na_u)

        # skip to region with dynamical matrices
        self.step_to('Phonon Calculation')

        # read dynamical matrices
        dyns = []
        found, qline = self.step_to('K point')
        while found:
            # Parse q-point
            qs = qline.split('=')[1].split()[:3]
            qpt = np.array([float(q) for q in qs])
            if np.allclose(np.zeros((3,)),qpt,rtol=1.e-6):
                dtype = np.float
            else:
                dtype = np.complex
            dyns.append([qpt, self.read_model(na_u,dtype=dtype)])

            found, qline = self.step_to('K point')

        return dyns

    def read_geom(self):
        """
        Returns a ``PHT_Geom`` by reading a the assigned output file.
        """
        if not hasattr(self,'fh'):
            # The file-handle has not been opened
            with self:
                return self.read_geom()

        Z = None
        xa = None
        cell = None

        self.step_to('Cartesian lattice vectors')

        # skip 1 line
        self.fh.readline()
        cell = np.zeros((3,3),np.float)
        for i in range(3):
            l = self.fh.readline().split()
            cell[i,0] = float(l[0])
            cell[i,1] = float(l[1])
            cell[i,2] = float(l[2])
            
        # Skip to fractional coordinates
        self.step_to('Final fractional coordinates')
                    
        # We skip 5 lines
        for i in range(5): self.fh.readline()
        
        l = 'STEP'
        Z = []
        fxyz = []
        l = self.fh.readline()
        while l[0] != '-':
            ls = l.split()
            Z.append(ls[1])
            fxyz.append([float(f) for f in ls[3:6]])
            l = self.fh.readline()
        fxyz = np.array(fxyz,np.float)
        fxyz.shape = (-1,3)
        # Correct the fractional coordinates
        xa = np.empty(fxyz.shape,np.float)
        xa[:,0] = fxyz[:,0] * np.sum(cell[:,0])
        xa[:,1] = fxyz[:,1] * np.sum(cell[:,1])
        xa[:,2] = fxyz[:,2] * np.sum(cell[:,2])
        
        if Z is None or cell is None or xa is None:
            raise ValueError('Could not read in cell information and/or coordinates')
        
        return PHT_Geom(cell,xa,Z=Z)
    

class PHT_Model(TBT_Model):
    """
    Inherits TBT_Model to create a phonon tight-binding model.

    See ``TBT_Model`` for specifications.
    """

    # TODO, do this in another way, with this object ALL energies
    # are transferred to Ry ** 2
    Ry = SIESTA_UNITS.Ry ** 2

    @staticmethod
    def read_output_periodic(output,V=[],geom=None,model=None):
        if model is None: model = PHT_Model
        return TBT_Model.read_output_periodic(output,V=V,geom=geom,model=model)

    @staticmethod
    def read_output(output,geom=None,model=None):
        if model is None: model = PHT_Model
        return TBT_Model.read_output(output,geom=geom,model=model)

    def correct_Newton(self):
        """
        Sometimes the dynamical matrix does not obey Newtons laws.

        We correct the dynamical matrix by imposing zero force.

        Correcting for Newton forces the matrix to be finalized.
        """
        from scipy.sparse import lil_matrix

        # We need to ensure that the block diagonal exists
        # to create the Newton-corrections
        # This will ONLY work if the matrix has not been 
        # finalized, OR if the block diagonal already exists
        for ja in range(self.na_u):
            jo = ja * 3
            for j in range(3):
                for i in range(3):
                    H, S = self[jo+j,jo+i]
                    if j == i:
                        self[jo+j,jo+i] = H, 1.
                    else:
                        self[jo+j,jo+i] = H, S

        # Create UC dynamical matrix
        d_sc, S_sc = self.tocsr() ; del S_sc
        d_sc = d_sc.tocoo()
        d_uc = lil_matrix((self.no_u,self.no_u),dtype=np.float)

        # Convert SC to UC
        for j, i, d in zip(d_sc.row,d_sc.col,d_sc.data):
            d_uc[j, i % self.no_u] += d
        del d_sc
        d_uc = d_uc.tocsr()

        # we need to correct the dynamical matrix found in GULP
        # This ensures that Newtons laws are obeyed, (i.e. 
        # action == re-action)
        ptbl = PeriodicTable()
        om = ptbl.atomic_mass(self.geom.Z)
        del ptbl

        for ja in range(self.na_u):

            # Create conversion to force-constant, and revert back
            # after correcting
            am = om[ja]
            MM = np.sqrt(am * om)
            jo = ja * 3

            for j in range(3):
                for i in range(3):
                    H, S = self[jo+j,jo+i]
                    self[jo+j,jo+i] = H - d_uc[jo+j,i::3].multiply(MM).sum() / am, S

        del d_uc
        
