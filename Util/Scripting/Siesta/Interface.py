"""
Extend the ASE classes to provide Siesta-specific functionality
"""

import ASE.Atom
import ASE.ListOfAtoms

#-------------------
class Atom(ASE.Atom):

    """Atom object for Siesta. Includes label"""

    def __init__(self, symbol=None, position=(0, 0, 0),
                 Z=None, mass=None, tag=0,
                 momentum=None, velocity=None,
                 magmom=0.0, label=None, valence_gs=None):
        """Atom(symbol, position, ...) -> atom object."""
	ASE.Atom.__init__(self, symbol, position,
                 Z, mass, tag, momentum, velocity,
                 magmom)

        if label is None: 
          self.label = self.symbol
        else:
          self.label = label

        if valence_gs is None: 
          self.valence_gs = []
        else:
          self.valence_gs = valence_gs

    def GetLabel(self):
        """Get label."""
        return self.label

    def SetLabel(self,label):
        """Set label."""
        self.label = label

#---------------------------------
class Crystal(ASE.ListOfAtoms):
       def dummy(self):
          print "Just checking in"


if __name__=="__main__":
  a = Atom("H",label="H_surf")
  b = Atom("O",label="O_bulk")
  c = Atom(Z=0,label="ON-0.50000",valence_gs=[[0,2],[1,3.5]])
  cryst = Crystal([a,b,c])
  print cryst
  cryst.dummy()

  import siesta

  pepe = siesta.Siesta()
  pepe.run(cryst)







