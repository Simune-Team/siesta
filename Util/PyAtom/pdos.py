from xml.sax import saxutils
from Numeric import *


#def normalize_whitespace(text):
#  "Remove redundant whitespace from a string"
#   return string.join(string.split(text), ' ')

class PDOS(saxutils.DefaultHandler):
   def __init__(self):
	self.inData = 0
	self.inOrbital = 0
	self.data = ""

   def startElement(self, name, attrs):
      if name == 'orbital':
         self.inOrbital = 1
	 n = attrs.get('n', None)
         l = attrs.get('l', None)
         m = attrs.get('m', None)
         z = attrs.get('z', None)
         print "Orbital", n, l, m, z
      if name == 'data':
         self.inData = 1

   def characters(self, ch):
         if self.inData:
            self.data = self.data + ch

   def endElement(self, name):
       if name == 'data':
         self.inData = 0
#
#        Process data ....
#         print self.data[0:30]

#------------------------------
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces

if __name__ == '__main__':
        # Create a parser
        parser = make_parser()
        # Tell the parser we are not interested in XML namespaces
        parser.setFeature(feature_namespaces, 0)

        # Create the handler
        dh = PDOS()

        # Tell the parser to use our handler
        parser.setContentHandler(dh)

        # Parse the input
        parser.parse('pdos.xml')