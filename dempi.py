#!/usr/bin/env python
#
# De-MPI source files
# Alberto Garcia, Dec 14, 2000
#
# The algorithm keeps track of
#
#  1. Whether we are within an MPI preprocessor block
#     (Since there are other kinds of preprocessor blocks..)
#  2. Whether we are in the serial or in the parallel part
# 
#   In addition, we need to keep track of any other #ifdef...#endif
#   blocks that might be opening *within* an MPI block. For that
#   we use a variable "depth".
#
#   This is a first version. It could be made more compact.

import sys
import string
import exceptions

class DeMPIError(exceptions.Exception):
    def __init__(self,errno,msg):
        self.errno=errno
        self.errmsg=msg

filename=sys.argv[1]

f=open(filename,"r")

# State flags
serial = 1                    
mpi_block = 0

line=f.readline() 
while line:

  t=string.split(line)
##  print t

  if len(t) == 0:
      if serial: print line[:-1]

  elif t[0]=="#else" and mpi_block and depth == 0:
       if serial:
           serial = 0
       else:
           serial = 1

  elif t[0]=="#endif":
       if mpi_block:
         if depth == 0:    # The MPI block is ending...
           if not serial: serial=1
           mpi_block = 0
         else:
           if serial: print line[:-1]
           depth = depth - 1
       else:
         if serial: print line[:-1]

  elif len(t) == 1:
      if serial: print line[:-1]

  elif t[0]=="#ifdef":
      if t[1] == "MPI":
        if mpi_block: raise DeMPIError(1,'Nested MPI block!!')
        mpi_block = 1
        depth = 0
        serial = 0
      else:
        if serial: print line[:-1]
        if mpi_block:
          depth = depth + 1

  elif t[0]=="#ifndef":
       if t[1]=="MPI":
         mpi_block = 1
         serial = 1
       else:
         if serial: print line[:-1]
         if mpi_block :
           depth = depth + 1

  else:
         if serial: print line[:-1]


  line = f.readline()

f.close()




