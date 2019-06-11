import os, shutil, sys
from linecache import getline
import linecache
from math import sqrt
from math import cos
from math import sin
from math import pow
from math import radians
import numpy as np



NAME_OF_INPUT_FILE="06-cubic.36-spg230.inp"


vector_info_temp=[]
vector_info=[]
with open(NAME_OF_INPUT_FILE) as f:
    for ind, line in enumerate(f,2):
        if line.rstrip()=="%LatticeVectors":
            for i in range(3):
                #print(getline(f.name, ind+i))
                vector_info_temp.append(getline(f.name, ind+i))
for i in range(3):
    vector_info.append(vector_info_temp[i].split())
lattice=np.array([(float(vector_info[0][0]),float(vector_info[0][2]),float(vector_info[0][4])),
                  (float(vector_info[1][0]),float(vector_info[1][2]),float(vector_info[1][4])),
                  (float(vector_info[2][0]),float(vector_info[2][2]),float(vector_info[2][4]))])                    


vector_parameter_temp=[]
vector_parameter=[]
with open(NAME_OF_INPUT_FILE) as f:
    for ind, line in enumerate(f,2):
        if line.rstrip()=="%LatticeParameters":
            for i in range(3):
                #print(getline(f.name, ind+i))
                vector_parameter_temp.append(getline(f.name, ind+i))
for i in range(3):
    vector_parameter.append(vector_parameter_temp[i].split())
a=float(vector_parameter[0][0])
b=float(vector_parameter[0][2])
c=float(vector_parameter[0][4])    


siesta_lattice=np.array([(lattice[0][0]*a,lattice[0][1]*b,lattice[0][2]*c),
                         (lattice[1][0]*a,lattice[1][1]*b,lattice[1][2]*c),
                         (lattice[2][0]*a,lattice[2][1]*b,(lattice[2][2]*c))])


import itertools
with open(NAME_OF_INPUT_FILE) as f,open('result.txt', 'w') as fout:
    it = itertools.dropwhile(lambda line: line.strip() != '%ReducedCoordinates', f)
    it = itertools.islice(it, 1, None)
    it = itertools.takewhile(lambda line: line.strip() != '%', it)
    fout.writelines(it)



atomic_info_temp=[]
with open ('result.txt') as la:
    while True:
        line = la.readline().strip()
        atomic_info_temp.append(line)
        print (line)
        if line == '':
            break


# In[73]:

atomic_info=[]
i_num=range(len(atomic_info_temp)-1)# range(96)
for i in i_num:
    atomic_info.append(atomic_info_temp[i].split())


# In[87]:

for i in i_num:
    print (str(atomic_info[i][2])+"   "+str(atomic_info[i][4])+"   "+str(atomic_info[i][6])+"  1  "+str(i+1)+"   H" )


# In[99]:

filename1 = str(NAME_OF_INPUT_FILE)+".fdf"
f1=open(filename1,'w')
f1.write("#--------------------------------------- \n")
f1.write("#Input Generated for Space Group Testing \n")
f1.write("#--------------------------------------- \n")
f1.write("SystemLabel  "+str(NAME_OF_INPUT_FILE) +" \n")
f1.write("NumberOfSpecies 1 \n")
f1.write("%block ChemicalSpeciesLabel \n")
f1.write(" 1  1  H \n")
f1.write("%endblock ChemicalSpeciesLabel \n")
f1.write("LatticeConstant  1.0 Ang   \n")
f1.write("%block LatticeVectors \n")
for j in range(3):
    f1.write("  "+str(siesta_lattice[j][0])+"   "+str(siesta_lattice[j][j])+"   "+str(siesta_lattice[j][2])+"  \n")
f1.write("%endblock LatticeVectors \n")
f1.write("%block AtomicCoordinatesAndAtomicSpecies \n")
for i in i_num:
    f1.write (str(atomic_info[i][2])+"   "+str(atomic_info[i][4])+"   "+str(atomic_info[i][6])+"  1  "+str(i+1)+"   H \n" )
f1.write("%endblock AtomicCoordinatesAndAtomicSpecies \n")
f1.write("MaxSCFIterations 	0  \n")
f1.write("MD.TypeOfRun    CG   \n")
f1.write("MD.NumCGsteps   0  \n")
f1.close() 
print ("Done!")



