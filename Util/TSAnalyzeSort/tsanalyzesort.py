#!/usr/bin/env python

# A small script for reading in the pivoting determined by the
# Cuthill-Mckee algorithm for optimizing the BTD block-sizes

# This code has been fully implemented by Nick Papior Andersen, 2013
# Contact the author prior to reconstructing the code at: nickpapior<at>gmail.com

# It will only print out the structure, so you have to pipe it to a
# file and copy paste the coordinates into the FDF file.

import sys
try:
    import argparse as opt
except:
    print("We use the argparse module for options, please update "+
          "your python installation to use this package.")
    sys.exit(1)

def run():
    """ The main program """
    
    parser = opt.ArgumentParser("Prints out a reorganized FDF coordinate block using the optimization obtained by TS.Analyze. Print to stdout.")
    
    parser.add_argument("-s","--struct",dest='FDF',type=str,required=True,
                        help="The FDF file which contains the AtomicCoordinatesAndAtomicSpecies block, does not follow %%include")
    parser.add_argument("-A","--analyze",dest='out',type=str,required=True,
                        help="The output of the transiesta run with TS.Analyze T")
    parser.add_argument("-b","--block",action="store_true",dest="add_block",
                        help="Add the block designation as well")

    # Parse the arguments of the options
    args = parser.parse_args()

    # Read in structure lines
    struct_lines = read_AtomCoordinates(args.FDF)
    # Read the pivoting of the atoms
    pvt = read_Pivot(args.out,len(struct_lines))
    # Check that the pivoting has the correct form
    check_Pivot(pvt)
    # Print to std out the pivoted structure
    if args.add_block:
        print("%block AtomicCoordinatesAndAtomicSpecies")
    print_new_AtomCoordinats(struct_lines,pvt)
    if args.add_block:
        print("%endblock AtomicCoordinatesAndAtomicSpecies")


_COORD_BLOCK = "AtomicCoordinatesAndAtomicSpecies".lower()
def read_AtomCoordinates(f):
    """
    Reads in the coordinates from the FDF file supplied
    
    It returns in a list the strings for each line
    with only the atoms.

    It cannot handle full lines of comments in the
    block.
    """
    found = False
    atms = []
    fh = open(f,'r')
    for line in fh.readlines():
        if found: 
            # if we are done reading the blocks then quit
            if _COORD_BLOCK in line.lower(): break
            if "endblock" in line.lower(): break
            atms.append(line.replace('\n','').replace('\r',''))
        if _COORD_BLOCK in line.lower():
            found = True
    fh.close()
    return atms

_OUT_PLACE = "transiesta: Central region..."
def read_Pivot(f,na):
    """
    Reads in the pivot table out-put by transiesta
    """
    found = False
    # Create the default pivoting array
    pvt = range(na)
    fh = open(f,'r')
    for line in fh.readlines():
        if found: 
            # If -> is in the line we know we have something
            if '->' in line:
                old,new = line.split('->')
                try:
                    old = int(old) - 1
                    new = int(new) - 1
                    pvt[new] = old
                except: pass
        if _OUT_PLACE in line:
            found = True
    fh.close()
    return pvt

def check_Pivot(pvt):
    """
    We check the pivoting indices to ensure that no number occurs
    more than once
    """
    ok = True
    for ia in range(len(pvt)):
        c = pvt.count(ia)
        if c != 1:
            ok = False
            if c == 0:
                print("Atom {0} has been removed due to the pivot table!".format(ia))
            else:
                print("Atom {0} is pivoted {1} times!".format(ia,c))
    if not ok:
        print("Please contact Nick Papior Andersen at nickpapior<at>gmail.com" +
              " with your FDF files for reviewing. Seems like a bug in the algorithm.")
        sys.exit(1)

def print_new_AtomCoordinats(struct,pvt):
    """
    Prints the new sorted structure
    """
    for i in range(len(pvt)):
        print(struct[pvt[i]])

if __name__ == "__main__":
    run()
    # Do a clean exit
    sys.exit(0)
    
