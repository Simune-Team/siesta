FFLAGS =  -fast -tune host
FC = f77  
CF = f77  
CFT = f77 

OBJECTS = atompla.o chkdim.o colinear.o dot.o io.o iodm.o denchar.o\
          dmna.o fdf.o length.o matvect.o\
          neighb.o parse.o paste.o phiatm.o planed.o prmem.o\
          ranger.o readpla.o redata.o rhoofr.o rlylm.o splint.o volcel.o\
          wrout.o

plden.x: $(OBJECTS)
	$(FC) $(FFLAGS) -o plden.x $(LFLAGS)    $(OBJECTS)


denchar.o:  denchar.f   siesta.h
redata.o:   redata.f    atom.h
phiatm.o:   phiatm.f    atom.h
ranger.o:   ranger.f    ranger.h
