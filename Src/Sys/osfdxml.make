# $Id: osfdxml.make,v 1.3 1999/01/31 12:05:54 emilio Exp $
#
FC=f77
FFLAGS= -fast -tune host
FFLAGS_DEBUG= -g3 -fast -tune host
LIBS= -ldxml
SYS=osfdxml
RANLIB=@echo No ranlib required for 
#
