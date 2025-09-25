.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
include ../../make.inc
#
module:  libwxml.a
#
#
OBJFILES= m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o \
          m_wxml_text.o m_wxml_core.o xmlf90_wxml.o
#
libwxml.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libwxml.a $(OBJFILES)
	ranlib libwxml.a
	#-$(RANLIB) libwxml.a
#
clean:
	rm -f *.o *.*d *.a

# DO NOT DELETE THIS LINE - used by make depend
xmlf90_wxml.o: m_wxml_core.o m_wxml_text.o
m_wxml_core.o: m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o

#
.f90.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f90
.F90.o:
	$(F90) $(FFLAGS) $(DEFS) $(INC) -c $*.F90
.f.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f
.F.o:
	$(F90) $(FFLAGS) $(DEFS) $(INC) -c $*.F
