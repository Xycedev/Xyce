
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
endif

# conditionally link in radiation-aware device models.
if RADMODELS
    RADLD = $(top_builddir)/src/DeviceModelPKG/SandiaModels/libSandiaModels.la
else
    RADLD = 
endif

# conditionally link in radiation-aware device models.
if NONFREEMODELS
    NONFREELD = $(top_builddir)/src/DeviceModelPKG/Xyce_NonFree/libNonFree.la
else
    NONFREELD = 
endif

# standalone testHarness1685 executable
check_PROGRAMS = testHarness1685
 
BUG1685TESTSOURCES = \
  $(srcdir)/testHarness1685.C 

# standalone testHarness1685  executable
testHarness1685_SOURCES = $(BUG1685TESTSOURCES)
testHarness1685_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD)
testHarness1685_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
