
AM_CPPFLAGS = @Xyce_INCS@ -DXYCEBUILDTIMESTAMP=\"`date +%Y%m%d%H%M`\"

noinst_LTLIBRARIES = libUtility.la

if BUILD_INTEL_FFT_SUPPORT
intelfft_SOURCES = \
  $(srcdir)/src/N_UTL_IntelFFT_Interface.C \
  $(srcdir)/include/N_UTL_IntelFFT_Interface.hpp
else
intelfft_SOURCES =
endif

if BUILD_FFTW_SUPPORT
fftw_SOURCES = \
  $(srcdir)/src/N_UTL_FFTW_Interface.C \
  $(srcdir)/include/N_UTL_FFTW_Interface.hpp
else
fftw_SOURCES =
endif

if THREADS 
threads_LIB = -lpthread 
threads_CFLAGS = -pthread 
else 
threads_LIB = 
threads_CFLAGS = 
endif 

libUtility_la_SOURCES = \
  $(intelfft_SOURCES) \
  $(fftw_SOURCES) \
  $(srcdir)/src/N_UTL_BreakPoint.C \
  $(srcdir)/src/N_UTL_Demangle.C \
  $(srcdir)/src/N_UTL_NoCase.C \
  $(srcdir)/src/N_UTL_Expression.C \
  $(srcdir)/src/N_UTL_ExpressionData.C \
  $(srcdir)/src/N_UTL_ExpressionInternals.C \
  $(srcdir)/src/N_UTL_APFT.C \
  $(srcdir)/src/N_UTL_FFTInterface.C \
  $(srcdir)/src/N_UTL_LogStream.C \
  $(srcdir)/src/N_UTL_ReportHandler.C \
  $(srcdir)/src/N_UTL_Marshal.C \
  $(srcdir)/src/N_UTL_JSON.C \
  $(srcdir)/src/N_UTL_SendCURL.C \
  $(srcdir)/src/N_UTL_ByThreadStreambuf.C \
  $(srcdir)/src/N_UTL_Misc.C \
  $(srcdir)/src/N_UTL_NetlistLocation.C \
  $(srcdir)/src/N_UTL_OptionBlock.C \
  $(srcdir)/src/N_UTL_Op.C \
  $(srcdir)/src/N_UTL_OpBuilder.C \
  $(srcdir)/src/N_UTL_Param.C \
  $(srcdir)/src/N_UTL_RandomNumbers.C \
  $(srcdir)/src/N_UTL_Registry.C \
  $(srcdir)/src/N_UTL_Timer.C \
  $(srcdir)/src/N_UTL_CPUTime.C \
  $(srcdir)/src/N_UTL_MallocUsed.C \
  $(srcdir)/src/N_UTL_FormatMemorySize.C \
  $(srcdir)/src/N_UTL_FormatTime.C \
  $(srcdir)/src/N_UTL_StatMetricTraits.C \
  $(srcdir)/src/N_UTL_Stats.C \
  $(srcdir)/src/N_UTL_Platform.C \
  $(srcdir)/src/N_UTL_PrintTable.C \
  $(srcdir)/src/N_UTL_PrintStats.C \
  $(srcdir)/src/N_UTL_WallTime.C \
  $(srcdir)/src/N_UTL_Version.C \
  $(srcdir)/src/N_UTL_Algorithm.C \
  $(srcdir)/src/N_UTL_NameLevelKey.C \
  $(srcdir)/include/N_UTL_TypeIndex.h \
  $(srcdir)/include/N_UTL_NameLevelKey.h \
  $(srcdir)/include/N_UTL_DeleteList.h \
  $(srcdir)/include/N_UTL_Interface_Enum_Types.h \
  $(srcdir)/include/N_UTL_Expression.h \
  $(srcdir)/include/N_UTL_ExpressionInternals.h \
  $(srcdir)/include/N_UTL_ExpressionData.h \
  $(srcdir)/include/N_UTL_FFTInterface.hpp \
  $(srcdir)/include/N_UTL_RandomNumbers.h \
  $(srcdir)/include/N_UTL_Registry.h \
  $(srcdir)/include/N_UTL_Demangle.h \
  $(srcdir)/include/N_UTL_IndentStreamBuf.h \
  $(srcdir)/include/N_UTL_TeeStreamBuf.h \
  $(srcdir)/include/N_UTL_PThread.h \
  $(srcdir)/include/N_UTL_ByThreadStreambuf.h \
  $(srcdir)/include/N_UTL_LogStream.h \
  $(srcdir)/include/N_UTL_ReportHandler.h \
  $(srcdir)/include/N_UTL_Marshal.h \
  $(srcdir)/include/N_UTL_JSON.h \
  $(srcdir)/include/N_UTL_SendCURL.h \
  $(srcdir)/include/N_UTL_Op.h \
  $(srcdir)/include/N_UTL_OpBuilder.h \
  $(srcdir)/include/N_UTL_APFT.h \
  $(srcdir)/include/N_UTL_DFTInterfaceDecl.hpp \
  $(srcdir)/include/N_UTL_FFTInterfaceDecl.hpp \
  $(srcdir)/include/N_UTL_FixedQueue.h \
  $(srcdir)/include/N_UTL_Functors.h \
  $(srcdir)/include/N_UTL_MachDepParams.h \
  $(srcdir)/include/N_UTL_Listener.h \
  $(srcdir)/include/N_UTL_Dump.h \
  $(srcdir)/include/N_UTL_Misc.h \
  $(srcdir)/include/N_UTL_NetlistLocation.h \
  $(srcdir)/include/N_UTL_SaveIOSState.h \
  $(srcdir)/include/N_UTL_BreakPoint.h \
  $(srcdir)/include/N_UTL_Factory.h \
  $(srcdir)/include/N_UTL_OptionBlock.h \
  $(srcdir)/include/N_UTL_Packable.h \
  $(srcdir)/include/N_UTL_Param.h \
  $(srcdir)/include/N_UTL_Timer.h \
  $(srcdir)/include/N_UTL_CPUTime.h \
  $(srcdir)/include/N_UTL_FormatMemorySize.h \
  $(srcdir)/include/N_UTL_FormatTime.h \
  $(srcdir)/include/N_UTL_MallocUsed.h \
  $(srcdir)/include/N_UTL_StatMetricTraits.h \
  $(srcdir)/include/N_UTL_Stats.h \
  $(srcdir)/include/N_UTL_Platform.h \
  $(srcdir)/include/N_UTL_PrintStats.h \
  $(srcdir)/include/N_UTL_PrintTable.h \
  $(srcdir)/include/N_UTL_WallTime.h \
  $(srcdir)/include/N_UTL_Graph.h \
  $(srcdir)/include/N_UTL_Xyce.h \
  $(srcdir)/include/N_UTL_IndentStreamBuf.h \
  $(srcdir)/include/N_UTL_NoCase.h \
  $(srcdir)/include/N_UTL_Version.h \
  $(srcdir)/include/N_UTL_Algorithm.h

BUILT_SOURCES = timestamp-generator

CLEANFILES = timestamp.h

noinst_PROGRAMS = ThreadTest

ThreadTest_SOURCES = $(srcdir)/unit_test/ThreadTest.C
ThreadTest_CFLAGS = $(threads_CFLAGS)
ThreadTest_LDADD = libUtility.la
ThreadTest_LDFLAGS = $(AM_LDFLAGS) $(threads_LIB) 

timestamp.h: timestamp-generator

timestamp-generator:
	echo "/* TIMESTAMP APPLIED */" > timestamp.h

