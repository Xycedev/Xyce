
# -- process files for this library ------------------------------------------

# class source list
set ( OutputInterface_SOURCES  
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_OutputFileBase.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_OutputPrn.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_OutputHDF5.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureManager.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureBase.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureRiseFallDelay.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureAverage.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureDuty.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureFrequency.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureMax.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureMin.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureOffTime.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureOnTime.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasurePeakToPeak.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureRMS.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureFindWhen.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureEquationEvaluation.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureDerivativeEvaluation.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureIntegralEvaluation.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureRelativeError.C
     ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_MeasureFourier.C
      CACHE INTERNAL "X_S" )

  
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      OutputInterface
      CACHE INTERNAL "X_L" )
