//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_OutputMgr.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.219.2.4 $
//
// Revision Date  : $Date: 2014/08/29 21:26:34 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMgr_h
#define Xyce_N_IO_OutputMgr_h

#include <iterator>
#include <list>
#include <set>
#include <string>
#include <vector>

#ifdef Xyce_USE_HDF5
#include <hdf5.h>
#endif

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_TOP_fwd.h>
#include <N_LAS_fwd.h>

#include <N_ANP_SweepParam.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_Objective.h>
#include <N_IO_Outputter.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>

class N_MPDE_Manager;

namespace Xyce {
namespace IO {

namespace OutputType {
enum OutputType {DC, TRAN, AC, AC_IC, HB_FD, HB_TD, HB_IC, HB_STARTUP, DCOP, HOMOTOPY, MPDE, MPDE_IC, SENS};
}


/***
 * The mapping from .PRINT <print-type> to outputters and then activiting the
 * outputters at the appropriate time in the appropriate context is the
 * trick here.
 *
 * The user communicates the output variables wanted via the .PRINT
 * <print-type>.  The .PRINT <print-type> creates a mapping from
 * print-type to output-type where each output-type represents a set of
 * parameters to the written when the print-type is active.  The
 * preparePrintParameters() function creates a mapping from .PRINT
 * <print-type> to outputters, using the output-types to identify the
 * parameters.  The application informs the output manager
 * what needs to be written when called via the ActiveOutput sentry
 * class and its functions.
 *
 * The ActiveOutput sentry class manages the stack of active
 * outputters.  The construction pushes an empty vector of outputters
 * onto the stack and the destructor pops the vector off, restoring the
 * previous vector.
 */


typedef std::map<std::string, Objective> ObjectiveMap;
typedef std::map<PrintType::PrintType, std::vector<Outputter::Interface *> > OutputterMap;
typedef std::vector<std::vector<Outputter::Interface *> > ActiveOutputterStack;
typedef std::map<OutputType::OutputType, std::vector<PrintParameters> > OutputParameterMap;
typedef std::set<Analysis::Analysis_Mode> EnabledAnalysisSet;
typedef std::vector<std::pair<std::string, std::string> > StringPairVector;

enum Domain {DOMAIN_TIME, DOMAIN_FREQUENCY};

//-----------------------------------------------------------------------------
// Class         : OutputMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
class OutputMgr : public Util::Listener<Analysis::StepEvent>
{
  friend void printLineDiagnostics(
    Parallel::Machine                                                   comm,
    const OutputMgr &                                                   output_manager,
    const Measure::Manager &                                            measure_manager,
    const std::set<std::string> &                                       node_names,
    const std::map<std::string, RefCountPtr<Device::InstanceBlock> > &  device_map,
    AliasNodeMap &                                                      alias_node_map);
  friend void deferredPrintLineDiagnostics(
    Parallel::Machine   comm,
    OutputMgr &         output_manager);

public:
  struct OutputterKey
  {
    OutputterKey(int analysis_mode, int output_type, int format)
      : analysisMode_(analysis_mode),
        outputType_(output_type),
        format_(format)
    {}

    const int   analysisMode_;
    const int   outputType_;
    const int   format_;
  };

  typedef Outputter::Interface *(*OutputterFactory)(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);
  typedef std::map<OutputterKey, OutputterFactory> OutputterFactoryMap;
  typedef std::map<std::string, std::pair<int, std::ostream *> > OpenPathStreamMap;

  OutputMgr(CmdParse &command_line, Util::Op::BuilderManager &op_builder_manager);
  ~OutputMgr();

private:
  OutputMgr(const OutputMgr & rhs);
  OutputMgr &operator=(const OutputMgr &rhs);

public:
  // registration functions:
  bool registerTopology(Topo::Topology * topology)
  {
    topology_ = topology;
    return true;
  }

  bool registerDeviceInterface(Device::DeviceInterface * device_interface)
  {
    deviceInterface_ = device_interface;
    return true;
  }

  bool registerAnalysisManager(Analysis::AnalysisManager *analysis_manager);

  bool getOutputIntervals(double & initialInterval, std::vector<std::pair< double, double > > & intervalPairs) const;

  bool registerPkgOptionsMgr( PkgOptionsMgr &pkgOpt);
  bool registerPkgOptionsMgrInitialConditions(PkgOptionsMgr &pkgOpt);

  bool registerMeasureMgr(Measure::Manager *measure_manager)
  {
    measureManager_ = measure_manager;
    return true;
  }

  bool registerFourierMgr(FourierMgr *fourier_manager)
  {
    fourierManager_ = fourier_manager;
    return true;
  }

  bool registerDCOptions(const Util::OptionBlock & option_block);
  bool registerTranOptions(const Util::OptionBlock & option_block);
  bool registerMPDETranOptions(const Util::OptionBlock & option_block);
  bool registerHBOptions(const Util::OptionBlock & option_block);
  bool registerSTEPOptions(const Util::OptionBlock & option_block);
  bool registerOutputOptions(const Util::OptionBlock & option_block);
  bool registerDeviceOptions(const Util::OptionBlock & option_block);

  bool parsePRINTBlock(const Util::OptionBlock & print_block);
  bool registerLoad(const Util::OptionBlock & option_block);
  bool registerSens(const Util::OptionBlock & option_block);
  bool registerSensOptions(const Util::OptionBlock & option_block);

  bool registerDCOPOptions(const Util::OptionBlock & option_block);
  bool registerIC(const Util::OptionBlock & option_block);
  bool registerNodeSet(const Util::OptionBlock & option_block);
  bool registerSave(const Util::OptionBlock & option_block);

  bool setOBJECTIVEParams(const Util::OptionBlock & option_block);

  void registerOutputter(int analysis_mode, int output_type, int format, OutputterFactory factory) {
    outputterFactoryMap_[OutputterKey(analysis_mode, output_type, format)] = factory;
  }

  void notify(const Analysis::StepEvent &step_event);

  void fixupPrintParameters(Parallel::Machine comm, PrintParameters &print_parameters);

public:
  void checkPrintParameters(Parallel::Machine comm);
  void prepareOutput(Parallel::Machine comm, Analysis::Analysis_Mode analysis_mode);
  void setStepSweep(const Analysis::SweepVector &step_sweep_parameters);
  void setDCSweep(const Analysis::SweepVector &dc_sweep_parameters);

  const std::vector<char> &getTypeRefs() const
  {
    return typeRefs_;
  }

  const NodeNamePairMap &getAllNodes() const
  {
    return allNodes_;
  }

  NodeNamePairMap &getAllNodes()
  {
    return allNodes_;
  }

  const NodeNamePairMap &getStateNodes() const
  {
    return stateNodes_;
  }

  const NodeNamePairMap &getStoreNodes() const
  {
    return storeNodes_;
  }

  void setOutputFilenameSuffix( std::string newSuffix )
  {
    filenameSuffix_ = newSuffix;
  };

  const std::string &getFilenameSuffix() const
  {
    return filenameSuffix_;
  }

  const std::string &getTitle() const
  {
    return title_;
  }

  void setTitle(const std::string &title)
  {
    title_ = title;
  }

  const std::string &getNetListFilename() const
  {
    return netListFilename_;
  }

  bool setupInitialConditions(
    Parallel::Machine                   comm,
    N_LAS_Vector &                      solnVec,
    N_LAS_Vector &                      flagVec);

  NodeNamePairMap &getICData( int &found, int &ic_type);

  void outputDCOP(
    Parallel::Machine                   comm,
    const N_LAS_Vector &                solution);

  // Runs specified output commands
  void output(
    Parallel::Machine                   comm,
    const double                        time,
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const int                           dcNumber,
    const int                           maxDC,
    const Analysis::SweepVector &       dcParamVec1,
    const N_LAS_Vector &                solnVecPtr,
    const N_LAS_Vector &                stateVecPtr,
    const N_LAS_Vector &                storeVecPtr,
    const std::vector<double> &         objectiveVec,
    const std::vector<double> &         dOdpVec,
    const std::vector<double> &         dOdpAdjVec,
    const std::vector<double> &         scaled_dOdpVec,
    const std::vector<double> &         scaled_dOdpAdjVec,
    bool                                skipPrintLineOutput = false);

  void outputAC(
    Parallel::Machine                   comm,
    double                              freq,
    const N_LAS_Vector &                freqDomainSolnVecReal,
    const N_LAS_Vector &                freqDomainSolnVecImaginary);

  void outputMPDE(
    Parallel::Machine                   comm,
    double                              time,
    const std::vector<double> &         fast_time_points,
    const N_LAS_Vector &                solution_vector);

  void outputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         param_values,
    const N_LAS_Vector &                solution_vector);

  void outputHB(
    Parallel::Machine                   comm,
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const std::vector<double> &         timePoints,
    const std::vector<double> &         freqPoints,
    const N_LAS_BlockVector &           timeDomainSolutionVec,
    const N_LAS_BlockVector &           freqDomainSolutionVecReal,
    const N_LAS_BlockVector &           freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &           timeDomainStoreVec,
    const N_LAS_BlockVector &           freqDomainStoreVecReal,
    const N_LAS_BlockVector &           freqDomainStoreVecImaginary);

  void finishOutput();

  void startStep(
    int                           step,
    int                           max_step);

  void steppingComplete();

  // void getMeasureValue (const std::string &name, double &result, bool &found) const;

  std::vector<std::string> getVariableNames();

  void setExternalNetlistParams(const StringPairVector &externalNetlistParams);

  void setAliasNodeMap(const AliasNodeMap &alias_node_map) {
    aliasNodeMap_ = alias_node_map;
  }

  const IO::AliasNodeMap &getAliasNodeMap() const {
    return aliasNodeMap_;
  }

  const Analysis::SweepVector &getStepParamVec() const
  {
    return outputState_.stepSweepVector_;
  }

  int getStepLoopNumber() const
  {
    return outputState_.stepLoopNumber_;
  }

  int getMaxParamSteps() const
  {
    return outputState_.stepMaxCount_;
  }

  const Analysis::SweepVector &getDCParamVec() const
  {
    return outputState_.dcSweepVector_;
  }

  int getDCLoopNumber() const
  {
    return dcLoopNumber_;
  }

  int getMaxDCSteps() const
  {
    return maxDCSteps_;
  }

  bool getTempSweepFlag() const
  {
    return tempSweepFlag_;
  }

  ParameterList getVariableList() const;

  double getPRINTDCvalue() const
  {
    return PRINTdcvalue_;
  }

  double getPRINTDCstart() const
  {
    return PRINTdcstart_;
  }

  double getPRINTDCstop() const
  {
    return PRINTdcstop_;
  }

  const Analysis::AnalysisManager *getAnaIntPtr() const
  {
    return analysisManager_;
  }

  const std::string &getPRINTDCname() const
  {
    return PRINTdcname_;
  }

  bool getPrintEndOfSimulationLine() const
  {
    return printEndOfSimulationLine_;
  }

  bool getOutputVersionInRawFile() const
  {
    return outputVersionInRawFile_;
  }


  void setEnableHomotopyFlag(bool value)
  {
    enableHomotopyFlag_ = value;
  }

  double getCircuitTime() const
  {
    return outputState_.circuitTime_;
  }

  void setCircuitTime(double time)
  {
    outputState_.circuitTime_ = time;
  }

  double getCircuitFrequency() const
  {
    return outputState_.circuitFrequency_;
  }

  void setCircuitFrequency(double frequency)
  {
    outputState_.circuitFrequency_ = frequency;
  }

  double getCircuitTemp() const
  {
    return outputState_.circuitTemp_;
  }

  bool getSTEPEnabledFlag() const
  {
    return STEPEnabledFlag_;
  }

  const ObjectiveMap &getObjectiveMap() const
  {
    return objectiveMap_;
  }

  double getTemperature() const
  {
    return outputState_.circuitTemp_;
  }

  double getTime() const
  {
    return outputState_.circuitTime_;
  }

  double getFrequency() const
  {
    return outputState_.circuitFrequency_;
  }

  double getStepSweep(size_t index) const
  {
    return outputState_.stepSweepVector_[index].currentVal;
  }

  double getDCSweep(size_t index) const
  {
    return outputState_.dcSweepVector_[index].currentVal;
  }

  const StringPairVector &getResponseFunctions() const {
    return responseFunctionsRequested_;
  }

  void pushActiveOutputters()
  {
    activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());
  }

  void popActiveOutputters()
  {
    activeOutputterStack_.pop_back();
  }

  void addActiveOutputter(PrintType::PrintType print_type, Analysis::Analysis_Mode analysis_mode)
  {
    OutputterMap::iterator find_it = outputterMap_.find(print_type);

    if (find_it == outputterMap_.end()) {

    }

    if (!activeOutputterStack_.empty() && find_it != outputterMap_.end()) 
    {
      for (OutputterMap::mapped_type::iterator 
          it = (*find_it).second.begin(), end = (*find_it).second.end(); it != end; ++it)
      {
        (*it)->setAnalysisMode(analysis_mode);
      }

      activeOutputterStack_.back().insert(activeOutputterStack_.back().end(), (*find_it).second.begin(), (*find_it).second.end());
    }
  }

  void addOutputPrintParameters(OutputType::OutputType output_type, const PrintParameters &print_parameters);

  Objective &getObjective(const std::string &varType) const
  {
    ObjectiveMap::const_iterator it = objectiveMap_.find(varType);
    return const_cast<Objective &>((*it).second);
  }

  std::pair<OutputParameterMap::const_iterator, bool> findOutputParameter(OutputType::OutputType output_type) const {
    OutputParameterMap::const_iterator it = outputParameterMap_.find(output_type);
    return std::pair<OutputParameterMap::const_iterator, bool>(it, it != outputParameterMap_.end());
  }

  void addOutputter(PrintType::PrintType print_type, Outputter::Interface *outputter) {
    outputterMap_[print_type].push_back(outputter);
  }

  const PrintParameters &getDefaultPrintParameters() const {
    return defaultPrintParameters_;
  }

  const Util::Op::BuilderManager &getOpBuilderManager() const {
    return opBuilderManager_;
  }

  Util::Op::BuilderManager &getOpBuilderManager() {
    return opBuilderManager_;
  }

  std::ostream *openFile(const std::string &path, std::ios_base::openmode mode);
  std::ostream *openFile(const std::string &path);
  std::ostream *openBinaryFile(const std::string &path);
  int closeFile(std::ostream *os);

  mutable std::set<std::string> deferredParameterCheck_;

  bool prepareHDF5Output(Parallel::Machine comm);
  bool updateHDF5Output(Parallel::Machine comm, const N_LAS_Vector &solnVecPtr);
  bool closeHDF5Output();

private:
  std::string                   title_;
  std::string                   netListFilename_;
  std::string                   filenameSuffix_;

  Util::Op::BuilderManager &    opBuilderManager_;
  Measure::Manager *            measureManager_;
  FourierMgr *                  fourierManager_;
  OutputterMap                  outputterMap_;
  OutputParameterMap            outputParameterMap_;
  ActiveOutputterStack          activeOutputterStack_;
  EnabledAnalysisSet            enabledAnalysisSet_;

  InitialConditionsData         icData_;

  struct OutputState
  {
    OutputState()
      : circuitTime_(0.0),
        circuitTemp_(27.0),
        circuitFrequency_(0.0),
        stepLoopNumber_(0),
        stepMaxCount_(0)
    {}

    double                      circuitTime_;           ///< transient circuit time:
    double                      circuitTemp_;           ///< circuit temperature
    double                      circuitFrequency_;      ///< ac current frequency
    int                         stepLoopNumber_;        ///< step current loop
    int                         stepMaxCount_;          ///< step max count

    Analysis::SweepVector       stepSweepVector_;
    Analysis::SweepVector       dcSweepVector_;
  };

  OutputState                   outputState_;

  std::vector<std::string>      stepParams_;
  std::vector<std::string>      dcParams_;

  Topo::Topology *              topology_;
  Device::DeviceInterface *     deviceInterface_;
  Analysis::AnalysisManager *   analysisManager_;

  // print statement vars
  bool                  enableHomotopyFlag_;
  bool                  enableSensitivityFlag_;
  int                   sensitivityOptions_;
  ParameterList         sensitivityVariableList_;

  PrintParameters       defaultPrintParameters_;

  double                PRINTdcstart_;
  double                PRINTdcstop_;
  double                PRINTdcvalue_;
  std::string           PRINTdcname_;

  // initial condition vars
  bool loadFlag_;
  bool outputOnceAlreadyFlag_;
  std::vector<Util::OptionBlock> ICblockVec_;
  std::vector<Util::OptionBlock> nodesetblockVec_;

  double initialOutputInterval_;
  std::vector<std::pair< double, double > > outputIntervalPairs_;

  bool                  STEPEnabledFlag_;

  bool                  printEndOfSimulationLine_;  // flag to indicate if user wants the "End of Xyce(TM)" line in the output.
  bool                  outputVersionInRawFile_;    // flag to indicate that Version should be output in the header of a RAW file.

  // response functions requested in the external params passed to Xyce.
  std::vector< std::pair< std::string, std::string> > variablesUsedInSimulation_ ;
  std::vector< std::pair< std::string, std::string> > responseFunctionsRequested_ ;
  std::vector< std::pair< std::string, std::string> > derivativeVariablesRequested_ ;
  std::vector< std::pair< std::string, std::string> > analysisComponentsRequested_ ;

  bool tempSweepFlag_;
  bool outputCalledBefore_;

  ObjectiveMap objectiveMap_;

  // dc loop information
  int dcLoopNumber_;
  int maxDCSteps_;

  bool isStarredPrintLineProcessed;

  std::vector<char>     typeRefs_;
  NodeNamePairMap       allNodes_;
  NodeNamePairMap       stateNodes_;
  NodeNamePairMap       storeNodes_;
  NodeNamePairMap       externNodes_;
  AliasNodeMap          aliasNodeMap_;

  OutputterFactoryMap   outputterFactoryMap_;
  OpenPathStreamMap     openPathStreamMap_;

  // HDF5 file support vars
  bool hdf5FileNameGiven_;
  bool hdf5HeaderWritten_;
  std::string hdf5FileName_;
  int hdf5IndexValue_;

#ifdef Xyce_USE_HDF5
  hid_t hdf5FileId_;
  hid_t hdf5PlistId_;
#endif  // Xyce_USE_HDF5

};

inline bool operator<(const OutputMgr::OutputterKey &lhs, const OutputMgr::OutputterKey &rhs) {
  return (lhs.analysisMode_ < rhs.analysisMode_)
    || (!(rhs.analysisMode_ < lhs.analysisMode_) && lhs.outputType_ < rhs.outputType_)
    || (!(rhs.analysisMode_ < lhs.analysisMode_ && rhs.outputType_ < lhs.outputType_) && lhs.format_ < rhs.format_);
}

void printLineDiagnostics(
  Parallel::Machine                                                     comm,
  const OutputMgr &                                                     output_manager,
  const Measure::Manager &                                              measure_manager,
  const std::set<std::string> &                                         node_names,
  const std::map<std::string, RefCountPtr<Device::InstanceBlock> > &    device_map,
  IO::AliasNodeMap &                                                    alias_node_map);

void deferredPrintLineDiagnostics(
  Parallel::Machine     comm,
  OutputMgr &           output_manager);

// if any macro-simulation level functions were defined, finish them and output results
void outputMacroResults(Parallel::Machine comm, const OutputMgr &output_manager);

NodeNamePairMap::const_iterator findNode(const std::string &name, const NodeNamePairMap &node_map, const AliasNodeMap &alias_map);

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::OutputMgr N_IO_OutputMgr;

#endif // Xyce_N_IO_OutputMgr_h
