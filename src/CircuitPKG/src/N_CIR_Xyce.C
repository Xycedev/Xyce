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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_CIR_Xyce.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.312.2.14 $
//
// Revision Date  : $Date: 2014/08/28 22:41:46 $
//
// Current Owner  : $Author: rlschie $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <ctime>

// ----------   Xyce Includes   ----------

#include <N_CIR_Xyce.h>

#include <N_DEV_fwd.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DAC.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternalSimulationData.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Print.h>
#include <N_DEV_Device.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_ADC.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_NetlistImportTool.h>

#include <N_IO_OutputMgr.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_OutputMacroResults.h>
#include <N_IO_OutputResponse.h>
#include <N_IO_PrintDeviceCount.h>
#include <N_IO_RestartMgr.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_NLS_Manager.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_ANP_AnalysisManager.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_TOP_Topology.h>
#include <N_TOP_TopologyMgr.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Timer.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Expression.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_PrintStats.h>
#include <N_UTL_Platform.h>
#include <N_UTL_JSON.h>
#include <N_UTL_SendCURL.h>

#include <N_UTL_Version.h>
#include <N_UTL_BreakPoint.h>

class N_LAS_QueryUtil;
class N_LAS_MultiVector;

namespace Xyce {
namespace Circuit {

#ifdef Xyce_TRACKING_URL
static const char *trackingURL = Xyce_TRACKING_URL;
#else
static const char *trackingURL = 0;
#endif

//--------  Global Declarations ------------
void report_handler(const char *message, unsigned report_mask);

namespace {

struct ADCDeviceInstanceParameterOp: public Device::DeviceInstanceOp
{
  ADCDeviceInstanceParameterOp(std::map<std::string, std::map<std::string, double> > &adc_device_parameter_map)
      : ADCDeviceParameterMap_(adc_device_parameter_map)
  {}

  virtual bool operator()(Device::DeviceInstance *instance) {
    Device::ADC::Instance *adc_instance = dynamic_cast<Device::ADC::Instance *>(instance);
      if (adc_instance) {
        adc_instance->getModel().getParam("LOWERVOLTAGELIMIT", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["lowerVoltageLimit"]);
        adc_instance->getModel().getParam("UPPERVOLTAGELIMIT", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["upperVoltageLimit"]);
        adc_instance->getModel().getParam("SETTLINGTIME", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["settlingTime"]);
      }

      return true;
    }

  std::map<std::string, std::map<std::string, double> > &ADCDeviceParameterMap_;
};

struct TimeVoltagePairsOp: public Device::DeviceInstanceOp
{
  TimeVoltagePairsOp(std::map<std::string, std::vector< std::pair<double,double> > >&time_voltage_map)
    : TimeVoltageMap_(time_voltage_map)
  {}

  virtual bool operator()(Device::DeviceInstance *instance) {
    Device::ADC::Instance &adc_instance = static_cast<Device::ADC::Instance &>(*instance);

    std::vector<std::pair <double,double> > TmpVec;
    adc_instance.getTVVEC(TmpVec);
    double current_time = adc_instance.getSolverState().currTime;
    double * solVector = adc_instance.getExternData().nextSolVectorRawPtr;

    double vPos = solVector[adc_instance.getLIPos()];
    double vNeg = solVector[adc_instance.getLINeg()];
    TmpVec.push_back(std::pair<double,double>(current_time, vPos - vNeg));

    TimeVoltageMap_[instance->getName().getEncodedName()] = TmpVec;

    return true;
  }

  std::map<std::string, std::vector< std::pair<double,double> > > &      TimeVoltageMap_;
};

} // namespace <unnamed>


//-----------------------------------------------------------------------------
// Function      : Simulator::Simulator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::Simulator(Parallel::Machine comm)
:
  comm_(comm),
  devIntPtr_(0),
  topPtr_(0),
  topMgrPtr_(0),
  netlistImportToolPtr_(0),
  outputManager_(0),
  measureManager_(0),
  fourierManager_(0),
  outputResponse_(0),
  lasSysPtr_(0),
  lasBuilderPtr_(0),
  analysisManager_(0),
  nlsMgrPtr_(0),
  parMgrPtr_(0),
  resMgrPtr_(0),
  rootStat_(Stats::createRootStat("Xyce", Stats::StatSet(Stats::STAT_ALL))),
  XyceTimerPtr_(0),
  ElapsedTimerPtr_(0),
  multiThreading_(false),
  numThreads_(0),
  initializeAllFlag_(false),
  commandLine()
{
  previousReportHandler_ = set_report_handler(report_handler);

  rootStat_.start();

  Device::registerDevices();
}

//-----------------------------------------------------------------------------
// Function      : ~Xyce
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::~Simulator()
{
  set_report_handler(previousReportHandler_);

  Stats::deleteRootStat(rootStat_);
}

//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : This passes a vector of pairs "key" "value" that will
//                 be substituted during the processing of the netlist.  This
//                 more easily allows Dakota to change any netlist parameter
//                 during netlist setup.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setNetlistParameters( const std::vector< std::pair< std::string, std::string > > & externalParams )
{
  externalNetlistParams_ = externalParams;
}

//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : Call through to the output manager to set the suffix to
//                 be used on the output file, as in circuit + suffix + prn
//                 This is useful in Dakota controlled runs to keep each
//                 simulation from overwritting the last one.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setOutputFileSuffix( const std::string newSuffix )
{
  if( outputManager_ )
  {
    outputManager_->setOutputFilenameSuffix( newSuffix );
  }
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setupParMgr_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool Simulator::setupParMgr_( int iargs, char **cargs)
{
  bool bsuccess = true;

  // Setup the Parallel Mgr. with a default load-balance based on the numProc
  // value.
  parMgrPtr_ = new N_PDS_Manager( isSerialFlag_, procZeroFlag_, iargs, cargs, comm_);
  if (comm_ == 0)
    comm_ = parMgrPtr_->getPDSComm()->comm();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doAllocations_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doAllocations_()

{
  std::string topotype = "Basic";

  // Allocate device manager:
  devIntPtr_  = Device::DeviceInterface::factory(commandLine);

  // Allocate Topology:
  topMgrPtr_ = Topo::Manager::instance();
  Topo::Manager & topMgr = *topMgrPtr_;

  topPtr_ = topMgrPtr_->createTopology(commandLine);

  // Allocate distribution mgr:
  netlistImportToolPtr_ = IO::NetlistImportTool::factory(commandLine, topMgr);

  // Allocate output mgr:
  opBuilderManager_ = new Util::Op::BuilderManager();
  outputManager_ = new IO::OutputMgr(commandLine, *opBuilderManager_);
  IO::registerOpBuilders(*opBuilderManager_, comm_, *outputManager_);
  outputResponse_ = new IO::OutputResponse(*outputManager_);
  measureManager_ = new IO::Measure::Manager(commandLine.getArgumentValue("netlist"));
  IO::registerOpBuilders(*opBuilderManager_, comm_, *measureManager_);
  fourierManager_ = new IO::FourierMgr(commandLine.getArgumentValue("netlist"));

  // Allocate restart mgr:
  resMgrPtr_ = IO::RestartMgr::factory(commandLine);

  // Linear Algebra allocations:
  lasSysPtr_       = new N_LAS_System();
  lasBuilderPtr_   = new N_LAS_Builder();

  // Allocate Time Integration:
  analysisManager_    = new Analysis::AnalysisManager(commandLine, rootStat_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *measureManager_);

  // Allocate nonlinear solver:
  nlsMgrPtr_ = new Nonlinear::Manager(commandLine);

  pkgOptionsMgrPtr_ = new IO::PkgOptionsMgr();

  bool bsuccess = true;
  bsuccess = bsuccess && (devIntPtr_ != 0);
  bsuccess = bsuccess && (topPtr_ != 0);
  bsuccess = bsuccess && (netlistImportToolPtr_ != 0 );
  bsuccess = bsuccess && (outputManager_ != 0);
  bsuccess = bsuccess && (lasSysPtr_ != 0);
  bsuccess = bsuccess && (lasBuilderPtr_ != 0);
  bsuccess = bsuccess && (analysisManager_ != 0);
  bsuccess = bsuccess && (nlsMgrPtr_ != 0);
  bsuccess = bsuccess && (parMgrPtr_ != 0);
  bsuccess = bsuccess && (pkgOptionsMgrPtr_ != 0);

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Simulator::doDeAllocations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doDeAllocations_()

{
  // de-allocate the device manager:

  delete devIntPtr_;
  delete nlsMgrPtr_;
  delete analysisManager_;
  delete outputManager_;
  delete measureManager_;
  delete fourierManager_;
  delete outputResponse_;
  delete opBuilderManager_;
  delete lasSysPtr_;
  delete lasBuilderPtr_;
  delete XyceTimerPtr_;
  delete ElapsedTimerPtr_;
  delete parMgrPtr_;
  delete topMgrPtr_;
  delete resMgrPtr_;
  delete pkgOptionsMgrPtr_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doRegistrations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doRegistrations_()

{
  bool bsuccess = true;
  bool bs1 = true;


  // ParallelDist Manager registrations
  bs1 = parMgrPtr_->registerTopology(topPtr_);
  bsuccess = bsuccess && bs1;

  // Device Manager registrations
  bs1 = devIntPtr_->registerLinearSystem(lasSysPtr_);    bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerNonlinearSolver(nlsMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerAnalysisManager(analysisManager_);  bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerParallelMgr(parMgrPtr_);     bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerOutputMgr(outputManager_);       bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerMeasureMgr(measureManager_);       bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Topology registrations:
  bs1 = topPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registerParallelMgr(parMgrPtr_);     bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registerAnalysisManager(analysisManager_);          bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Distribution manager registrations:
  bs1 = netlistImportToolPtr_->registerParallelServices(parMgrPtr_->getPDSComm());
  bsuccess = bsuccess && bs1;

  bs1 = netlistImportToolPtr_->registerDevMgr(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = netlistImportToolPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Restart manager registrations:
  bs1 = resMgrPtr_->registerTopology(topPtr_);           bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerAnalysisManager(analysisManager_);          bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerParallelServices(parMgrPtr_);bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Output manager registrations:
  bs1 = outputManager_->registerTopology(topPtr_);           bsuccess = bsuccess && bs1;
  bs1 = outputManager_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = outputManager_->registerAnalysisManager(analysisManager_);          bsuccess = bsuccess && bs1;
  bs1 = outputManager_->registerPkgOptionsMgr( *pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = measureManager_->registerPkgOptionsMgr( *pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = fourierManager_->registerPkgOptionsMgr( *pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  outputManager_->registerMeasureMgr(measureManager_);
  outputManager_->registerFourierMgr(fourierManager_);

  // Analysis manager registrations:
  bs1 = analysisManager_->registerLinearSystem(lasSysPtr_); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerNLSManager(nlsMgrPtr_);   bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerParallelServices(parMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerOutputMgr(outputManager_);    bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerElapsedTimer(ElapsedTimerPtr_); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerRestartMgr(resMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerDeviceInterface(devIntPtr_);  bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerTopology(topPtr_);            bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerApplicationBuilder(lasBuilderPtr_); bsuccess = bsuccess && bs1;

  // Nonlinear Solver registrations:
  bs1 = nlsMgrPtr_->registerLinearSystem(lasSysPtr_);   bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerAnalysisManager(analysisManager_); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerOutputMgr(outputManager_);      bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerTopology(topPtr_); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerParallelMgr(parMgrPtr_); bsuccess = bsuccess && bs1;

  // Linear Solver registrations:
  bs1 = lasSysPtr_->registerPDSManager(parMgrPtr_);    bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerBuilder( lasBuilderPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerQueryUtil(
    (N_LAS_QueryUtil *) topPtr_->get_LinSolvUtil() );  bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerAnalysisManager(analysisManager_);    bsuccess = bsuccess && bs1;

  bs1 = lasBuilderPtr_->registerPDSManager(parMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = lasBuilderPtr_->registerSystem(lasSysPtr_);     bsuccess = bsuccess && bs1;
  bs1 = lasBuilderPtr_->registerQueryUtil(
    (N_LAS_QueryUtil *) topPtr_->get_LinSolvUtil() );
  bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setUpTopology_
// Purpose       : This function handles a lot of the initial setup.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool Simulator::setUpTopology_()
{
  std::string netListFile;

  if( procZeroFlag_ )
  {
    if (commandLine.getArgumentValue("netlist") != "")
      netListFile = commandLine.getArgumentValue("netlist");
    else
      netListFile = "xyce.in";

    FILE * testFile;
    if ( (testFile = fopen(netListFile.c_str(), "r")) == 0)
    {
      Report::UserError() << netListFile << " file not found";
    }
    else
    {
      fclose( testFile );
    }
  }

  Report::safeBarrier(comm_);

  Xyce::lout() << "***** Reading and parsing netlist..." << std::endl;

  netlistImportToolPtr_->constructCircuitFromNetlist(netListFile, externalNetlistParams_, *outputManager_, *measureManager_);

  Report::safeBarrier(comm_);

  if ( commandLine.argExists("-syntax") || commandLine.argExists("-count") )
  {
    if (commandLine.argExists("-syntax"))
    {
      Xyce::lout() << "***** Netlist syntax OK\n";
    }

    Xyce::lout() << std::endl;

    Xyce::lout() << "***** Device Type Counts ...\n" << std::endl;

    IO::printDeviceCount(Xyce::lout(), devIntPtr_->getDeviceCountMap());

    Xyce::lout() << std::endl;

    reportTotalElapsedTime ();

    Xyce_exit(0);
  }

  outputManager_->setAliasNodeMap(netlistImportToolPtr_->getAliasNodeMap());

  delete netlistImportToolPtr_;
  netlistImportToolPtr_ = 0;

  Xyce::lout() << "***** Setting up topology...\n" << std::endl;

  // topology query's device manager to see if any devices are bad (i.e. a resistor with zero resistance)
  // if so, a list of nodes to be supernoded is created
  topPtr_->verifyNodesAndDevices();

#ifdef Xyce_PARALLEL_MPI
  // create a union of the supernode list on all processors
  topPtr_->mergeOffProcTaggedNodesAndDevices();
#endif

// combine nodes into supernodes and remove now redundant devices (i.e. those only connected to 1 processor )
  topPtr_->removeTaggedNodesAndDevices();

  // if "-remeasure" was on the command line, then we don't need to
  // instantiate the devices.
  if (commandLine.argExists("-remeasure"))
  {
    measureManager_->remeasure(*parMgrPtr_->getPDSComm(), commandLine.getArgumentValue("netlist"), commandLine.getArgumentValue("-remeasure"), analysisManager_->getTransientFlag(), *outputManager_);
    Xyce::lout() << "***** Remeasure analysis complete\n" << std::endl;
    return false;
  }
  topPtr_->instantiateDevices();

  IO::deferredPrintLineDiagnostics(comm_, *outputManager_);

  Report::safeBarrier(comm_);

#ifdef Xyce_PARALLEL_MPI
  devIntPtr_->setGlobalFlags();
#endif

  // Setup of indices including global reordering.
  topPtr_->setupGlobalIndices();

  Report::safeBarrier(comm_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setUpMatrixStructure_
// Purpose       : This function needs to set up the various linear algebra
//                 entities which are owned by the LAS system class.  This
//                 includes, most importantly, the Jacobian matrix and the
//                 right hand side vector.  It should also set the solution
//                 vector size and the state vector size.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool Simulator::setUpMatrixStructure_()
{
  lasBuilderPtr_->generateParMaps();
  lasBuilderPtr_->generateGraphs();

  lasSysPtr_->initializeSystem();

  topPtr_->registerLIDswithDevs();
  topPtr_->registerJacLIDswithDevs();

  devIntPtr_->setupExternalDevices();

  int lasSize = lasSysPtr_->getGlobalSolutionSize();
  Xyce::lout() << "***** Number of Unknowns = " << lasSize << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doInitializations_
// Purpose       : This function calls "initializeAll" functions in all the
//                 packages which have them.  Packages that have the LAS system
//                 class registered with them will have one of these functions.
//
// Special Notes : This is called once the size of the linear system is known,
//                 as various packages will need to allocate vectors, matrices,
//                 etc.
//
//                 These probably can be done in any order, but to be safe,
//                 make sure for now that the time integrator's initializeAll
//                 function is called first.  The other two are essentially
//                 secondary registrations, while the TIA one includes a lot of
//                 allocations.
//
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::doInitializations_()
{
  bool bsuccess = true;
  bool bs1 = true;

  bs1 = analysisManager_->initializeAll();  bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->initializeAll();  bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->initializeAll();  bsuccess = bsuccess && bs1;

  if( resMgrPtr_->isRestart() ) resMgrPtr_->restoreRestartData();

  topPtr_->generateICLoader();

  initializeAllFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSolvers_
// Purpose       : This function runs the solvers.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::runSolvers_()
{
  return analysisManager_->run();
}


//-----------------------------------------------------------------------------
// Function      : Simulator::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool Simulator::run(int argc, char *argv[])
{
  bool bsuccess = true;

  try
  {
    bool bs1 = initialize(argc, argv);
    bsuccess = bsuccess && bs1;
  }
  catch (std::exception &x)
  {
    Xyce::lout() << "Exception " << x.what() << std::endl;
    bsuccess = false;
  }

  if (!bsuccess)
  {
    reportTotalElapsedTime ();
    Xyce::lout() << "Xyce Initialization Phase failed." << std::endl;
    Xyce_exit(-1);
  }

  try
  {
    bool bs1 = runSimulation();
    bsuccess = bsuccess && bs1;

    bs1 = finalize();
    bsuccess = bsuccess && bs1;
  }
  catch (std::exception &x)
  {
    Xyce::lout() << "Exception " << x.what() << std::endl;
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSimulation
// Purpose       : Main simulation driver.
// Special Notes : Not private as this is also called from N_DAK_DakotaInterface
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::runSimulation()
{
  return runSolvers_();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR MIXED-SIGNAL and other external applications
//
//-----------------------------------------------------------------------------
// Function      : Simulator::initialize
// Purpose       : capture all "initialization-type" activities in one
//                 method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
bool Simulator::initialize( int argc, char **argv)
{
  argc_ = argc;
  argv_ = argv;

  // Setup the Parallel Manager first to give us a N_PDS_Comm object so the
  // reporting package will work.
  bool bsuccess = setupParMgr_(argc, argv);
  bool b2;

  // Start the solvers timer.
  ElapsedTimerPtr_ = new Util::Timer( *(parMgrPtr_->getPDSComm()) );

  const time_t now=time(NULL);
  char timeDate[40];
  strftime(timeDate,40,"%x %X %Z",localtime(&now));

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD
  N_PDS_Comm * commPtr  = parMgrPtr_->getPDSComm();
  size_t dirsize = 512;
  char directory[dirsize];  for (int idir=0;idir<dirsize;++idir)  directory[idir] = 0;
#endif
#endif

  // register parallel mgr to allow parallel support
  commandLine.registerParallelMgr( parMgrPtr_ );

  initializeLogStream(parMgrPtr_->getPDSComm()->procID(), parMgrPtr_->getPDSComm()->numProc());

  // read in command line arguments
  int status = commandLine.parseCommandLine(argc, argv );
  if (status != 0)
  {
    Xyce_exit(status == -1 ? 0 : status);
  }

  // Set the output stream of the "-l" flag exists
  if (commandLine.argExists("-l"))
  {
    openLogFile(commandLine.getArgumentValue("-l"), commandLine.argExists("-per-processor"));
  }

  if( procZeroFlag_ )
  {
    Xyce::lout() << "\n"
                 << "*****\n"
                 << "***** Welcome to the Xyce(TM) Parallel Electronic Simulator\n"
                 << "*****\n"
                 << "***** This is version " << Util::Version::getFullVersionString() << "\n\n\n"
                 << "***** Executing netlist " << (commandLine.getArgumentValue("netlist") == "" ? "MISSING!" : commandLine.getArgumentValue("netlist")) << "\n\n";

    // Don't bother doing anything else if we weren't given a netlist!
    if (commandLine.getArgumentValue("netlist") == "")
    {
      Xyce::lout() << "  Usage: " << argv[0] << " netlist\n" << std::endl;
      return false;
    }

    // check if a parameter file was specified for param substitution during parsing
    if (commandLine.argExists(std::string("-prf")))
    {
      std::string parameterFile = commandLine.getArgumentValue("-prf");
      readExternalParamsFromFile( parameterFile, externalNetlistParams_ );
#ifdef Xyce_DEBUG
      // debug print out externNetlistParams_
      std::vector<std::pair<std::string,std::string> >::iterator currentPair = externalNetlistParams_.begin();
      std::vector<std::pair<std::string,std::string> >::iterator endPair = externalNetlistParams_.end();
      while( currentPair != endPair )
      {
        dout() << "\"" << currentPair->first << "\" = \"" << currentPair->second << "\"" << std::endl;
        currentPair++;
      }
#endif
    }
  }

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD

  if( procZeroFlag_ )
  {
    (void) getcwd(directory, dirsize);
  }
  commPtr->bcast(directory, dirsize, 0);
  if( !procZeroFlag_ )
  {
    chdir(directory);
  }

#endif
#endif

  // Allocate all the various packages:

  b2 = doAllocations_();
  bsuccess = bsuccess && b2;

  if (b2)
  {
    if (DEBUG_CIRCUIT)
    {
      dout() << "Allocation was successful.";
    }
  }
  else
  {
    Report::DevelFatal() << "Allocation was NOT successful.";
  }

  // Register the external package pointers:
  b2 = doRegistrations_();
  bsuccess = bsuccess && b2;

  if (b2)
  {
    if (DEBUG_CIRCUIT)
      dout() << "Registration was successful.";
  }
  else
    Report::DevelFatal() << "Registration was NOT successful.";

  b2 = setUpTopology_();
  if( !b2 )
  {
    return false;
  }
  bsuccess = bsuccess && b2;

  Xyce::lout() << "***** Device Count Summary ..." << std::endl;

  {
    IO::DeviceCountMap global_device_count_map;

    IO::gatherGlobalDeviceCount(comm_, global_device_count_map, devIntPtr_->getDeviceCountMap());

    auditJSON_ << Util::nameValuePair("StartTime", timeDate) << Util::JSON::sep
               << Util::nameValuePair("DeviceCount", global_device_count_map) << Util::JSON::sep
               << Util::nameValuePair("DeviceTotalCount", std::accumulate(global_device_count_map.begin(), global_device_count_map.end(), 0, IO::DeviceCountMapSum()));

    IO::printDeviceCount(Xyce::lout(), global_device_count_map);
  }

  if( commandLine.argExists( "-norun" ) || commandLine.argExists( "-namesfile" )  )
  {
    Xyce::lout() << "\n***** Syntax and topology analysis complete" << std::endl;

    if( commandLine.argExists( "-namesfile" ) )
    {
      setUpMatrixStructure_();
      topPtr_->outputNameFile(true);
    }

    return false;
  }
  else
  {
    Xyce::lout() << "\n***** Setting up matrix structure..." << std::endl;
    b2 = setUpMatrixStructure_();
    bsuccess = bsuccess && b2;

    Xyce::lout() << "***** Initializing...\n" << std::endl;
    b2 = doInitializations_();
    bsuccess = bsuccess && b2;

    // optional diagnostic output file:
    topPtr_->outputNameFile();

    outputManager_->checkPrintParameters(comm_);

    // if we loaded a parameter file from the command line, then the output manager
    // should scan the results as well as such files can specify response functions
    // than need to be reported by the output manager.
    if (commandLine.argExists(std::string("-prf")))
    {
      outputManager_->setExternalNetlistParams( externalNetlistParams_ );
    }
    if (commandLine.argExists(std::string("-rsf")))
    {
      std::string responseFile = commandLine.getArgumentValue("-rsf");
      outputResponse_->setResponseFilename( responseFile );
    }
    // Start the solvers timer.
    XyceTimerPtr_ = new Util::Timer( *(parMgrPtr_->getPDSComm()) );
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getDeviceNames
// Purpose       : get all names of devices of specified type in the netlist
// Special Notes : "deviceType" takes a string of the form the device would
//                 have when instantiated in a netlist, e.g. "R" for a resistor
//                 or "YXYGRA" for a Xygra device
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/25/08
//-----------------------------------------------------------------------------
bool Simulator::getDeviceNames(const std::string &modelGroupName, std::vector<std::string> &deviceNames)
{
  Device::EntityTypeId model_group = devIntPtr_->getModelGroup(modelGroupName);

  if (!model_group.defined())
    model_group = devIntPtr_->getModelGroup(Device::InstanceName(modelGroupName).getDeviceType());

  if (!model_group.defined())
    return false;

  Device::Device *device = devIntPtr_->getDevice(model_group);
  if (!device)
    return false;

  Device::getDeviceInstanceNames(*device, std::back_inserter(deviceNames));

  return true;
}

//----------------------------------------------------------------------------
// Function       : Simulator::getDACDeviceNames
// Purpose        : Gets the (stripped) names of the DAC devices
//                  in the circuit.
// Special Notes  :
// Scope          :
// Creator        : Lisa Maynes
// Creation Date  : 06/13/2003
//----------------------------------------------------------------------------
bool Simulator::getDACDeviceNames(std::vector< std::string >& dacNames)
{
  dacNames.clear();

  Device::Device *device = devIntPtr_->getDevice(Device::DAC::Traits::modelGroup());
  if (device)
    Device::getDeviceInstanceNames(*device, std::back_inserter(dacNames));

  return true;
}


//----------------------------------------------------------------------------
// Function       : getADCMap
// Purpose        : Gets the (stripped) names of the ADC devices
//                 in the circuit(as key of map) and map of parameters
//                 (keyed by parameter name) for each device
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::getADCMap(std::map<std::string, std::map<std::string, double> >&ADCMap)
{
  ADCDeviceInstanceParameterOp op(ADCMap);

  Device::Device *device = devIntPtr_->getDevice(Device::ADC::Traits::modelGroup());
  if (device)
    device->forEachInstance(op);

  return true;
}

//----------------------------------------------------------------------------
// Function       : updateTimeVoltagePairs
// Purpose        : Update the DAC devices in a circuit by adding the set
//                  of time and voltage pairs built up on the "digital side"
//                  since the last update and by removing the time-voltage
//                  pairs for times that pre-date the given simulation time.
// Special Notes  : The current method for locating DAC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 06/09/2003
//----------------------------------------------------------------------------
bool Simulator::updateTimeVoltagePairs(const std::map< std::string, std::vector<std::pair<double,double> > *> & timeVoltageUpdateMap)
{
  for (std::map<std::string, std::vector< std::pair<double,double> >* >::const_iterator it = timeVoltageUpdateMap.begin(), end = timeVoltageUpdateMap.end(); it != end; ++it) {
    const std::string &dacName = (*it).first;
    const std::vector<std::pair<double,double> > &tv_pair_vector = *(*it).second;

    Device::DAC::Instance *dacInstancePtr = getDACInstance_(dacName);

    if (dacInstancePtr) {
      // Update the time-voltage pairs for the given DAC instance.
      if (!dacInstancePtr->updateTVVEC(tv_pair_vector))
        Report::UserWarning0() << "Failed to update the time-voltage pairs for the DAC " << dacName;
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : getTimeVoltagePairs
// Purpose        : get a map of all time-voltage pairs from all ADC instances
//
// Special Notes  : The current method for locating ADC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          : public
// Creator        : Tom Russo
// Creation Date  : 05/10/2004
//----------------------------------------------------------------------------
bool Simulator::getTimeVoltagePairs(std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap)
{
  Device::Device *device = devIntPtr_->getDevice(Device::ADC::Traits::modelGroup());
  if (device) {
    TimeVoltagePairsOp op(timeVoltageUpdateMap);
    timeVoltageUpdateMap.clear();

    device->forEachInstance(op);
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : setADCWidths
// Purpose        : Update the ADC devices in a circuit by informing them
//                  of the width of their bitvector output on the
//                  "digital side"
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::setADCWidths(const std::map<std::string, int> &ADCWidthMap)
{
  for (std::map<std::string, int>::const_iterator it = ADCWidthMap.begin(), end = ADCWidthMap.end(); it != end; ++it) {
    const std::string &adcName = (*it).first;
    const int width = (*it).second;

    Device::ADC::Instance *adcInstancePtr = getADCInstance_(adcName);

    if (adcInstancePtr) {
      // Update the time-voltage pairs for the given ADC instance.
      if (!adcInstancePtr->setBitVectorWidth(width))
        Report::UserWarning0() << "Failed to update the width for ADC " << adcName;
    }
  }

  return true;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulateUntil
// Purpose       : To continue the existing analog circuit simulation
//                 until either the given <requestedUntilTime> is reached
//                 or the simulation termination criterion is met.
//                 Return a Boolean indicating whether the simulation
//                 run was successful. (Note that the run is successful
//                 even when the given <requestedUntilTime> is not reached,
//                 so long as the run completed normally.)
// Special Notes : The time variables are in units of seconds.
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool Simulator::simulateUntil(
  double        requestedUntilTime,
  double &      completedUntilTime)
{
  bool bsuccess = false;
  double currentTimeBeforeSim = analysisManager_->getTime();
  double finalTime=analysisManager_->getFinalTime();
  double initialTime = analysisManager_->getInitialTime();

  // Silence the "Percent Complete" noise, we don't want it when using this
  // interface

  analysisManager_->silenceProgress();

#ifdef Xyce_DEBUG_CIRCUIT
  dout() << "simulateUntil: ";
  dout() << "finalTime = " << finalTime
            << ", currentTimeBeforeSim = " << currentTimeBeforeSim << std::endl;
#endif

  if (currentTimeBeforeSim >= finalTime)
  {
    // We have already simulated as far as the netlist requested.
    bsuccess = true;
    completedUntilTime = currentTimeBeforeSim;
#ifdef Xyce_DEBUG_CIRCUIT
    dout() << "Case1: completedUntilTime = " << completedUntilTime;
#endif
  }
  else
  {
    // We are not already past the end of the netlist time
    analysisManager_->setPauseTime(Xycemin(requestedUntilTime,finalTime));
#ifdef Xyce_DEBUG_CIRCUIT
    dout() << "simulateUntil currentTimeBeforeSim = " << currentTimeBeforeSim << "  initialTime = " << initialTime << std::endl;
#endif
    if (currentTimeBeforeSim > initialTime)
    {
      analysisManager_->resumeSimulation();
    }

#ifdef Xyce_DEBUG_CIRCUIT
    dout() << "simulateUntil: Case2: requestedUntilTime = " << requestedUntilTime
                 << ", pauseTime = " << analysisManager_->getPauseTime() << std::endl;
#endif

    bsuccess = runSolvers_();
    completedUntilTime = analysisManager_->getTime();
#ifdef Xyce_DEBUG_CIRCUIT
    dout() << "simulateUntil: Case2: completedUntilTime = " << completedUntilTime << std::endl;
#endif
  }

#ifdef Xyce_DEBUG_CIRCUIT
  dout() << std::endl;
#endif

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::finalize
// Purpose       : To clean up after driving Xyce with the SIMBUS
//                 simulation backplane. This includes the following:
//                    Free any dynamically allocated memory...
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Lisa Maynes, CoMeT Solutions, Inc.
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool Simulator::finalize()
{
  bool bsuccess = true;

  Xyce::lout() << "\n***** Solution Summary *****"  << std::endl;

  analysisManager_->printLoopInfo(0, 0);

  // The Stats will replace this bit of ugliness.  But for now I'll write a function to just get the value.
  Analysis::StatCounts analysis_stat_counts = analysisManager_->getAnalysisObject().getStatCounts() - analysisManager_->getAnalysisObject().getStatCounts(0);

  IO::outputMacroResults(comm_,
                               outputManager_->getObjectiveMap(),
                               *measureManager_,
                               *fourierManager_,
                               outputManager_->getNetListFilename(),
                               outputManager_->getResponseFunctions(),
                               outputResponse_->getResponseFilename(),
                               outputManager_->getStepLoopNumber());

  rootStat_.stop();

  Xyce::lout() << std::endl
               << "***** Total Simulation Solvers Run Time: " << XyceTimerPtr_->elapsedTime() << " seconds" << std::endl
               << "***** Total Elapsed Run Time:            " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl
               << "*****" << std::endl
               << "***** End of Xyce(TM) Simulation" << std::endl
               << "*****" << std::endl;

  const char *xyce_no_tracking = ::getenv("XYCE_NO_TRACKING");
  if (Parallel::rank(comm_) == 0 && trackingURL && !xyce_no_tracking)
  {
    const time_t now=time(NULL);
    char timeDate[40];
    strftime(timeDate,40,"%x %X %Z",localtime(&now));

    auditJSON_ << Util::JSON::sep
               << Util::nameValuePair("Hostname", hostname()) << Util::JSON::sep
               << Util::nameValuePair("Domainname", domainname()) << Util::JSON::sep
               << Util::nameValuePair("Username", username()) << Util::JSON::sep
               << Util::nameValuePair("Hardware", hardware()) << Util::JSON::sep
               << Util::nameValuePair("OSname", osname()) << Util::JSON::sep
               << Util::nameValuePair("OSversion", osversion()) << Util::JSON::sep
               << Util::nameValuePair("Version", Util::Version::getFullVersionString()) << Util::JSON::sep
               << Util::nameValuePair("Processors", Parallel::size(comm_)) << Util::JSON::sep
               << Util::nameValuePair("PrimaryAnalysis", Analysis::analysisModeName(analysisManager_->getAnalysisMode())) << Util::JSON::sep
               << Util::nameValuePair("EndTime", timeDate) << Util::JSON::sep
               << Util::nameValuePair("RuntimeStats", analysis_stat_counts);
    Util::JSON json;
    json << Util::nameValuePair("audit", auditJSON_);

    const std::string &json_string = json.str();

    Util::sendTrackingData(trackingURL, 0, json_string);
  }

  if (Parallel::rank(comm_) == 0 && Parallel::size(comm_) > 1) {
    pout() << std::endl
           << "Timing summary of processor " << Parallel::rank(comm_) << std::endl;
    Stats::printStatsTable(pout(), rootStat_, Stats::METRICS_ALL, false);
  }

  Xyce::lout() << std::endl
               << "Timing summary of " << Parallel::size(comm_) << " processor" << (Parallel::size(comm_) == 1 ? "" : "s") << std::endl;
  Stats::printStatsTable(Xyce::lout(), rootStat_, Stats::METRICS_ALL, false, comm_);

  // Close the output stream:
  closeLogFile();

  bsuccess = doDeAllocations_();

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::reportTotalElapsedTime ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/01/2009
//---------------------------------------------------------------------------
void Simulator::reportTotalElapsedTime()
{
  Xyce::lout() <<  "\n***** Total Elapsed Run Time: " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulationComplete
// Purpose       : Simply report whether we've reached the end of the
//                 simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//---------------------------------------------------------------------------
bool Simulator::simulationComplete()
{
  return analysisManager_->isSimulationComplete();
}

//
// new mixed-signal functions:
// These are provisional!

//---------------------------------------------------------------------------
// Function      : Simulator::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
bool Simulator::provisionalStep
  (double maxTimeStep,
   double &timeStep,
   std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap)
{
  bool bsuccess=true;

  bool b1 = analysisManager_->provisionalStep(maxTimeStep, timeStep);
  bsuccess = bsuccess && b1;

  b1=getTimeVoltagePairs(timeVoltageUpdateMap);

  bsuccess = bsuccess && b1;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/25/2009
//---------------------------------------------------------------------------
double Simulator::getFinalTime()
{
  double ft=0.0;
  if (analysisManager_!=0)
  {
    ft = analysisManager_->getFinalTime();
  }
  return ft;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
double Simulator::getTime()
{
  double t=0.0;
  if (analysisManager_!=0)
  {
    t = analysisManager_->getTime();
  }
  return t;
}

//---------------------------------------------------------------------------
// Function      : Simulator::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
void Simulator::acceptProvisionalStep()
{
  analysisManager_->acceptProvisionalStep ();
}

//---------------------------------------------------------------------------
// Function      : Simulator::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
void Simulator::rejectProvisionalStep()
{
  analysisManager_->rejectProvisionalStep();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR Two-level Functions:

//---------------------------------------------------------------------------
// Function      : Simulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/2006
//---------------------------------------------------------------------------
bool Simulator::simulateStep
      ( const N_DEV_SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  dout() << "\nsimulateStep: " << std::endl;
#endif

  // Apply the input voltages to their appropriate sources:
  std::map<std::string,double>::const_iterator iterM = inputMap.begin();
  std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
  int i=0;
  for (; iterM != endM;++iterM, ++i)
  {
    bool found = true;
    std::string name = iterM->first;
    double val  = iterM->second;
    devIntPtr_->setParam (name,val);
  }

  analysisManager_->setExternalSolverState (solState);
  bsuccess = analysisManager_->runStep (solState.tiInfo, tlError);

  // calculate the conductance:
  nlsMgrPtr_->obtainConductances(
        inputMap,
        outputVector,
        jacobian
    );

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 07/16/2009
//---------------------------------------------------------------------------
bool Simulator::simulateStep
      (const N_DEV_ExternalSimulationData & ext_data,
       const std::map<std::string,double> & inputMap,
       std::vector<double> & outputVector,
       std::vector< std::vector<double> > & jacobian,
       N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  dout() << "\nsimulateStep: " << std::endl;
#endif

  // Create N_DEX_SolverState object
  N_DEV_SolverState state;

  if (ext_data.is_transient)
  {
    Xyce::lout() << "Xychron (mixed-model): Transient Xyce Step" << std::endl;

    state.currTimeStep = ext_data.current_time_step_size;
    state.lastTimeStep = ext_data.previous_time_step_size;
    state.currTime = ext_data.current_time;
    state.finalTime = ext_data.final_time;
    state.timeStepNumber = ext_data.time_step_number;
    state.startingTimeStep = ext_data.current_time_step_size;
    state.tranopFlag = false;
    state.transientFlag = true;
    state.dcopFlag = false;
    state.dcsweepFlag = false;
    // tiInfo
    state.tiInfo.nextTimeStep = ext_data.current_time_step_size;
    state.tiInfo.currTimeStep = ext_data.current_time_step_size;
    state.tiInfo.nextTime = ext_data.current_time;
    state.tiInfo.finalTime = ext_data.final_time;
    state.tiInfo.dcopFlag = false;
    state.tiInfo.tranopFlag = false;
    state.tiInfo.transientFlag = true;
    state.tiInfo.dcsweepFlag = false;
    state.tiInfo.timeStepNumber = ext_data.time_step_number;
    state.tiInfo.timeIntMode = 6;

    // New changes
    // pdt = 1/dt
    state.pdt = 1.0/ext_data.current_time_step_size;
    // bptol - ignore for now
    state.doubleDCOPEnabled = false;
    state.doubleDCOPStep = 1;
    // newtonIteration - doens't matter for particular devices.  For volt lim

    if (state.timeStepNumber == 0)
    {
      state.initTranFlag = true;
      state.tiInfo.initTranFlag = true;
      state.tiInfo.beginIntegrationFlag = true;
    }
    else
    {
      state.initTranFlag = false;
      state.tiInfo.initTranFlag = false;
      state.tiInfo.beginIntegrationFlag = false;
    }

    // tia pdt - set this
    state.tiInfo.pdt = state.pdt;
    state.tiInfo.currentOrder = 1;
  }
  else
  {
    state.tranopFlag = true;  // for a DCOP only, this should be false
    state.transientFlag = false;
    state.dcopFlag = false;
    state.dcsweepFlag = false;
    state.tiInfo.dcopFlag = true;
    state.tiInfo.tranopFlag = true;
    state.tiInfo.transientFlag = false;
    state.tiInfo.dcsweepFlag = false;
    state.tiInfo.timeIntMode = 0;
  }

  bsuccess = this->simulateStep(state,inputMap,outputVector,jacobian,tlError);

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startupSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool Simulator::startupSolvers ()
{
  bool bsuccess = true;
  bsuccess = analysisManager_->startupSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::finishSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool Simulator::finishSolvers ()
{
  bool bsuccess = true;
  bsuccess = analysisManager_->finishSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::homotopyStepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
void Simulator::homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals)
{
  analysisManager_->homotopyStepSuccess ( paramNames, paramVals);
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::homotopyStepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//---------------------------------------------------------------------------
void Simulator::homotopyStepFailure ()
{
  analysisManager_->homotopyStepFailure ();
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::stepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void Simulator::stepSuccess(Analysis::CurrentMode analysis)
{
  analysisManager_->stepSuccess(analysis);
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::stepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void Simulator::stepFailure(Analysis::CurrentMode analysis)
{
  analysisManager_->stepFailure(analysis);
  return;
}


//---------------------------------------------------------------------------
// Function      : Simulator::getInitialQnorm
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool Simulator::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = analysisManager_->getInitialQnorm (tle);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getBreakPoints
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool Simulator::getBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = analysisManager_->getBreakPoints (breakPointTimes);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::updateStateArrays
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool Simulator::updateStateArrays ()
{
  bool bsuccess = true;
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::setInternalParam
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool Simulator::setInternalParam (std::string & name, double val)
{
  devIntPtr_->setParam (name, val);
  return true;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
bool Simulator::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  bool bsuccess = analysisManager_->startTimeStep(tiInfo);
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 08/28/2009
//---------------------------------------------------------------------------
bool Simulator::startTimeStep (const N_DEV_ExternalSimulationData & ext_data)
{
  N_TIA_TimeIntInfo tiInfo;

  if (ext_data.is_transient)
  {
    tiInfo.nextTimeStep = ext_data.current_time_step_size;
    tiInfo.currTimeStep = ext_data.current_time_step_size;
    tiInfo.nextTime = ext_data.current_time;
    tiInfo.finalTime = ext_data.final_time;
    tiInfo.dcopFlag = false;
    tiInfo.tranopFlag = false;
    tiInfo.transientFlag = true;
    tiInfo.dcsweepFlag = false;
    tiInfo.timeStepNumber = ext_data.time_step_number;
    tiInfo.timeIntMode = 6;

    if (ext_data.time_step_number == 0)
    {
      tiInfo.initTranFlag = true;
      tiInfo.beginIntegrationFlag = true;
    }
    else
    {
      tiInfo.initTranFlag = false;
      tiInfo.beginIntegrationFlag = false;
    }

    // tia pdt - set this
    tiInfo.pdt = 1.0/ext_data.current_time_step_size;
    tiInfo.currentOrder = 1;
  }
  else
  {
    tiInfo.dcopFlag = true;
    tiInfo.tranopFlag = true;
    tiInfo.transientFlag = false;
    tiInfo.dcsweepFlag = false;
    tiInfo.timeIntMode = 0;
  }

  bool bsuccess = analysisManager_->startTimeStep(tiInfo);
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
bool Simulator::endTimeStep (N_DEV_ExternalSimulationData & ext_data)
{
  // We could opt to obtain selected time integration data here or
  // allow it to be obtained higher up, eg in charon::sc::CircuitDriver
  bool bsuccess = true;
  N_TIA_TimeIntInfo tiInfo;
  analysisManager_->getTimeIntInfo(tiInfo);

  // Copy N_TIA_TimeIntInfo into our ext_data container
  ext_data.currentOrder         = tiInfo.currentOrder        ;
  ext_data.numberOfSteps        = tiInfo.numberOfSteps       ;
  ext_data.usedOrder            = tiInfo.usedOrder           ;
  ext_data.nscsco               = tiInfo.nscsco              ;
  ext_data.pdt                  = tiInfo.pdt                 ;
  ext_data.nextTimeStep         = tiInfo.nextTimeStep        ;
  ext_data.currTimeStep         = tiInfo.currTimeStep        ;
  ext_data.currentTime          = tiInfo.currentTime         ;
  ext_data.nextTime             = tiInfo.nextTime            ;
  ext_data.finalTime            = tiInfo.finalTime           ;
  ext_data.startingTimeStep     = tiInfo.startingTimeStep    ;
  ext_data.bpTol                = tiInfo.bpTol               ;
  ext_data.dcopFlag             = tiInfo.dcopFlag            ;
  ext_data.acopFlag             = tiInfo.acopFlag            ;
  ext_data.inputOPFlag          = tiInfo.inputOPFlag         ;
  ext_data.tranopFlag           = tiInfo.tranopFlag          ;
  ext_data.transientFlag        = tiInfo.transientFlag       ;
  ext_data.dcsweepFlag          = tiInfo.dcsweepFlag         ;
  ext_data.timeStepNumber       = tiInfo.timeStepNumber      ;
  ext_data.initTranFlag         = tiInfo.initTranFlag        ;
  ext_data.beginIntegrationFlag = tiInfo.beginIntegrationFlag;
  ext_data.doubleDCOPStep       = tiInfo.doubleDCOPStep      ;
  ext_data.doubleDCOPEnabled    = tiInfo.doubleDCOPEnabled   ;
  ext_data.stepLoopIter         = tiInfo.stepLoopIter        ;
  ext_data.timeIntMode          = tiInfo.timeIntMode         ;
  ext_data.sweepSourceResetFlag = tiInfo.sweepSourceResetFlag;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
void Simulator::enable_lte_analysis()
{
  N_TIA_TIAParams & tiaParams = analysisManager_->getTIAParams();
  tiaParams.errorAnalysisOption = 0; // use local truncation error estimates
}

//---------------------------------------------------------------------------
// Function      : Simulator::readExternalParamsFromFile
// Purpose       :
// Special Notes : Used to read parameters "tag" = "value" from an external
//                 file.  Any "tag"s found while reading the netlist during
//                 parsing will be replaced by "value".
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical Systems Modeling
// Creation Date : 07/24/2012
//---------------------------------------------------------------------------
void Simulator::readExternalParamsFromFile( std::string filename,
  std::vector< std::pair< std::string, std::string > > & paramList )
{
  // at this stage just support the Dakota params.in format of "value" = "tag".
  // we could support other formats as well.

  const std::string allowedChars("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_.$");
  const std::string whiteSpace(" \t=\n\r");    // note we will treat "=" as whitespace
  const std::string commentChars("*#;");  // if we find any of these then the rest of the line is a comment

  // attempt to open the params file
  std::ifstream paramFile(filename.c_str(), std::ios::in);
  if( paramFile )
  {
    // loop over the file trying to gather names / values and response functions required.
    // general format:
    // white space tag/value white space or = tag/value
    //
    std::string aLine;
    getline(paramFile, aLine);
    while ( paramFile.good() )
    {
      // getline worked, so try to parse the line
      // procedure (1) discard any comments
      //           (2) find tag/value boundaries by looking for allowedChars followed by whitespace
      //           (3) sort out order of tag/value pair
      //           (4) store pair in paramList
      //           (5) try to get another line

      // don't bother with lines that are blank or essentially blank (i.e. just \r\n)
      // shortest line we could have is x=1 or 3 chars.
      if( aLine.length() > 2 )
      {
        std::string::size_type commentLoc = aLine.find_first_of( commentChars, 0 );
        if( commentLoc != std::string::npos )
        {
          // chop off comment to end of line.
          aLine.erase( commentLoc, aLine.length()-commentLoc );
        }
        //check overall line length again.  This could have been just a comment line
        if( aLine.length() > 2 )
        {
          std::string::size_type word1Start = aLine.find_first_of( allowedChars, 0 );
          // check that we found a valid word1Start otherwise stop trying to parse this line.
          if( word1Start != std::string::npos )
          {
            std::string::size_type word1End   = aLine.find_first_of( whiteSpace, word1Start );
            // check that we found a valid word1End otherwise stop trying to parse this line.
            if( word1End != std::string::npos )
            {
              std::string::size_type word2Start = aLine.find_first_of( allowedChars, word1End );
              // check that we found a valid word2Start otherwise stop trying to parse this line.
              if( word2Start != std::string::npos )
              {
                std::string::size_type word2End   = aLine.find_first_of( whiteSpace, word2Start );
                // check that we found a valid word2End
                if( word2End == std::string::npos )
                  word2End = aLine.length();
                // if we get here then we have valid start,end indicies for word1 and word2

                std::string word1=aLine.substr( word1Start, (word1End - word1Start) );
                std::string word2=aLine.substr( word2Start, (word2End - word2Start) );

                // sort out tag/value ordering.
                // if word1=number assume format is value = tag
                // otherwise assume format is tag = value
                std::stringstream converter;
                converter << word1;
                double testvalue;
                converter >> testvalue;
                if( converter.fail() )
                {
                  // couldn't convert tag1 to a double value so assume format is word1=tag word2=value
                  paramList.push_back( std::pair<std::string,std::string>(word1,word2) );
                }
                else
                {
                  // tag1 was successfully converted to a value so assume format is word1=value word2=tag
                  paramList.push_back( std::pair<std::string,std::string>(word2,word1) );
                }

              }  // if ( word2Start != std::string::npos )
            }    // if( word1End != std::string::npos )
          }      // if( word1Start != std::string::npos )
        }        // if( aLine.length() > 2 ) -- check after removing comments
      }          // if( aLine.lenght() > 2 ) -- outer check
      // try to get another line
      getline(paramFile, aLine);
    }
  }
  else
  {
    Report::UserWarning() << "Could not open parameter file: " + filename + ". Attempting to continue.";
  }

  // for debug purposes.  output the params as read
  dout() << "Parameters read from \"" << filename << "\"" << std::endl;
  std::vector< std::pair< std::string, std::string > >::iterator listitr = paramList.begin();
  std::vector< std::pair< std::string, std::string > >::iterator enditr = paramList.end();
  while( listitr != enditr )
  {
    dout() << "  " << listitr->first << " , " << listitr->second << std::endl;
    listitr++;
  }
}

//---------------------------------------------------------------------------
// Function      : Simulator::checkResponseVars
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what .measure lines are to be used as response functions
//                 to pass back to Dakota.  This call checks the measure manager
//                 in the I/O package has set up measure objects for each label.
// Scope         : public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/24/2014
//---------------------------------------------------------------------------
bool
Simulator::checkResponseVars(
  const std::vector<std::string> &      resStrings) const
{
  if (!outputManager_)
  {
    Report::DevelFatal0() << "checkResponseVars outputManager_ is null";
  }

  for (std::vector< std::string >::const_iterator it = resStrings.begin(), end = resStrings.end(); it != end; ++it)
  {
    if (!measureManager_->find(*it))
    {
      return false;
      break;
    }
  }

  return true;
}

//---------------------------------------------------------------------------
// Function      : Simulator::obtainResponses
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what .measure lines are to be used as response functions
//                 to pass back to Dakota.  This call obtains the responses for
//                 each labelled measure.
// Scope         : public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/24/2014
//---------------------------------------------------------------------------
void
Simulator::obtainResponses(
  const std::vector<std::string> &      resStrings,
  std::vector<double> &                 resVector) const
{
  if (!outputManager_)
  {
    Report::DevelFatal0() << "obtainResponses outputManager_ is null";
  }

  resVector.resize(resStrings.size());

  std::vector<double>::iterator d_it = resVector.begin();
  for (std::vector<std::string>::const_iterator it = resStrings.begin(), end = resStrings.end(); it != end; ++it, ++d_it)
  {
    bool found = true;
    measureManager_->getMeasureValue(*it, *d_it, found);

    if (!found)
      Report::UserError() << "Response " << *it << " was not found by .MEASURE";
  }
}

//---------------------------------------------------------------------------
// Function      : Simulator::initializeTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void Simulator::initializeTransientModel()
{
  analysisManager_->initializeTransientModel();
}

//---------------------------------------------------------------------------
// Function      : Simulator::evalTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
bool Simulator::evalTransientModel
    (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * CurrSolVectorPtr,
     N_LAS_Vector * LasSolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * CurrStaVectorPtr,
     N_LAS_Vector * LasStaVectorPtr,
     N_LAS_Vector * StaDerivVectorPtr,
     N_LAS_Vector * StoVectorPtr,
     N_LAS_Vector * CurrStoVectorPtr,
     N_LAS_Vector * LasStoVectorPtr,
     N_LAS_Vector * stoLeadCurrQCompVectorPtr,
     N_LAS_Vector * QVectorPtr,
     N_LAS_Vector * FVectorPtr,
     N_LAS_Vector * BVectorPtr,
     N_LAS_Vector * dFdxdVpVectorPtr,
     N_LAS_Vector * dQdxdVpVectorPtr,
     N_LAS_Matrix * dQdxMatrixPtr,
     N_LAS_Matrix * dFdxMatrixPtr
    )
{
  return (
      analysisManager_->evalTransientModel(
      t,
      SolVectorPtr,
      CurrSolVectorPtr,
      LasSolVectorPtr,

      StaVectorPtr,
      CurrStaVectorPtr,
      LasStaVectorPtr,
      StaDerivVectorPtr,
      StoVectorPtr,
      CurrStoVectorPtr,
      LasStoVectorPtr,
      stoLeadCurrQCompVectorPtr,
      QVectorPtr,
      FVectorPtr,
      BVectorPtr,
      dFdxdVpVectorPtr,
      dQdxdVpVectorPtr,
      dQdxMatrixPtr,
      dFdxMatrixPtr
      )
    );
}

//---------------------------------------------------------------------------
// Function      : Simulator::evalTransientModelState
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//---------------------------------------------------------------------------
bool Simulator::evalTransientModelState
    (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * StoVectorPtr
    )
{
  return (
      analysisManager_->evalTransientModelState(
      t,
      SolVectorPtr,
      StaVectorPtr,
      StoVectorPtr
      )
    );
}

//---------------------------------------------------------------------------
// Function      : Simulator::getMapsAndGraphs
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void Simulator::getMapsAndGraphs
  (
   RCP<N_PDS_ParMap> & x_map,
   RCP<N_PDS_ParMap> & x_map_ognd,
   RCP<N_PDS_ParMap> & s_map,
   RCP<N_PDS_ParMap> & store_map,
   RCP<Epetra_CrsGraph> & dQdx_graph,
   RCP<Epetra_CrsGraph> & dQdx_graph_ognd,
   RCP<Epetra_CrsGraph> & dFdx_graph,
   RCP<Epetra_CrsGraph> & dFdx_graph_ognd
  )
{
  lasBuilderPtr_->getSolutionMaps(x_map,x_map_ognd);
  lasBuilderPtr_->getStateMap(s_map);
  lasBuilderPtr_->getStoreMap(store_map);
  lasBuilderPtr_->getdQdxGraphs(dQdx_graph,dQdx_graph_ognd);
  lasBuilderPtr_->getdFdxGraphs(dFdx_graph,dFdx_graph_ognd);
}


//---------------------------------------------------------------------------
// Function      : Simulator::getVariableNames
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey
// Creation Date : 07/28/09
//---------------------------------------------------------------------------
std::vector<std::string> Simulator::getVariableNames()
{
  return outputManager_->getVariableNames();
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getXygraInstancePtr_
// Purpose       : Returns the pointer to a named Xygra device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Device::DAC::Instance *
Simulator::getDACInstance_(const std::string &deviceName)
{
  // See if we've looked this up before.
  if (dacDeviceMap_.empty())
  {
    Device::Device *device = devIntPtr_->getDevice(Device::DAC::Traits::modelGroup());
    if (device)
      Device::mapDeviceInstances(*device, dacDeviceMap_);
  }

  std::map<std::string, Device::DAC::Instance *>::iterator mapIter = dacDeviceMap_.find(deviceName);
  if (mapIter == dacDeviceMap_.end())
    return 0;

  return (*mapIter).second;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getXygraInstancePtr_
// Purpose       : Returns the pointer to a named Xygra device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Device::ADC::Instance *
Simulator::getADCInstance_(const std::string &deviceName)
{
  // See if we've looked this up before.
  if (adcDeviceMap_.empty())
  {
    Device::Device *device = devIntPtr_->getDevice(Device::ADC::Traits::modelGroup());
    if (device)
      Device::mapDeviceInstances(*device, adcDeviceMap_);
  }

  std::map<std::string, Device::ADC::Instance *>::iterator mapIter = adcDeviceMap_.find(deviceName);
  if (mapIter == adcDeviceMap_.end())
    return 0;

  return (*mapIter).second;
}



void
report_handler(
  const char *  message,
  unsigned      report_mask)
{
  // if ( !comm_->isSerial() && !(report_mask & MSG_SYMMETRIC))
  //   os << "P" << comm_->procID() << " - ";

  std::ostringstream oss;
  Util::word_wrap(oss, message, 78, " ", "");

  // If symmetric then all processors are getting the same message, only write to p0 and not to ~p0 backlog.  If
  // asymetric then one processor is getting the message, write to per processor stream which writes to per processor
  // log file and to backlog.
  if (report_mask & Report::MSG_SYMMETRIC)
    Xyce::lout() << oss.str();
  else
    pout() << oss.str();

  // If fatal error also send the message to the standard error file:
  // Also save it for output on proc 0 if running in parallel
  if (report_mask & Report::MSG_TERMINATE)
  {
    std::cerr << oss.str() << std::endl;
    Report::abort();
  }
}

} // namespace Circuit
} // namespace Xyce
