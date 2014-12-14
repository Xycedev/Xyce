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
// Filename       : $RCSfile: N_CIR_Xyce.h,v $
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
// Revision Number: $Revision: 1.118.2.8 $
//
// Revision Date  : $Date: 2014/08/29 21:26:33 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_CIRCUIT_h
#define Xyce_N_CIR_CIRCUIT_h

#include <string>
#include <map>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_ERH_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_NLS_fwd.h>

#include <N_PDS_Manager.h>
#include <N_UTL_ReportHandler.h>
#include <N_UTL_Stats.h>
#include <N_UTL_JSON.h>
#include <N_DEV_ADC.h>
#include <N_DEV_DAC.h>

class N_LAS_System;
class N_LAS_Builder;
class N_LAS_Vector;
class N_LAS_Matrix;

class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

class N_LOA_Loader;
class N_LOA_CktLoader;
class N_LOA_NonlinearEquationLoader;
class N_LOA_LoaderMgr;

class Epetra_CrsGraph;

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

#include <N_PDS_fwd.h>
#include <N_IO_CmdParse.h>

namespace Xyce {
namespace Circuit {

//-----------------------------------------------------------------------------
// Class         : Xyce
// Purpose       : This is the main "top level" class for Xyce.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class Simulator
{
  // Methods
 public:

  // Default constructor
    Simulator(Parallel::Machine comm = 0);

  // Default destructor (virtual so we derive other classes from this properly)
  virtual ~Simulator();

  // These are all the API calls that we are suppose to be making available
  // for external programs and/or other objects

  //---------------------------------------------------------------------------
  // Function      : setNetlistParameters
  // Purpose       : This passes a vector of pairs "key", "value" that will
  //                 be substituted during the processing of the netlist.  This
  //                 more easily allows Dakota to change any netlist parameter
  //                 during netlist setup.
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical and MEMS Modeling
  // Creation Date : 10/9/2008
  //---------------------------------------------------------------------------
  void setNetlistParameters( const std::vector< std::pair< std::string, std::string > > & externalParams );


  //---------------------------------------------------------------------------
  // Function      : setNetlistParameters
  // Purpose       : Call through to the output manager to set the suffix to
  //                 be used on the output file, as in circuit + suffix + prn
  //                 This is useful in Dakota controlled runs to keep each
  //                 simulation from overwritting the last one.
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical and MEMS Modeling
  // Creation Date : 10/9/2008
  //---------------------------------------------------------------------------
  void setOutputFileSuffix( const std::string newSuffix );

  //---------------------------------------------------------------------------
  // Function      : run
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 02/19/01
  //---------------------------------------------------------------------------
  bool run(int argc, char *argv[]);

  //---------------------------------------------------------------------------
  // Function      : initialize
  // Purpose       : To initialize Xyce to be driven by the SAX
  //                 simulation backplane. This includes the following:
  //                    Set up and register the Parallel Manager.
  //                    Parse the command-line arguments.
  //                    Redirect the output stream of processor 1,
  //                    if requested.
  //                    Read in the Netlist.
  //                    Allocate and register the external packages.
  //                    Set up the representation of the circuit topology.
  //                    Set up the matrix structures.
  //                    Initialize the solvers.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  bool initialize(int argc, char **argv);

  //---------------------------------------------------------------------------
  // Function      : runSimulation
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 5/27/00
  //---------------------------------------------------------------------------
  bool runSimulation();

  bool getDeviceNames(const std::string & modelGroupName, std::vector<std::string> & deviceNames);

  //---------------------------------------------------------------------------
  // Function      : getDACDeviceNames
  // Purpose       : Gets the (stripped) names of the DAC devices
  //                 in the circuit.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/13/03
  //---------------------------------------------------------------------------
  bool getDACDeviceNames(std::vector< std::string >& dacNames);

  //---------------------------------------------------------------------------
  // Function      : getADCMap
  // Purpose       : Gets the (stripped) names of the ADC devices
  //                 in the circuit(as key of map) and map of parameters
  //                 (keyed by parameter name) for each device
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo, SNL, Component Information and Models
  // Creation Date : 05/06/2004
  //---------------------------------------------------------------------------
  bool getADCMap(std::map<std::string,std::map<std::string,double> >& ADCMap);

  //---------------------------------------------------------------------------
  // Function      : updateTimeVoltagePairs
  // Purpose       : Update the DAC devices in a circuit by adding the set
  //                 of time and voltage pairs built up on the "digital side"
  //                 since the last update and by removing the time-voltage
  //                 pairs for times that pre-date the given simulation time.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/10/03
  //---------------------------------------------------------------------------
  bool updateTimeVoltagePairs(
        std::map< std::string, std::vector< std::pair<double,double> >* > const&
        timeVoltageUpdateMap);

  //---------------------------------------------------------------------------
  // Function      : getTimeVoltagePairs
  // Purpose       : query the DAC devices in a circuit for the set
  //                 of time and voltage pairs
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo, SNL ComponentInformation and Models
  // Creation Date : 05/10/2004
  //---------------------------------------------------------------------------
  bool getTimeVoltagePairs(
        std::map< std::string, std::vector< std::pair<double,double> > > &
        timeVoltageUpdateMap);

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
  bool setADCWidths(std::map< std::string, int > const& ADCWidthMap);

  //---------------------------------------------------------------------------
  // Function      : simulateUntil
  // Purpose       : To continue the existing analog circuit simulation
  //                 until either the given <requestedUntilTime> is reached
  //                 or the simulation termination criterion is met.
  //                 Return a Boolean indicating whether the simulation
  //                 run was successful. (Note that the run is successful
  //                 even when the given <requestedUntilTime> is not reached,
  //                 so long as the run completed normally.)
  // Special Notes : The time variables are in units of seconds.
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  bool simulateUntil(double requestedUntilTime, double& completedUntilTime);

  //---------------------------------------------------------------------------
  // Function      : finalize
  // Purpose       : To clean up after driving Xyce with the SIMBUS
  //                 simulation backplane. This includes the following:
  //                    Free any dynamically allocated memory...
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/29/03
  //---------------------------------------------------------------------------
  bool finalize();

  void reportTotalElapsedTime ();

  void readExternalParamsFromFile( std::string filename, std::vector< std::pair< std::string, std::string > > & paramList );

  bool checkResponseVars(const std::vector< std::string >& resStrings) const;
  void obtainResponses(const std::vector< std::string >& resStrings, std::vector< double >& resVector) const;

  // report on whether simulation is finished or not
  bool simulationComplete();

  // new mixed-signal functions:
  bool provisionalStep
    (double maxtimeStep,
     double &timeStep,
     std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap);

  void acceptProvisionalStep();
  void rejectProvisionalStep();

  double getFinalTime();
  double getTime();

  // 2-level, power node functions:
  bool startupSolvers ();
  bool finishSolvers ();

  void homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess(Analysis::CurrentMode analysis);
  void stepFailure(Analysis::CurrentMode analysis);
  bool getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);
  bool updateStateArrays ();
  bool startTimeStep (const N_TIA_TimeIntInfo & tiInfo);
  bool startTimeStep (const N_DEV_ExternalSimulationData & ext_data);
  bool endTimeStep (N_DEV_ExternalSimulationData & ext_data);
  void enable_lte_analysis();
  bool setInternalParam (std::string & name, double val);

  bool getInitialQnorm (N_TIA_TwoLevelError & tle);

  bool simulateStep
      ( const N_DEV_SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError);

  // 05/27/09 ModelEvaluator Interface functions
  void initializeTransientModel();

  bool evalTransientModel (
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
    );

  bool evalTransientModelState (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * StoVectorPtr
    );

  void getMapsAndGraphs (
    RCP<N_PDS_ParMap> & x_map,
    RCP<N_PDS_ParMap> & x_map_ognd,
    RCP<N_PDS_ParMap> & s_map,
    RCP<N_PDS_ParMap> & store_map,
    RCP<Epetra_CrsGraph> & dFdx_graph,
    RCP<Epetra_CrsGraph> & dFdx_graph_ognd,
    RCP<Epetra_CrsGraph> & dQdx_graph,
    RCP<Epetra_CrsGraph> & dQdx_graph_ognd
    );

  std::vector<std::string> getVariableNames();

  bool simulateStep
      ( const N_DEV_ExternalSimulationData & ext_data,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError);

  private:

    bool setupParMgr_(int iargs, char **cargs);
    bool doAllocations_();
    bool doInitializations_();
    bool doRegistrations_();
    bool doDeAllocations_();
    bool setUpTopology_();
    bool setUpMatrixStructure_();
    bool runSolvers_();

    Device::ADC::Instance *getADCInstance_(const std::string &deviceName);
    Device::DAC::Instance *getDACInstance_(const std::string &deviceName);

  private:
    Parallel::Machine                   comm_;

  protected:
    // We put the device interface pointer as a protected method so we can
    // access it from derived classes without kludges.

    // Device manager
    Device::DeviceInterface *     devIntPtr_;

  private:
    std::map<std::string, Device::DAC::Instance *> dacDeviceMap_;
    std::map<std::string, Device::ADC::Instance *> adcDeviceMap_;

    Topo::Manager *                     topMgrPtr_;                     ///< Topology manager
    Topo::Topology *                    topPtr_;
    N_LAS_System *                      lasSysPtr_;                     ///< Linear algebra system
    N_LAS_Builder *                     lasBuilderPtr_;                 ///< Linear algebra system builder
    Analysis::AnalysisManager *         analysisManager_;               ///< Analysis manager
    N_NLS_Manager *                     nlsMgrPtr_;                     ///< Time integration manager
    N_PDS_Manager *                     parMgrPtr_;                     ///< Parallel distribution manager
    IO::NetlistImportTool *             netlistImportToolPtr_;          ///< Netlist import tool.
    Util::Op::BuilderManager *          opBuilderManager_;              ///< Op builder manager
    IO::OutputMgr *                     outputManager_;                 ///< Output manager
    IO::Measure::Manager *              measureManager_;                ///< Measure manager
    IO::FourierMgr *                    fourierManager_;                ///< Fourier manager
    IO::OutputResponse *                outputResponse_;                ///< Response output
    IO::RestartMgr *                    resMgrPtr_;                     ///< Restart manager
    Stats::Stat                         rootStat_;                      ///< Stats collection root
    Util::JSON                          auditJSON_;                     ///< Audit JSON structure
    Util::Timer *                       XyceTimerPtr_;                  ///< Xyce solver timing utility
    Util::Timer *                       ElapsedTimerPtr_;               ///< Elapsed time from beginning of run
    IO::PkgOptionsMgr *                 pkgOptionsMgrPtr_;              ///< Package options manager

    int                                 argc_;
    char **                             argv_;

    // if the user is providing an external file with parameters in it
    std::vector< std::pair< std::string, std::string> > externalNetlistParams_;

    bool isSerialFlag_;
    bool procZeroFlag_;
    bool multiThreading_;
    int numThreads_;

    bool initializeAllFlag_;

    IO::CmdParse commandLine;

    REH previousReportHandler_;
};

} // namespace Circuit
} // namespace Xyce

typedef Xyce::Circuit::Simulator N_CIR_Xyce;

#endif // Xyce_N_CIR_CIRCUIT_h
