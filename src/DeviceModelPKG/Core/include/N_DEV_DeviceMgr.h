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
// Filename       : $RCSfile: N_DEV_DeviceMgr.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.374.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:49 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceMgr_h
#define Xyce_N_DEV_DeviceMgr_h

#include <map>
#include <set>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_NLS_fwd.h>

#include <N_DEV_ArtificialParameters.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>

class N_LAS_Matrix;
class N_LAS_System;
class N_LAS_Vector;
class N_TIA_TwoLevelError;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceMgr : public Util::Listener<Analysis::StepEvent>
{
  friend class ArtificialParameters::ArtificialParameter;

public:
  typedef std::vector<Device *> DeviceVector;
  typedef std::vector<DeviceEntity *> EntityVector;
  typedef std::vector<DeviceInstance *> InstanceVector;
  typedef std::vector<DeviceModel *> ModelVector;
  typedef std::map<ModelTypeId, ModelVector> ModelTypeModelVectorMap;
  typedef std::map<ModelTypeId, InstanceVector> ModelTypeInstanceVectorMap;
  typedef std::map<std::string, DeviceEntity *, LessNoCase> DeviceEntityMap;
  typedef std::map<std::string, ModelTypeId, LessNoCase> ModelTypeNameModelTypeIdMap;
  typedef std::map<std::string, ArtificialParameters::ArtificialParameter *, LessNoCase> ArtificialParameterMap;

  static DeviceMgr * factory(IO::CmdParse & cp);

private:
  DeviceMgr(IO::CmdParse & cp);                       ///< Only the factory can create a device manager

  DeviceMgr(const DeviceMgr &);                       ///< No copying
  DeviceMgr &operator=(const DeviceMgr &);            ///< No assignment

public:
  ~DeviceMgr();

  void notify(const Analysis::StepEvent &event);

  // registration functions:
  bool registerLinearSystem(N_LAS_System * tmp_system_ptr);
  bool registerAnalysisManager(N_ANP_AnalysisManager * tmp_anaIntPtr);
  bool registerOutputMgr(IO::OutputMgr * output_manager);
  bool registerMeasureMgr(IO::Measure::Manager * measure_manager);
  bool registerParallelMgr(N_PDS_Manager * tmp_pdsMgrPtr);
  bool registerNonlinearSolver (Nonlinear::Manager * tmp_nlsMgrPtr);
  bool registerICLoads( std::vector<std::pair<int,double> > * icLoads );
  bool registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr );

  // this function is called from the output manager (through the
  // device interface) to inform the device package of the devices for
  // which lead currents have been requested. The device manager will
  // take care of doing isolated F and Q loads for these devices so
  // the lead currents can be calculated
  bool setLeadCurrentRequests(const std::set<std::string> & deviceNames );

  // MPDE related registrations:
  std::vector<double> getFastSourcePeriod (std::vector<std::string>& sourceNames);
  std::vector<double> registerFastSources (std::vector<std::string> & sourceNames);
  void deRegisterFastSources (std::vector<std::string> & sourceNames);
  void deactivateSlowSources();
  void activateSlowSources();

  void setMPDEFlag( bool flagVal );
  void setBlockAnalysisFlag( bool flagVal );
  void setFastTime( double timeVal );

  // Initialization function, to be called after all registrations are
  // finished, and the linear system class is completely set up.
  bool initializeAll();

  // Device accessor functions:

  bool addDeviceModel(const ModelBlock & MB);

  std::pair<ModelTypeId, ModelTypeId> getModelType(const InstanceBlock &instance_block);

  bool verifyDeviceInstance(InstanceBlock & IB);

  DeviceInstance * addDeviceInstance(InstanceBlock & IB);

  bool deleteDeviceInstance (const std::string & name);

  int getHomotopyBlockSize() const;

  bool output ();
  bool finishOutput ();

  void  dotOpOutput ();

  // Load functions:
  bool setInitialGuess (N_LAS_Vector * solVectorPtr);
  bool loadDeviceMask();

  void debugOutput1();
  void debugOutput2();

  void getAnalyticSensitivities(
      std::string & name, 
      std::vector<double> & dfdpVec, 
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec);

  bool analyticSensitivitiesAvailable (std::string & name);
  bool setParam(std::string & name, double val);
  double getParamAndReduce(const std::string & name) const;
  bool getParamAndReduce(const std::string & name, double & val) const;
  double getParamNoReduce(const std::string & name) const;

  bool findParam(const std::string & name) const;

  bool updateTemperature(double val);

  bool updateSources();

  bool resetRHSLoadFlags (int index);

  const Nonlinear::Manager &getNlsMgrPtr() const 
  {
    return *nlsMgrPtr_;
  }

    bool getLinearSystemFlag ()
    {
      return linearSystemFlag_;
    }

    bool getVoltageLimiterFlag ()
    {
      return devOptions_.voltageLimiterFlag;
    }

    bool getPDESystemFlag ()
    {
      return solState_.PDESystemFlag;
    }

  // setup initial conditions on devices
  bool setICs (N_LAS_Vector * tmpSolVectorPtr,
               N_LAS_Vector * tmpCurrSolVectorPtr,
               N_LAS_Vector * tmpLastSolVectorPtr,
               N_LAS_Vector * tmpStaVectorPtr,
               N_LAS_Vector * tmpCurrStaVectorPtr,
               N_LAS_Vector * tmpLasStaVectorPtr,
               N_LAS_Vector * tmpStaDerivVectorPtr,
               N_LAS_Vector * tmpStoVectorPtr,
               N_LAS_Vector * tmpCurrStoVectorPtr,
               N_LAS_Vector * tmpLastStoVectorPtr,
               N_LAS_Vector * tmpQVectorPtr,
               N_LAS_Vector * tmpFVectorPtr,
               N_LAS_Vector * tmpBVectorPtr,
               N_LAS_Vector * tmpdFdxdVpVectorPtr,
               N_LAS_Vector * tmpdQdxdVpVectorPtr);

  // time integration stuff:
  bool   getBreakPoints     ( std::vector<Util::BreakPoint> & breakPointTimes );
  double getMaxTimeStepSize ();
  void  declareCurrentStepAsBreakpoint();

  // two-level newton and pde-continuation
  int  enablePDEContinuation ();
  bool disablePDEContinuation ();
  void getNumInterfaceNodes (std::vector<int> & numInterfaceNodes);
  bool loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iElectrode, const N_LAS_Vector * dxdvPtr);
  bool raiseDebugLevel (int increment);

  bool calcPDESubProblemInfo ();

  // load functions:
  bool loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                         N_LAS_Vector * tmpStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpStoVectorPtr,
                         N_LAS_Matrix * tmpdQdxMatrixPtr,
                         N_LAS_Matrix * tmpdFdxMatrixPtr);

  bool loadDAEVectors   (N_LAS_Vector * tmpNextSolVectorPtr,
                         N_LAS_Vector * tmpCurrSolVectorPtr,
                         N_LAS_Vector * tmpLastSolVectorPtr,
                         N_LAS_Vector * tmpNextStaVectorPtr,
                         N_LAS_Vector * tmpCurrStaVectorPtr,
                         N_LAS_Vector * tmpLastStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpNextStoVectorPtr,
                         N_LAS_Vector * tmpCurrStoVectorPtr,
                         N_LAS_Vector * tmpLastStoVectorPtr,
                         N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
                         N_LAS_Vector * tmpQVectorPtr,
                         N_LAS_Vector * tmpFVectorPtr,
                         N_LAS_Vector * tmpBVectorPtr,
                         N_LAS_Vector * tmpdFdxdVpVectorPtr,
                         N_LAS_Vector * tmpdQdxdVpVectorPtr);

  bool updateState      (N_LAS_Vector * nextSolVectorPtr,
                         N_LAS_Vector * currSolVectorPtr,
                         N_LAS_Vector * lastSolVectorPtr,
                         N_LAS_Vector * nextStaVectorPtr,
                         N_LAS_Vector * currStaVectorPtr,
                         N_LAS_Vector * lastStaVectorPtr,
                         N_LAS_Vector * nextStoVectorPtr,
                         N_LAS_Vector * currStoVectorPtr,
                         N_LAS_Vector * lastStoVectorPtr
                         );

  bool loadBVectorsforAC (N_LAS_Vector * bVecRealPtr,
                          N_LAS_Vector * bVecImagPtr);

  bool getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec, std::vector<int>& bMatPosEntriesVec);

  // voltlim doesn't work with MPDE, but does work for DCOP and the
  // various initial conditions used in MPDE.  Hence the need for these
  // functions.
  //    void setVoltageLimiterFlag ();
  void unsetVoltageLimiterFlag ();

  void setVoltageLimiterFlag ( bool flagVal );

  void addGlobalPar(const Util::Param &);
  const double *findGlobalPar( const std::string & parName) const;
  double getGlobalPar( const std::string & parName ) const;


  // functions related to options registration
  bool registerOptions(const Util::OptionBlock & OB) {
    bool result = devOptions_.registerOptions(OB);
    devOptions_.applyCmdLineOptions(commandLine_);
    return result;
  }
  bool registerSensParams(const Util::OptionBlock & OB);
  bool registerTimeOptions(const Util::OptionBlock & OB){return true;}
  bool setTranAnalysisParams(const Util::OptionBlock & OB);
  bool setDCAnalysisParams(const Util::OptionBlock & OB);
  bool setOPAnalysisParams(const Util::OptionBlock & OB);
  bool setSTEPAnalysisParams(const Util::OptionBlock & OB);
  bool setMPDEAnalysisParams(const Util::OptionBlock & OB);
  bool setHBAnalysisParams(const Util::OptionBlock & OB);
  bool setACAnalysisParams(const Util::OptionBlock & OB);
  bool setMORAnalysisParams(const Util::OptionBlock & OB);

  // convergence:  allow devices to signal back to the solvers that
  // they've played some game that invalidates normal convergence tests,
  // and so the solution should be considered unconverged no matter how
  // small the various norms are.
  bool allDevsConverged();

  // Similar to allDevsConverged, but specific to "inner" devices of
  // 2-level solves.  They are handled slightly differently.
  bool innerDevsConverged();

  // Functions needed for power node (2-level) algorithm):

  // for the parallel case, we need to give all the processors a copy
  // of the device so all the parallel synchronized calls such are
  // called by all processors together.
  bool setupExternalDevices();

  const std::map<std::string,int> &getDeviceCountMap() {
    return localDeviceCountMap_;
  }

  void addDeviceToCount(const std::string & device_name, int num_devs = 1)
  {
    localDeviceCountMap_[device_name] += num_devs;
  }
 
  void addDevicesToCount(const std::map<std::string,int> & device_map);

  DeviceEntity *getDeviceEntity(const std::string &full_param_name) const;

  void homotopyStepSuccess(const std::vector<std::string> & paramNames, const std::vector<double> & paramVals);
  void homotopyStepFailure ();

  void stepSuccess(Analysis::CurrentMode analysis);
  void stepFailure(Analysis::CurrentMode analysis);

  void acceptStep();

  bool getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec );
  bool getInnerLoopErrorSums (std::vector<N_TIA_TwoLevelError> & tleVec);

  bool updateStateArrays();
  bool startTimeStep ();
  void setExternalSolverState (const SolverState & ss);

  int restartDataSize(bool pack);

  // Output restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  // Load restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  // needed for parallel only:
  void setGlobalFlags ();

  bool setupSolverInfo() 
  {
    return setupSolverInfo_();
  }

  Device *getDevice(EntityTypeId model_type_id) 
  {
    EntityTypeIdDeviceMap::iterator it = deviceMap_.find(model_type_id);
    return it == deviceMap_.end() ? 0 : (*it).second;
  }

  EntityTypeId getModelGroup(const std::string &model_or_device_type_name);

  void addArtificialParameter(const std::string &name, ArtificialParameters::ArtificialParameter *artificial_parameter) {
    artificialParameterMap_[name] = artificial_parameter;
    passThroughParamsMap_[name] = 1;
  }

private:
  bool getParam(const std::string &name, double &value) const;
  Device &getDeviceByModelType(const EntityTypeId model_type);

  bool setupSolverInfo_ ();
  bool setupRawVectorPointers_ ();
  bool setupRawMatrixPointers_ ();

  bool updateIntermediateVars_();
  bool updatePrimaryState_();
  bool updateSecondaryState_();

  bool updateDependentParameters_();

  // Do the actual solve/calculation for the external devices
  void updateExternalDevices_();

  // add external devices for processors that don't own it
  ExternDevice::Instance * addExtDeviceInstance_(const InstanceBlock & IB);

private:
  IO::CmdParse &              commandLine_;           ///< Command line
  DeviceOptions               devOptions_;            ///< user-defined options:
  EntityTypeIdDeviceMap       deviceMap_;

  bool                        sensFlag_;              ///< .SENS present in netlist
  bool                        linearSystemFlag_;      ///< True if all devices in netlist have isLinearDevice() true
  bool                        firstDependent_;        ///< True until updateDependentParameters_ is called.
  bool                        parameterChanged_;      ///< Only used locally in updateDependentParameters_, don't know if stateful of just a member for fun
  bool                        breakPointInstancesInitialized;
  double                      timeParamsProcessed_;   ///< Time updateDependentParameters was called

  ExternData                  externData_;
  MatrixLoadData              matrixLoadData_;        ///< temporary jacobian load structures:
  SolverState                 solState_;              ///< real time solver data:
  Globals &                   globals_;               ///< global variables
  SolverState                 solStateExternal_;
  bool                        externalStateFlag_;

  N_LAS_Vector *              numJacSolVectorPtr_;
  N_LAS_Vector *              numJacStaVectorPtr_;
  N_LAS_Vector *              numJacStoVectorPtr_;
  N_LAS_Vector *              diagonalVectorPtr_;
  N_LAS_System *              lasSysPtr_;
  N_ANP_AnalysisManager *     anaIntPtr_;
  IO::OutputMgr *             outputMgrPtr_;
  IO::Measure::Manager *      measureManager_;
  N_PDS_Manager *             pdsMgrPtr_;
  Nonlinear::Manager *        nlsMgrPtr_;
  IO::PkgOptionsMgr *         pkgOptMgrPtr_;
  std::vector<std::pair<int, double> > *                icLoads_;

  ModelTypeNameModelTypeIdMap modelTypeMap_;          ///< Model type name to model
  ModelTypeNameModelTypeIdMap modelGroupMap_;         ///< Model type name to model group

  std::map<std::string, int>     localDeviceCountMap_;
  std::multimap<int, DeviceInstance *> solDevInstMap_;

  DeviceVector devicePtrVec_;
  DeviceVector pdeDevicePtrVec_;

  InstanceVector instancePtrVec_;
  InstanceVector bpInstancePtrVec_;                   ///< instances with breakpoints functions
  InstanceVector pdeInstancePtrVec_;
  InstanceVector nonPdeInstancePtrVec_;
  InstanceVector plotFileInstancePtrVec_;

  ModelTypeInstanceVectorMap  modelGroupInstanceVector_;
  ModelTypeInstanceVectorMap  modelTypeInstanceVector_;

  std::vector<SourceInstance *> indepSourceInstancePtrVec_;
  std::map<std::string, SourceInstance *, LessNoCase> indepSourcePtrMap_;
  // this is used to store the contents of the indepSourceInstancePtrVec_
  // during an mpde initialization where we'll remove slow sources from
  // the that vector so that they don't get updated
  std::vector<SourceInstance *> indepSourceInstanceBackupPtrVec_;

  std::set<std::string> devicesNeedingLeadCurrentLoads_;

  std::map<std::string, int, LessNoCase> passThroughParamsMap_;

  mutable DeviceEntityMap     parameterDeviceCache_;                  ///< Full parameter name to device entity cache

  // vector of pointers to devices under test.
  InstanceVector              testJacDevicePtrVec_;

  ModelVector                 modelVector_;
  ModelTypeModelVectorMap     modelGroupModelVector_;
  ModelTypeModelVectorMap     modelTypeModelVector_;

  EntityVector                dependentPtrVec_;

  std::vector<int>            numInterfaceNodes_;
  bool                        calledBeforeCSPI;

  // sensitivities:
  DeviceSensitivities *       devSensPtr_;

  // MPDE fast source list
  std::vector<std::string>    fastSourceNames_;

  // Device mask flag:
  bool                        nonTrivialDeviceMaskFlag;

  // .OP output flags.
  bool                        dotOpOutputFlag;

  // solution variable names vector:
  std::vector<std::string>      nameVec_;

  ArtificialParameterMap        artificialParameterMap_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/00
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerLinearSystem (N_LAS_System * tmp_system_ptr)
{
  lasSysPtr_ = tmp_system_ptr;

  if (lasSysPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerParallelMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/29/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerParallelMgr (N_PDS_Manager * tmp_pdsMgrPtr )
{
  pdsMgrPtr_              = tmp_pdsMgrPtr;

  return pdsMgrPtr_ != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerNonlinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/23/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerNonlinearSolver (Nonlinear::Manager * tmp_nlsMgrPtr)
{
  nlsMgrPtr_              = tmp_nlsMgrPtr;

  if (nlsMgrPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/04
//-----------------------------------------------------------------------------
//inline void DeviceMgr::setVoltageLimiterFlag ()
//{
//  devOptions_.voltageLimiterFlag = true;
//}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::unsetVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/04
//-----------------------------------------------------------------------------
inline void DeviceMgr::unsetVoltageLimiterFlag ()
{
  devOptions_.voltageLimiterFlag = false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/18/14
//-----------------------------------------------------------------------------
inline void DeviceMgr::setVoltageLimiterFlag ( bool flagVal )
{
  devOptions_.voltageLimiterFlag = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerICLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/03/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerICLoads( std::vector< std::pair<int,double> > * icLoads )
{
  return (icLoads_ = icLoads) != 0;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceMgr N_DEV_DeviceMgr;

#endif // Xyce_N_DEV_DeviceMgr_h
