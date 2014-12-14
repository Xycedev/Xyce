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
// Filename       : $RCSfile: N_DEV_DeviceMgr.C,v $
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
// Revision Number: $Revision: 1.633.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:49 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_InstanceName.h>
#include <N_DEV_Op.h>
#include <N_DEV_Print.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_Source.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_IO_CmdParse.h>
#include <N_IO_Op.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Op.h>

// These are slowly being removed as DeviceMgr should not depend on Devices.  Devices should plug in to DeviceMgr.
#include <N_DEV_Bsrc.h>
#include <N_DEV_ISRC.h>
#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_MOSFET_B3SOI.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_Vsrc.h>


namespace Xyce {
namespace Device {

namespace {

// structs related to options registration
struct OptionsReg : public IO::PkgOptionsReg
{
  OptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr->registerOptions(options);}

  DeviceMgr * devMgr;
};

struct SensOptionsReg : public IO::PkgOptionsReg
{
  SensOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr->registerSensParams(options);}

  DeviceMgr * devMgr;
};

struct TimeOptionsReg : public IO::PkgOptionsReg
{
  TimeOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr->registerTimeOptions(options);}

  DeviceMgr * devMgr;
};

// .TRAN
struct TransAnalysisReg : public IO::PkgOptionsReg
{
  TransAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return Mgr->setTranAnalysisParams(options);}

  DeviceMgr * Mgr;
};

// .DC
struct DCAnalysisReg : public IO::PkgOptionsReg
{
  DCAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return Mgr->setDCAnalysisParams(options);}

  DeviceMgr * Mgr;
};

// .OP
struct OPAnalysisReg : public IO::PkgOptionsReg
{
  OPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return Mgr->setOPAnalysisParams(options);}

  DeviceMgr * Mgr;
};

// .STEP
struct STEPAnalysisReg : public IO::PkgOptionsReg
{
  STEPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return Mgr->setSTEPAnalysisParams(options);}

  DeviceMgr * Mgr;
};

// .MPDE
struct MPDE_AnalysisReg : public IO::PkgOptionsReg
{
  MPDE_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr_->setMPDEAnalysisParams(options);}

  DeviceMgr * devMgr_;
};

// .HB
struct HB_AnalysisReg : public IO::PkgOptionsReg
{
  HB_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr_->setHBAnalysisParams(options);}

  DeviceMgr * devMgr_;
};

// .AC
struct AC_AnalysisReg : public IO::PkgOptionsReg
{
  AC_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr_->setACAnalysisParams(options);}

  DeviceMgr * devMgr_;
};

// .MOR
struct MOR_AnalysisReg : public IO::PkgOptionsReg
{
  MOR_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

  bool operator()(const Util::OptionBlock & options)
  {return devMgr_->setMORAnalysisParams(options);}

  DeviceMgr * devMgr_;
};

// ----------------------------------------------------------------------------
// Function      : setupIOName
//
// Purpose       : This function takes the device instance name and creates
//                 an appropriate "outputName" to be used for file outputs.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/26/2013
// ----------------------------------------------------------------------------
const std::string setupIOName(const InstanceName &entity_name)
{
  return entity_name.getDeviceName();
}

DeviceEntity *findDeviceEntity(EntityTypeIdDeviceMap::const_iterator begin, EntityTypeIdDeviceMap::const_iterator end, const std::string entity_name) {
  for (EntityTypeIdDeviceMap::const_iterator it = begin; it != end; ++it) {
    DeviceEntity *device_entity = (*it).second->findInstance(InstanceName(entity_name));
    if (device_entity)
      return device_entity;

    device_entity = (*it).second->findModel(ModelName(entity_name));
    if (device_entity)
      return device_entity;
  }

  return 0;
}

} // namespace <unnamed>

struct DeviceGlobalParameterOpBuilder : public Util::Op::Builder
{
  typedef std::list<Util::Param> ParameterList;

  DeviceGlobalParameterOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  virtual ~DeviceGlobalParameterOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceMgrGlobalParameterOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "GLOBAL_PARAMETER") {
      new_op  = new DeviceMgrGlobalParameterOp(param_string, deviceManager_, param_string);
      new_op->addArg(param_string);
    }
    else
    {
      const double *result = deviceManager_.findGlobalPar(param_tag);
      if (result)
      {
        // Refactor: [DGB] Should use the address.
        new_op = new DeviceMgrGlobalParameterOp(param_tag, deviceManager_, param_tag);
      }
    }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

struct DeviceEntityOpBuilder : public Util::Op::Builder
{
  typedef std::list<Util::Param> ParameterList;

  DeviceEntityOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  virtual ~DeviceEntityOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceEntityParameterOp>();
  }
  
  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    const DeviceEntity *device_entity = deviceManager_.getDeviceEntity(param_tag);
    if (device_entity)
    {
      std::string param_name = Util::paramNameFromFullParamName(param_tag);
      new_op = new DeviceEntityParameterOp(param_tag, *device_entity, param_name);
    }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

namespace {

typedef std::list<Util::Param> ParameterList;

void parameterNameAndArgs(std::string &name, std::vector<std::string> &args, ParameterList::const_iterator &it) 
{
  const std::string &param_tag = (*it).tag();

  if (param_tag[0] == 'V' || param_tag[0] == 'I' || param_tag[0] == 'N')
  {
    std::ostringstream oss;
    oss << param_tag << "(";
    int arg_count = (*it).getImmutableValue<int>();
    for (int i = 0; i < arg_count; ++i)
    {
      ++it;
      if (i != 0)
        oss << ",";
      oss << (*it).tag();
      args.push_back((*it).tag());
    }
    oss << ")";
    name = oss.str();
  }
}

} // namespace <unnamed>

struct DeviceMgrOpBuilder : public Util::Op::Builder
{
  DeviceMgrOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  virtual ~DeviceMgrOpBuilder()
  {}

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (deviceManager_.findParam(param_tag))
    {
      new_op = new DeviceMgrParameterOp(param_tag, deviceManager_, param_tag);
    }
    // else
    // {
    //   // Get arguments and build pretty name
    //   std::vector<std::string> args;
    //   std::string name;
    //   parameterNameAndArgs(name, args, it);

    //   // this is confusing.  While the solution vector and state/store vector's
    //   // use maps with modified device names(specifically where the device
    //   // type is always first as in "D:subcircuitname:devciename" as apposed to
    //   // "subcircuitname:Ddevicename", the device manager does not use the modified
    //   // device name to find a device.  So, when we set up the device name below
    //   // use the "nodeName" rather than the "modifiedDeviceName".
    //   std::string store_name = args[0] + ":DEV_" + param_tag;
    //   // have to repeat this check for spaces as in I(YPDE NAME)
    //   std::string::size_type space = store_name.find_first_of(" ");
    //   if (space != std::string::npos)
    //   {
    //     if (space == 4 && store_name.substr(0, 4) == "YPDE")
    //     {
    //       store_name.replace(4, 1, ":");
    //     }
    //   }

    //   if (deviceManager_.findParam(store_name))
    //   {
    //     new_op = new DeviceMgrParameterOp(param_tag, deviceManager_, store_name);
    //   }

    //   std::string ppde("YPDE!" + store_name);
    //   if (!new_op && deviceManager_.findParam(ppde))
    //   {
    //     {
    //       new_op = new DeviceMgrParameterOp(param_tag, deviceManager_, ppde);
    //     }
    //   }
    //   if (new_op)
    //     new_op->addArg(args[0]);
    // }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/15/06
//-----------------------------------------------------------------------------
bool DeviceMgr::registerOutputMgr (IO::OutputMgr * tmp_outputMgrPtr)
{
  outputMgrPtr_ = tmp_outputMgrPtr;

  outputMgrPtr_->getOpBuilderManager().addBuilder(new DeviceGlobalParameterOpBuilder(*this));
  outputMgrPtr_->getOpBuilderManager().addBuilder(new DeviceEntityOpBuilder(*this));
  
  return outputMgrPtr_ != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerMeasureMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/15/06
//-----------------------------------------------------------------------------
bool DeviceMgr::registerMeasureMgr (IO::Measure::Manager * measure_manager)
{
  measureManager_ = measure_manager;

  return measureManager_ != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerAnalysisManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::registerAnalysisManager (N_ANP_AnalysisManager * tmp_anaIntPtr)
{
  anaIntPtr_              = tmp_anaIntPtr;

  if (anaIntPtr_)
    static_cast<Xyce::Util::Notifier<Analysis::StepEvent> &>(*anaIntPtr_).subscribe(*this);
  
  return anaIntPtr_ != 0;
}
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::factory
// Purpose       : factory function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr * DeviceMgr::factory(IO::CmdParse & cp)
{
   return new DeviceMgr(cp);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::DeviceMgr
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::DeviceMgr(IO::CmdParse &command_line)
  : commandLine_(command_line),
    devOptions_(),
    icLoads_(NULL),
    devSensPtr_(new DeviceSensitivities(*this, devOptions_)),
    calledBeforeCSPI (false),
    sensFlag_(false),
    linearSystemFlag_(true),
    firstDependent_(true),
    externalStateFlag_(false),
    parameterChanged_(false),
    breakPointInstancesInitialized(false),
    timeParamsProcessed_(0.0),
    nonTrivialDeviceMaskFlag(false),
    dotOpOutputFlag(false),
    numJacStaVectorPtr_(0),
    numJacSolVectorPtr_(0),
    numJacStoVectorPtr_(0),
    diagonalVectorPtr_(0),
    globals_(solState_.globals_)
{
  addArtificialParameter("MOSFET:GAINSCALE", new ArtificialParameters::MOSFETGainScaleParam());
  addArtificialParameter("MOSFET:GAIN", new ArtificialParameters::MOSFETGainScaleParam());
  addArtificialParameter("MOSFET:NLTERMSCALE", new ArtificialParameters::MOSFETNLTermScaleParam());
  addArtificialParameter("MOSFET:NLTERM", new ArtificialParameters::MOSFETNLTermScaleParam());
  addArtificialParameter("MOSFET:L", new ArtificialParameters::MOSFETLParam());
  addArtificialParameter("MOSFET:W", new ArtificialParameters::MOSFETWParam());
  addArtificialParameter("MOSFET:SIZESCALE", new ArtificialParameters::MOSFETSizeScaleParam());
  addArtificialParameter("MOSFET:TOX", new ArtificialParameters::MOSFETTOXParam());
  addArtificialParameter("BJT:BF", new ArtificialParameters::BJTBFParam());
  addArtificialParameter("BJT:NF", new ArtificialParameters::BJTNFParam());
  addArtificialParameter("BJT:NR", new ArtificialParameters::BJTNRParam());
  addArtificialParameter("BJT:EXPORD", new ArtificialParameters::BJTExpOrdParam());
  addArtificialParameter("DIODE:N", new ArtificialParameters::DiodeNParam());
  addArtificialParameter("VSRCSCALE", new ArtificialParameters::VsrcScaleParam());
  addArtificialParameter("PDEALPHA", new ArtificialParameters::PDEAlphaParam());
  addArtificialParameter("PDEBETA", new ArtificialParameters::PDEBetaParam());
  addArtificialParameter("PDECHARGEALPHA", new ArtificialParameters::PDEChargeAlphaParam());
  addArtificialParameter("GSTEPPING", new ArtificialParameters::GSteppingParam());
  addArtificialParameter("GMIN", new ArtificialParameters::GMinParam());
  addArtificialParameter("TEMP", new ArtificialParameters::TempParam());

  passThroughParamsMap_[ "MOSFET_ALL:GAINSCALE"   ] = 1;
  passThroughParamsMap_[ "MOSFET_ALL:NLTERMSCALE" ] = 1;
  passThroughParamsMap_[ "MOSFET1:GAINSCALE"      ] = 1;
  passThroughParamsMap_[ "MOSFET1:NLTERMSCALE"    ] = 1;

  devOptions_.setupDefaultOptions(command_line);
  devOptions_.applyCmdLineOptions(command_line);

  solState_.debugTimeFlag = true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::~DeviceMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::~DeviceMgr()
{
  delete numJacStaVectorPtr_;
  delete numJacSolVectorPtr_;
  delete numJacStoVectorPtr_;
  delete diagonalVectorPtr_;

  delete externData_.numJacRHSVectorPtr;
  delete externData_.numJacFVectorPtr;
  delete externData_.numJacQVectorPtr;
  delete externData_.perturbVectorPtr;
  delete externData_.numJacLoadFlagPtr;

  delete externData_.tmpdIdXPtr;
  delete externData_.tmpdQdXPtr;

  for (EntityTypeIdDeviceMap::iterator it = deviceMap_.begin(), end = deviceMap_.end(); it != end; ++it)
    delete (*it).second;

  for (ArtificialParameterMap::iterator it = artificialParameterMap_.begin(), end = artificialParameterMap_.end(); it != end; ++it)
    delete (*it).second;

  delete icLoads_;
  delete devSensPtr_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 10/20/2008
//-----------------------------------------------------------------------------
bool DeviceMgr::registerPkgOptionsMgr(IO::PkgOptionsMgr *pkgOptPtr)
{
  pkgOptMgrPtr_ = pkgOptPtr;

  std::string netListFile("");

  if (commandLine_.getArgumentValue(std::string("netlist")) != "")
  {
    netListFile = commandLine_.getArgumentValue(std::string("netlist"));
  }

  pkgOptMgrPtr_->submitRegistration(
      "DEVICE", netListFile, new OptionsReg(this));

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new SensOptionsReg(this));

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new TimeOptionsReg(this));

  // different analysis types.
  pkgOptMgrPtr_->submitRegistration(
      "TRAN", netListFile, new TransAnalysisReg(this));

  pkgOptMgrPtr_->submitRegistration(
      "DC", netListFile, new DCAnalysisReg(this));

  pkgOptMgrPtr_->submitRegistration(
      "OP", netListFile, new OPAnalysisReg(this));

  pkgOptMgrPtr_->submitRegistration(
      "STEP", netListFile, new STEPAnalysisReg(this));

  // MPDE specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MPDE", netListFile, new MPDE_AnalysisReg(this));

  // HB Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "HB", netListFile, new HB_AnalysisReg(this));

  // AC Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "AC", netListFile, new AC_AnalysisReg(this));

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR", netListFile, new MOR_AnalysisReg(this));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerSensParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::registerSensParams (const Util::OptionBlock & OB)
{
  sensFlag_ = true;

  if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
  {
    dout() << "DeviceMgr::registerSensParams called!" <<std::endl;
  }

  return devSensPtr_->registerSensParams (OB);
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerLeadCurrentRequests
// Purpose       : this function is called from the output manager (through the
//                 device interface) to inform the device package of the devices
//                 for which lead currents have been requested.  The device
//                 manager will take care of doing isolated F and Q loads for
//                 these devices so the lead currents can be calculated
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 03/20/13
//-----------------------------------------------------------------------------
bool DeviceMgr::setLeadCurrentRequests(const std::set<std::string> & deviceNames)
{
  // this is called prior to fully constructing the devices.  So for now
  // save the list
  devicesNeedingLeadCurrentLoads_ = deviceNames;

  if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
  {
    if (! devicesNeedingLeadCurrentLoads_.empty())
    {
      dout() << "DeviceMgr::registerLeadCurrentRequests Devices for which lead currents were requested: ";
      for (std::set<std::string>::iterator it = devicesNeedingLeadCurrentLoads_.begin(), end = devicesNeedingLeadCurrentLoads_.end(); it != end; ++it)
      {
        dout() << (*it) << "  ";
      }
      dout() << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setTranAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setTranAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.TRANspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setDCAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setDCAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.DCspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setOPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setOPAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.OPspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setSTEPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setSTEPAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.STEPspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setMPDEAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.MPDEspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setHBAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setHBAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.HBspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setACAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setACAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.ACspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMORAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/30/12
//-----------------------------------------------------------------------------
bool DeviceMgr::setMORAnalysisParams (const Util::OptionBlock & OB)
{
  solState_.MORspecified = true;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getFastSourcePeriod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
std::vector<double> DeviceMgr::getFastSourcePeriod(std::vector<std::string> &sourceNames)
{
  int numFastSrcs = sourceNames.size();

  // Setup return of source periods
  std::vector<double> srcPeriods(numFastSrcs);

  // Now loop over them, and mark them.
  for(int i = 0; i < numFastSrcs; ++i)
  {
    std::map<std::string, SourceInstance *, LessNoCase>::iterator iterFS = indepSourcePtrMap_.find(sourceNames[i]);
    if (iterFS != indepSourcePtrMap_.end())
    {
      SourceInstance * SIPtr = iterFS->second;
      srcPeriods[i] = SIPtr->period();
    }
    else
    {
#ifndef Xyce_PARALLEL_MPI
      Report::UserError message;
      message << "Unable to find source: " <<  fastSourceNames_[i] << "\n"
              << "Potential names are: ";
      for (std::map<std::string,SourceInstance*, LessNoCase>::const_iterator it = indepSourcePtrMap_.begin(); it != indepSourcePtrMap_.end(); ++it)
        message << (*it).first << " ";
#endif
    }
  }

  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_MAX, srcPeriods);

// #ifdef Xyce_PARALLEL_MPI
//   // Collect all the periods from all the processors, assuming periods are positive values.
//   std::vector<double> tmpSrcPeriods( srcPeriods );
//   pdsMgrPtr_->getPDSComm()->maxAll( &tmpSrcPeriods[0], &srcPeriods[0], numFastSrcs );
// #endif

  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
std::vector<double> DeviceMgr::registerFastSources(std::vector<std::string> &sourceNames)
{
  int numFastSrcs = sourceNames.size();
  //Setup return of source periods
  std::vector<double> srcPeriods;

  if (numFastSrcs > 0)
  {
    srcPeriods.resize(numFastSrcs, 0.0);

      // Default case, the sources are explicitely listed on the .option line
    fastSourceNames_.resize(numFastSrcs);
    std::copy(sourceNames.begin(), sourceNames.end(), fastSourceNames_.begin());

    // Now loop over them, and mark them.
    for(int i = 0; i < numFastSrcs; ++i)
    {
      std::map<std::string,SourceInstance*, LessNoCase>::iterator iterFS = indepSourcePtrMap_.find(fastSourceNames_[i]);
      if (iterFS != indepSourcePtrMap_.end())
      {
        SourceInstance * SIPtr = iterFS->second;
        SIPtr->setFastSourceFlag (true);
        srcPeriods[i] = SIPtr->period();
      }
      else
      {
#ifndef Xyce_PARALLEL_MPI
      Report::UserError message;
      message << "Unable to find source: " << fastSourceNames_[i] << "\n"
              << "Potential names are: ";
      for (std::map<std::string,SourceInstance*, LessNoCase>::const_iterator it = indepSourcePtrMap_.begin(); it != indepSourcePtrMap_.end(); ++it)
        message << (*it).first << " ";
#endif
      }
    }

  }
  else
  {
    // tscoffe/tmei 09/16/08
    // Special case:  Use all sources
    // Compute the total number of fast sources for all processors.
    // NOTE:  In parallel, this will not work correctly if more than one processor has fast sources.
    int myNumFastSrcs = indepSourceInstancePtrVec_.size();
    pdsMgrPtr_->getPDSComm()->sumAll( &myNumFastSrcs, &numFastSrcs, 1 );
    srcPeriods.resize(numFastSrcs, -1.0);
    for (int i=0 ; i<myNumFastSrcs ; ++i) {
      indepSourceInstancePtrVec_[i]->setFastSourceFlag(true);
      srcPeriods[i] = indepSourceInstancePtrVec_[i]->period();
    }
  }

  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_MAX, srcPeriods);

  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deRegisterFastSources
// Purpose       : reverses the effect of registerFastSources
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void DeviceMgr::deRegisterFastSources (std::vector<std::string> &sourceNames)
{
  int numFastSrcs = sourceNames.size();

  if (numFastSrcs > 0)
  {
      // Default case, the sources are explicitely listed on the .option line
    fastSourceNames_.resize(numFastSrcs);
    std::copy(sourceNames.begin(), sourceNames.end(), fastSourceNames_.begin());

    // Now loop over them, and mark them.
    for(int i = 0; i < numFastSrcs; ++i)
    {
      std::map<std::string,SourceInstance*, LessNoCase>::iterator iterFS = indepSourcePtrMap_.find(fastSourceNames_[i]);
      if (iterFS != indepSourcePtrMap_.end())
      {
        SourceInstance * SIPtr = iterFS->second;
        SIPtr->setFastSourceFlag (false);
      }
      else
      {
#ifndef Xyce_PARALLEL_MPI
        Report::DevelFatal message;
        message.in("DeviceMgr::deRegisterFastSources");
        message << "Unable to find source: " <<  fastSourceNames_[i] << std::endl
                << "Potential names are: ";
        for (std::map<std::string,SourceInstance*, LessNoCase>::const_iterator it = indepSourcePtrMap_.begin(); it != indepSourcePtrMap_.end(); ++it)
          message << (*it).first << " ";
#endif
      }
    }
  }
  else
  {
    // Special case:  Use all sources
    numFastSrcs = indepSourceInstancePtrVec_.size();
    for (int i=0 ; i<numFastSrcs ; ++i) {
      indepSourceInstancePtrVec_[i]->setFastSourceFlag(false);
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deactivateSlowSources
// Purpose       : traverse fast source list and remove any slow sources from
//                 the deviceArray
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::deactivateSlowSources()
{
  // first back-up a copy of the deviceArray so we can edit out
  // the slow sources
  indepSourceInstanceBackupPtrVec_.resize(indepSourceInstancePtrVec_.size());
  std::copy(indepSourceInstancePtrVec_.begin(), indepSourceInstancePtrVec_.end(), indepSourceInstanceBackupPtrVec_.begin());

  // erase the existing list of sources
  indepSourceInstancePtrVec_.clear();

  // now copy back only those that are fast sources
  std::vector<SourceInstance*>::iterator iter;
  std::vector<SourceInstance*>::iterator begin =indepSourceInstanceBackupPtrVec_.begin();
  std::vector<SourceInstance*>::iterator end =indepSourceInstanceBackupPtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    if ((*iter)->getFastSourceFlag())
    {
      indepSourceInstancePtrVec_.push_back(*iter);
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::activateSlowSources
// Purpose       : restore any slow sources to the device array.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::activateSlowSources()
{
  // restore the independent source list from backup
  indepSourceInstancePtrVec_.clear();
  indepSourceInstancePtrVec_.resize(indepSourceInstanceBackupPtrVec_.size());
  std::copy(indepSourceInstanceBackupPtrVec_.begin(), indepSourceInstanceBackupPtrVec_.end(), indepSourceInstancePtrVec_.begin());
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setMPDEFlag(bool flagVal)
{
  solState_.mpdeOnFlag  = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setBlockAnalysisFlag(bool flagVal)
{
  solState_.blockAnalysisFlag = flagVal;
  devOptions_.setBlockAnalysisFlag(flagVal);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setFastTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setFastTime(double timeVal)
{
  solState_.currFastTime = timeVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::initializeAll
// Purpose       : This function, via the LAS system class, sets up
//                 the pointers to the various linear algebra entities.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool DeviceMgr::initializeAll()
{
  bool bsuccess = true;

  externData_.lasSysPtr = lasSysPtr_;

  // nullify ptrs that are passed in at each step  (see the loadDAEVectors function args)
  externData_.nextSolVectorPtr = 0;
  externData_.currSolVectorPtr = 0;
  externData_.lastSolVectorPtr = 0;
  externData_.daeQVectorPtr    = 0;
  externData_.daeFVectorPtr    = 0;
  externData_.daeBVectorPtr    = 0;
  externData_.dFdxdVpVectorPtr = 0;
  externData_.dQdxdVpVectorPtr = 0;
  externData_.nextStaVectorPtr = 0;
  externData_.currStaVectorPtr = 0;
  externData_.lastStaVectorPtr = 0;
  externData_.storeLeadCurrQCompPtr = 0;
  externData_.nextStaDerivVectorPtr = 0;
  externData_.nextStoVectorPtr =
  externData_.currStoVectorPtr = 0;
  externData_.lastStoVectorPtr = 0;

  if (DEBUG_DEVICE) {
    // get f vector pointer:
    externData_.fVectorPtr  = lasSysPtr_->getFVector();
    bsuccess = bsuccess && (externData_.fVectorPtr != 0);

    // create Jdxp vector pointer:
    externData_.JdxpVectorPtr = lasSysPtr_->getJDXPVector();
    bsuccess = bsuccess && (externData_.JdxpVectorPtr != 0);
  }

#ifdef Xyce_DEBUG_VOLTLIM
  // get test matrix: (old DAE)
  externData_.JTestMatrixPtr = lasSysPtr_->getJacTestMatrix();
  bsuccess = bsuccess && (externData_.JTestMatrixPtr != 0);

  // get test matrix:
  externData_.FTestMatrixPtr = lasSysPtr_->getdFdxTestMatrix();
  bsuccess = bsuccess && (externData_.FTestMatrixPtr != 0);
  externData_.QTestMatrixPtr = lasSysPtr_->getdQdxTestMatrix();
  bsuccess = bsuccess && (externData_.QTestMatrixPtr != 0);

  // get dxVoltlim vector pointer:
  externData_.dxVoltlimVectorPtr = lasSysPtr_->getDxVoltlimVector() ;
  bsuccess = bsuccess && (externData_.dxVoltlimVectorPtr != 0);

  // create Jdx2 vector pointer: (old DAE)
  externData_.Jdx2VectorPtr = lasSysPtr_->getJDX2Vector();
  bsuccess = bsuccess && (externData_.Jdx2VectorPtr != 0);

  // create Jdx2 vector pointer:
  externData_.Fdx2VectorPtr = lasSysPtr_->getFDX2Vector();
  bsuccess = bsuccess && (externData_.Fdx2VectorPtr != 0);
  externData_.Qdx2VectorPtr = lasSysPtr_->getQDX2Vector();
  bsuccess = bsuccess && (externData_.Qdx2VectorPtr != 0);
#endif

  // get flag solution pointer-pointer:
  externData_.flagSolVectorPtr = lasSysPtr_->getFlagSolVector();
  bsuccess = bsuccess && (externData_.flagSolVectorPtr != 0);

  // get device mask pointer.
  externData_.deviceMaskVectorPtr = lasSysPtr_->getDeviceMaskVector ();
  bsuccess = bsuccess && (externData_.deviceMaskVectorPtr != 0);

  // create the temporary numerical jacobian vectors
  if (devOptions_.numericalJacobianFlag || devOptions_.testJacobianFlag || sensFlag_)
  {
    numJacStaVectorPtr_ = lasSysPtr_->builder().createStateVector();
    numJacSolVectorPtr_ = lasSysPtr_->builder().createVector();
    numJacStoVectorPtr_ = lasSysPtr_->builder().createStoreVector();

    externData_.numJacRHSVectorPtr = lasSysPtr_->builder().createVector();
    externData_.numJacFVectorPtr   = lasSysPtr_->builder().createVector();
    externData_.numJacQVectorPtr   = lasSysPtr_->builder().createVector();
    externData_.perturbVectorPtr   = lasSysPtr_->builder().createVector();
    externData_.numJacLoadFlagPtr  = lasSysPtr_->builder().createVector();
  }

  externData_.tmpdIdXPtr = lasSysPtr_->builder().createVector();
  externData_.tmpdQdXPtr = lasSysPtr_->builder().createVector();

  // create a diagonal vector to be used for 2-level
  diagonalVectorPtr_  = lasSysPtr_->builder().createVector();

  externData_.initializeAllFlag = true;

  // For Homotopy on block gainscale
  solState_.InitializeHomotopyBlockSize(devOptions_.numGainScaleBlocks);

#ifdef Xyce_SIZEOF
  int size = sizeof(*this);
  dout() << "Size of device package after initializeAll  = " << size << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::resetForStepAnalysis
// Purpose       :
// Special Notes : Some "resetForStep" functions (only HB so far) will
//                 call dev->initializeAll.  So, this must be called first.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------

void DeviceMgr::notify(const Analysis::StepEvent &event) 
{
  if (event.state_ == Analysis::StepEvent::STEP_STARTED)
  {
    delete numJacStaVectorPtr_;  numJacStaVectorPtr_ = 0;
    delete numJacSolVectorPtr_;  numJacSolVectorPtr_ = 0;
    delete numJacStoVectorPtr_;  numJacStoVectorPtr_ = 0;
    delete externData_.numJacRHSVectorPtr;  externData_.numJacRHSVectorPtr = 0;
    delete externData_.numJacFVectorPtr;  externData_.numJacFVectorPtr = 0;
    delete externData_.numJacQVectorPtr;  externData_.numJacQVectorPtr = 0;
    delete externData_.perturbVectorPtr;  externData_.perturbVectorPtr = 0;
    delete externData_.numJacLoadFlagPtr;  externData_.numJacLoadFlagPtr = 0;
    delete externData_.tmpdIdXPtr;  externData_.tmpdIdXPtr = 0;
    delete externData_.tmpdQdXPtr;  externData_.tmpdQdXPtr = 0;
    delete diagonalVectorPtr_;  diagonalVectorPtr_ = 0;

    solState_.ltraTimeIndex = 0;
    solState_.ltraTimeHistorySize = 10;
    solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::createDevice
// Purpose       : This function creates a single device based on the passed
//                 index.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Device &DeviceMgr::getDeviceByModelType(const ModelTypeId model_type_id)
{
  Device *device = 0;

  if (model_type_id.defined())
  {
    EntityTypeIdDeviceMap::const_iterator it = deviceMap_.find(model_type_id);
    if (it == deviceMap_.end())
    {
      FactoryBlock factory_block(devOptions_, solState_, matrixLoadData_, externData_, commandLine_);

      const Configuration *configuration = Configuration::findConfiguration(model_type_id);
      device = configuration->createDevice(factory_block);
      deviceMap_[model_type_id] = device;
      devicePtrVec_.push_back(device);

      if (device->isPDEDevice())
      {
        pdeDevicePtrVec_.push_back(device);
      }
    }
    else
      device = (*it).second;
  }

  return *device;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getModelGroup
// Purpose       : This function returns the device type index for a given
//                 named string.  This assumes that the device names used
//                 in the simulation will obey the spice3f5 netlist language
//                 convention.
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/**
 * Return the ModelGroup of the device associated with the model type name or device type name.
 *
 * @param model_type_name model type name or device type name
 *
 * @return
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355
 * @date   Mon Sep 23 07:53:04 2013
 */
EntityTypeId DeviceMgr::getModelGroup(const std::string &model_or_device_type_name)
{
  return Configuration::getModelGroup(model_or_device_type_name);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::addDeviceModel(const ModelBlock & model_block)
{
  ModelTypeId model_type;
  ModelTypeId model_group = Configuration::getModelGroup(model_block.getType());

  if (!model_block.getName().empty()) {
    model_type = Configuration::getModelType(model_block.getType(), model_block.getLevel());

    if (!model_type.defined())
      Report::UserError() << "There is no device " << model_block.getType() << " of level " << model_block.getLevel() << " to define model " << model_block.getName();
  }

  if (!model_type.defined())
    model_type = model_group;

  if (!model_type.defined())
    return false;

  FactoryBlock factory_block(devOptions_, solState_, matrixLoadData_, externData_, commandLine_);
  Device &device = getDeviceByModelType(model_type);
  DeviceModel *device_model = device.addModel(model_block, factory_block);

  modelTypeMap_[model_block.getName()] = model_type;
  modelGroupMap_[model_block.getName()] = model_group;

  // add the various model vectors:
  if (device_model != 0)
  {
    modelVector_.push_back(device_model);
    modelGroupModelVector_[model_group].push_back(device_model);
    modelTypeModelVector_[model_type].push_back(device_model);
  }

  return device_model != 0;
}


std::pair<ModelTypeId, ModelTypeId>
DeviceMgr::getModelType(
  const InstanceBlock &         instance_block)
{
  // If the ModelName string is not a null string, use it to get the
  // device type index.  Otherwise, use the name of the device itself to determine the device.
  ModelTypeId model_type;
  ModelTypeId model_group;
  if (instance_block.getModelName().empty())
  {
    model_type = getModelGroup(instance_block.getInstanceName().getDeviceType());
    model_group = model_type;
  }
  else
  {
    model_type = modelTypeMap_[instance_block.getModelName()];
    model_group = modelGroupMap_[instance_block.getModelName()];
  }

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }
    return std::pair<ModelTypeId, ModelTypeId>(ModelTypeId(), ModelTypeId());
  }
  else if (instance_block.bsourceFlag)
  {
    // This is an E, F, G, or H source that is to be treated as a B source,
    // set its type now to BSRC.
    model_type = Bsrc::Traits::modelType();
  }

  // Check if this is a simple resistor, but with a resistance of zero.  If so,
  // then change the type to RESISTOR3.  This variant of the resistor acts like
  // a voltage souce with zero voltage difference.
  if ((devOptions_.checkForZeroResistance) && (model_type == Resistor::Traits::modelType()))
  {
    const double zeroResistanceValue = devOptions_.zeroResistanceTol;
    // loop over the parameters
    std::vector<Param>::const_iterator currentParam = instance_block.params.begin();
    std::vector<Param>::const_iterator endParam = instance_block.params.end();
    while(currentParam != endParam)
    {
      if ((currentParam->uTag() == "R"))
      {
        if (currentParam->given())
        {
          std::vector<std::string> variables, specials;

          // check if this is a time-dependent, or variable-dependent expression.
          // If it is, then skip.
          const Param * devPar = &(*(currentParam));
          const Util::Param * tmpPar = (dynamic_cast<const Util::Param*> (devPar));
          // only check if this is an expression-type parameter.
          if (tmpPar->getType() == Util::EXPR)
          {
            Util::Expression tmpExp = tmpPar->getValue<Util::Expression>();
            tmpExp.get_names(XEXP_VARIABLE, variables);
            tmpExp.get_names(XEXP_SPECIAL, specials);
          }

          if (specials.empty() && variables.empty())
          {
            if (fabs(currentParam->getImmutableValue<double>()) < devOptions_.zeroResistanceTol)
            {
              // change device type to be the level 3 resistor
              model_type = Resistor3::Traits::modelType();
            }
          }
        }
        break;
      }
      currentParam++;
    }
  }

  return std::pair<ModelTypeId, ModelTypeId>(model_type, model_group);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::verifyDeviceInstance
// Purpose       : This function verifies a device instance prior to
//                 instantiating.
//
//                 Theoretically, we could do this in addDeviceInstance() and
//                 not make a device if it fails some verification criteria (a
//                 resistor with zero resistance is the primary case here).
//                 However, later unlinking of redundant devices that are
//                 connected to just one node is difficult after
//                 addDeviceInstance() is called because the device instance
//                 pointer can be placed in many other containers
//                 It could be done, but this is a simpler first step to having the
//                 device manager be in charge of device verification -- rather
//                 than have it in toplogy or IO.
//
// Special Notes : return true if this device is ok to instantiate, false otherwise
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
bool DeviceMgr::verifyDeviceInstance(InstanceBlock & instance_block)
{
  ModelTypeId model_type;
  ModelTypeId model_group;
  std::pair<ModelTypeId, ModelTypeId> x = getModelType(instance_block);
  model_type = x.first;
  model_group = x.second;

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine model type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name" << instance_block.getModelName();
    }

    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceInstance
// Purpose       : addDeviceInstance will create a new instance of the
//                 designated device type.  This version of the function
//                 accepts a parameter list as one of the arguments,
//                 so it is assumed that a parameter instance will
//                 also have to be allocated for it.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceInstance * DeviceMgr::addDeviceInstance(InstanceBlock & instance_block)
{
  // If the ModelName string is not a null string, use it to get the
  // device type index.  Otherwise, use the name of the device itself to determine the device.
  ModelTypeId model_type;
  ModelTypeId model_group;
  std::pair<ModelTypeId, ModelTypeId> x = getModelType(instance_block);
  model_type = x.first;
  model_group = x.second;

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }

    DeviceInstance *instance = 0;
    return instance;
  }

  // Add an instance of this type.
  Device &device = getDeviceByModelType(model_type);
  DeviceInstance *instance = device.addInstance(instance_block, FactoryBlock(devOptions_, solState_, matrixLoadData_, externData_, commandLine_));

  std::string outputName = setupIOName(instance->getName());

  if ((devicesNeedingLeadCurrentLoads_.find(outputName) != devicesNeedingLeadCurrentLoads_.end()) ||
      (devOptions_.calculateAllLeadCurrents))
  {
    instance->enableLeadCurrentCalc();

    if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
    {
      dout() << "DeviceMgr::addDeviceInstance Enabling lead current load for device \""
        << instance->getName()
        << "\" ->  \""
        << outputName
        << "\"" << std::endl;
    }
  }
  else if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
  {
    dout() << "DeviceMgr::addDeviceInstance Cannot enable lead current load for device \""
           << instance->getName()
           << "\" ->  \""
           << outputName
           << "\""
           << std::endl;
  }

  localDeviceCountMap_[device.getDefaultModelName()]++;

  linearSystemFlag_ = linearSystemFlag_ && device.isLinearDevice();

  solState_.PDESystemFlag = solState_.PDESystemFlag || device.isPDEDevice();

  // Set up the instance vectors.  These are the main containers used in the load procedures.
  instancePtrVec_.push_back(instance);

  // set up the list of pde device instances
  // and the list of non-pde devices instances.
  if (device.isPDEDevice())
  {
    pdeInstancePtrVec_.push_back(instance);
  }
  else
  {
    nonPdeInstancePtrVec_.push_back(instance);
  }

  modelGroupInstanceVector_[model_group].push_back(instance);
  modelTypeInstanceVector_[model_type].push_back(instance);

  // set up the independent source map.
  if (model_type == Vsrc::Traits::modelType() || model_type == ISRC::Traits::modelType())
  {
    indepSourcePtrMap_[instance_block.getInstanceName().getEncodedName()] = dynamic_cast<SourceInstance*>(instance);
    indepSourceInstancePtrVec_.push_back(dynamic_cast<SourceInstance*>(instance));
  }

  if (instance->plotfileFlag ())
  {
    plotFileInstancePtrVec_.push_back(instance);
  }

  // Set up the vector of devices subject to the jacobian test.
  if (instance->getName() == devOptions_.testJacDeviceName)
  {
    testJacDevicePtrVec_.push_back(instance);
  }

  return instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deleteDeviceInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::deleteDeviceInstance(const std::string & name)
{
  Report::DevelFatal().in("DeviceMgr::deleteDeviceInstance") << "Not ready with the new boilerplate-free device package";

  bool bsuccess = true;

  // bool tmpBool = true;
  // EntityTypeId type = getModelGroup(name);

  // DeviceMap::iterator it = deviceMap_.find(type);
  // if (it != deviceMap_.end()) {
  //   // bsuccess &= (*it).second->deleteInstance(name);
  //   DeviceEntity *entity = (*it).second->findEntity(name);
  //   DeviceInstance *instance = dynamic_cast<DeviceInstance *>(entity);
  //   if (instance)
  //     bsuccess = (*it).second->deleteInstance(instance);
  // }

  // note as this is written it ignores lots of other containers
  // this may be used for clean up at the end of a run, but
  // it is not sufficient to use during simulation setup.
  //
  // need to remove pointer to the instance "name" from other arrays
  // candidate lists are:
  //     InstanceVector instancePtrVec_;
  //     InstanceVector bpInstancePtrVec_; // instances with breakpoints functions
  //     InstanceVector pdeInstancePtrVec_;
  //     InstanceVector nonPdeInstancePtrVec_;
  //     InstanceVector mosfetInstancePtrVec_;
  //     InstanceVector vsrcInstancePtrVec_;
  //     InstanceVector bjtInstancePtrVec_;
  //     std::map<std::string, VsrcInstance*> vsrcInstancePtrMap_;
  //
  //     InstanceVector plotFileInstancePtrVec_;
  //
  //     std::map<std::string,SourceInstance*> indepSourcePtrMap_;
  //     std::vector<SourceInstance*> indepSourceInstancePtrVec_;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput1()
{
   // dump the fvector and the jdxp vector to files.
  if (devOptions_.debugLevel > 3 && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int outIter = solState_.newtonIter + 1;

    // f-vector
    char fn_fv[256]; for (int ich = 0; ich < 256; ++ich) fn_fv[ich] = 0;
    sprintf(fn_fv, "fvector.%03d.txt", outIter);

    // Note: this needs to change sign to match Spice.
    (externData_.fVectorPtr)->scale(-1.0);
    (externData_.fVectorPtr)->writeToFile(fn_fv);
    (externData_.fVectorPtr)->scale(-1.0);

    // jdxp-vector
    char fn_jdxp[256]; for (int ich = 0; ich < 256; ++ich) fn_jdxp[ich] = 0;
    sprintf(fn_jdxp, "Jdxp.%03d.txt", outIter);
    (externData_.JdxpVectorPtr)->writeToFile(fn_jdxp);

#ifdef Xyce_DEBUG_VOLTLIM
    // the voltlim dx vector.
    char fn_dxvl[256]; for (int ich = 0; ich < 256; ++ich) fn_dxvl[ich] = 0;
    sprintf(fn_dxvl, "dxVL.%03d.txt", outIter);
    (externData_.dxVoltlimVectorPtr)->writeToFile(fn_dxvl);

    // jdx2-vector
    char fn_jdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_jdx2[ich] = 0;
    sprintf(fn_jdx2, "Jdx2.%03d.txt", outIter);
    (externData_.Jdx2VectorPtr)->writeToFile(fn_jdx2);
#endif

  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput2
// Purpose       : new-dae version
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/19/08
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput2()
{
  // dump the fvector and the jdxp vector to files.
  if (devOptions_.debugLevel > -1 && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int newtonIter = solState_.newtonIter + 1;
    int outputStepNumber = 0;

    if (solState_.tranopFlag)
    {
      outputStepNumber = 0;
    }
    else if (solState_.initTranFlag)
    {
      outputStepNumber = solState_.timeStepNumber+1;
    }
    else
    {
      outputStepNumber = solState_.timeStepNumber+1;
    }

#ifdef Xyce_DEBUG_VOLTLIM
    // the voltlim dx vector.
    char fn_dxvl[256]; for (int ich = 0; ich < 256; ++ich) fn_dxvl[ich] = 0;
    sprintf(fn_dxvl, "dxVL.%03d.%03d.txt", outputStepNumber, newtonIter);
    (externData_.dxVoltlimVectorPtr)->writeToFile(fn_dxvl);

    // Fdx2-vector
    char fn_fdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_fdx2[ich] = 0;
    sprintf(fn_fdx2, "Fdx2.%03d.%03d.txt", outputStepNumber, newtonIter);
    (externData_.Fdx2VectorPtr)->writeToFile(fn_fdx2);

    // Qdx2-vector
    char fn_qdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_qdx2[ich] = 0;
    sprintf(fn_qdx2, "Qdx2.%03d.%03d.txt", outputStepNumber, newtonIter);
    (externData_.Qdx2VectorPtr)->writeToFile(fn_qdx2);
#endif

  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setInitialGuess
// Purpose       : This is a function call that sets the initial guess for
// devices that have initial guesses.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setInitialGuess (N_LAS_Vector * solVectorPtr)
{
  bool bsuccess = true;

  if (solVectorPtr != 0)
  {
    externData_.nextSolVectorPtr = solVectorPtr;

    // if two-level, and just the inner problem, only load the PDE devices.
    InstanceVector::iterator iter;
    InstanceVector::iterator begin;
    InstanceVector::iterator end;
    begin = instancePtrVec_.begin ();
    end = instancePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      bool tmpBool = (*iter)->setInitialGuess ();
      bsuccess = bsuccess && tmpBool;
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::analyticSensitivitiesAvailable
// Purpose       : 
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/30/2014
//-----------------------------------------------------------------------------
bool DeviceMgr::analyticSensitivitiesAvailable (std::string & name)
{
  // This assumes that the only available analytic sensitivities are for 
  // non-artificial parameters.
  DeviceEntity * device_entity = getDeviceEntity(name);

  int entityFound = (device_entity != 0);
  int available=0;
  if (entityFound)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);
    available = device_entity->analyticSensitivityAvailable (paramName);
  }
  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_LOR, &available, 1);

  return static_cast<bool>(available);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getAnalyticSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/30/2014
//-----------------------------------------------------------------------------
void DeviceMgr::getAnalyticSensitivities
  (std::string & name, 
   std::vector<double> & dfdpVec, 
   std::vector<double> & dqdpVec,
   std::vector<double> & dbdpVec,
   std::vector<int> & FindicesVec,
   std::vector<int> & QindicesVec,
   std::vector<int> & BindicesVec)
{
  // If not artificial, then search for the appropriate natural param(s).
  DeviceEntity * device_entity = getDeviceEntity(name);

  int entityFound = (device_entity != 0);

  if (entityFound)
  {
    bool found;
    std::string paramName = Util::paramNameFromFullParamName(name);
    found = device_entity->getAnalyticSensitivity (paramName, 
        dfdpVec, dqdpVec, dbdpVec, 
        FindicesVec, QindicesVec, BindicesVec);
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setParam
//
// Purpose       : This function sets named parameters (name) to a
//                 specified value (val).
//
// Special Notes : Used for continuation calculations, as well as possibly
//                 intrusive sensitivity/optimization calculations.  It is
//                 assumed that this is called after everything (devices,
//                 solvers, etc.) is set up.
//
//                 The specified parameter can be either a natural or
//                 artificial parameter.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setParam(std::string & name, double val)
{
  bool bsuccess = true, success = true;

  if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
  {
    std::string netListFile("");
    if (commandLine_.getArgumentValue(std::string("netlist")) != "")
    {
      netListFile = commandLine_.getArgumentValue(std::string("netlist"));
    }

    dout() << netListFile << "\t\t";
    dout() << "DeviceMgr::setParam.  ";
    dout() << name;
    dout() << "  " << val;
    dout() << std::endl;
  }

  ArtificialParameterMap::iterator artificial_param_it = artificialParameterMap_.find(name);
  if (artificial_param_it != artificialParameterMap_.end())
  {
    (*artificial_param_it).second->setValue(*this, val);
  }
  else
  {
    GlobalParameterMap::iterator global_param_it = globals_.global_params.find(name);
    if (global_param_it != globals_.global_params.end())
    {
      if ((*global_param_it).second != val)
      {
        (*global_param_it).second = val;
        for (EntityVector::iterator it = dependentPtrVec_.begin(), end = dependentPtrVec_.end(); it != end; ++it)
        {
          if ((*it)->updateGlobalParameters(globals_.global_params));
          {
            (*it)->processParams();
            (*it)->processInstanceParams();
          }
        }
      }
    }
    else
    {
      // If not artificial, then search for the appropriate natural param(s).
      DeviceEntity * device_entity = getDeviceEntity(name);

      int entityFound = (device_entity != 0);
      if (entityFound)
      {
        bool found;
        std::string paramName = Util::paramNameFromFullParamName(name);
        if (paramName == "")
        {
          found = device_entity->setDefaultParam(val);
        }
        else
        {
          found = device_entity->setParam(paramName, val);
        }

        if (found)
        {
          device_entity->processParams (); // if this "entity" is a model, then need to
          // also do a  "processParams" on the related instances.
          device_entity->processInstanceParams();
        }

        entityFound = found;
      }

      Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_LOR, &entityFound, 1);

      if (entityFound == 0)
      {
        if (DEBUG_DEVICE)
        {
          Report::DevelWarning() << "DeviceMgr::setParam.  Unable to find parameter " << name;
        }
        else
        {
          Report::UserError() << "Unable to find parameter " << name;
        }
      }
    }
  }

  // Certain parameters should be passed through to the inner solve,
  // if there is an inner circuit problem.  The names of parameters
  // that should be passed through are stored in the map.
  if (passThroughParamsMap_.find(name) != passThroughParamsMap_.end())
  {
    ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
    if (model_type_it != modelTypeInstanceVector_.end()) {
      for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

        bool bs1 = extern_device.setInternalParam(name, val);
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParam
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
//                 If the parameter is not found, this function sets val to
//                 0.0, and returns a "false" boolean.  It does not invoke the
//                 error handler.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/13/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getParam(const std::string &   name, double & value) const
{
  value = 0.0;
  bool entityFound = false;

  ArtificialParameterMap::const_iterator artificial_param_it = artificialParameterMap_.find(name);
  if (artificial_param_it != artificialParameterMap_.end())
  {
    value = (*artificial_param_it).second->getValue(*this);
    entityFound = true;
  }
  else
  {
    GlobalParameterMap::const_iterator global_param_it = globals_.global_params.find(name);
    if (global_param_it != globals_.global_params.end())
    {
      value = (*global_param_it).second;
      entityFound = true;
    }
    else
    {
      const DeviceEntity *device_entity = getDeviceEntity(name);
      entityFound = device_entity != 0;

      if (entityFound)
      {
        std::string paramName = Util::paramNameFromFullParamName(name);
        entityFound = device_entity->getParam(paramName, value);
        if (paramName == "")
          entityFound = true;
      }

      if (!entityFound)
        measureManager_->getMeasureValue(name, value, entityFound);
    }
  }

  return entityFound;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::findParam
// Purpose       : Returns the current value of a named parameter.
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::findParam(const std::string & name) const 
{
  double value = 0.0;
  return getParam(name, value);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParamNoReduce
// Purpose       : Returns the current value of a named parameter.
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
double DeviceMgr::getParamNoReduce(const std::string & name) const 
{
  double value = 0.0;
  return getParam(name, value);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParam
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
bool DeviceMgr::getParamAndReduce(const std::string & name, double & value) const
{
  bool found = getParam(name, value);
  int found_count = found ? 1 : 0;
  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_SUM, &found_count, 1);

  if (found_count) {
    Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_SUM, &value, 1);
    value /= found_count;
  }
  
  return found_count != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParamAndReduce
// Purpose       : Returns the current value of a named parameter.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
double DeviceMgr::getParamAndReduce(const std::string & name) const
{
  double value = 0.0;
  bool found = getParamAndReduce(name, value);

  if (!found)
  {
    if (DEBUG_DEVICE)
    {
      Report::DevelWarning() << 
        "DeviceMgr::getParamAndReduce.  Unable to find parameter " << name;
    }
    else
    {
      Report::UserError() << "Unable to find parameter " << name;
    }
  }

  return value;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateState
// Purpose       : This should be called prior to loadDAEVectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceMgr::updateState (
     N_LAS_Vector * nextSolVectorPtr,
     N_LAS_Vector * currSolVectorPtr,
     N_LAS_Vector * lastSolVectorPtr,
     N_LAS_Vector * nextStaVectorPtr,
     N_LAS_Vector * currStaVectorPtr,
     N_LAS_Vector * lastStaVectorPtr,
     N_LAS_Vector * nextStoVectorPtr,
     N_LAS_Vector * currStoVectorPtr,
     N_LAS_Vector * lastStoVectorPtr
    )
{
  bool bsuccess = true;
  bool tmpBool = true;

  tmpBool = setupSolverInfo_();
  bsuccess = bsuccess && tmpBool;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = nextSolVectorPtr;
  externData_.currSolVectorPtr = currSolVectorPtr;
  externData_.lastSolVectorPtr = lastSolVectorPtr;
  externData_.nextStaVectorPtr = nextStaVectorPtr;
  externData_.currStaVectorPtr = currStaVectorPtr;
  externData_.lastStaVectorPtr = lastStaVectorPtr;
  externData_.nextStoVectorPtr = nextStoVectorPtr;
  externData_.currStoVectorPtr = currStoVectorPtr;
  externData_.lastStoVectorPtr = lastStoVectorPtr;

#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
#endif

  // Now reset the relevant RAW pointers:
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  externData_.currSolVectorRawPtr = &((*externData_.currSolVectorPtr)[0]);
  externData_.lastSolVectorRawPtr = &((*externData_.lastSolVectorPtr)[0]);
  externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  externData_.currStaVectorRawPtr = &((*externData_.currStaVectorPtr)[0]);
  externData_.lastStaVectorRawPtr = &((*externData_.lastStaVectorPtr)[0]);
  externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);
  externData_.currStoVectorRawPtr = &((*externData_.currStoVectorPtr)[0]);
  externData_.lastStoVectorRawPtr = &((*externData_.lastStoVectorPtr)[0]);


  updateDependentParameters_();

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM)
  {
    // call all the intermediate vars loads:
    std::vector<Device*>::iterator iter;
    std::vector<Device*>::iterator begin = pdeDevicePtrVec_.begin ();
    std::vector<Device*>::iterator end = pdeDevicePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      tmpBool = (*iter)->updateState (externData_.nextSolVectorRawPtr,
                                      externData_.nextStaVectorRawPtr, 
                                      externData_.nextStoVectorRawPtr);
      bsuccess = bsuccess && tmpBool;
    }
  }
  else
  {
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->updateState (externData_.nextSolVectorRawPtr,
                                                 externData_.nextStaVectorRawPtr, 
                                                 externData_.nextStoVectorRawPtr);
    }
  }

  updateExternalDevices_();

#ifdef Xyce_PARALLEL_MPI
  externData_.nextStaVectorPtr->importOverlap();
  externData_.nextStoVectorPtr->importOverlap();
#endif

  Report::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEMatrices
// Purpose       : This function loads the various DAE related matrices to
//                 set up the following expression:
//
//                 residual:  f(x) = dQ/dt + F(x) - B(t) = 0
//
//                 jacobian:  J(x) = d(dQ/dt)dx + dFdx
//                                 = d(dQdx)dt + dFdx
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadJacobianMatrix
//
//                 Note that this function, unlike the loadJacobianMatrix
//                 function *requires* that vector/matrix pointers
//                 be passed in to the function.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load matrices
//                 are zeroed out.  That is, dFdx, dQdx are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEMatrices
  (N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpStateVectorPtr,
   N_LAS_Vector * tmpStateDerivVectorPtr,
   N_LAS_Vector * tmpStoreVectorPtr,
   N_LAS_Matrix * tmpdQdxMatrixPtr,
   N_LAS_Matrix * tmpdFdxMatrixPtr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // copy over the passed pointers:
  (externData_.nextSolVectorPtr) = tmpSolVectorPtr;
  bool resetRawMatrixPointers = true;
  if (
      (externData_.dQdxMatrixPtr == tmpdQdxMatrixPtr) &&
      (externData_.dFdxMatrixPtr == tmpdFdxMatrixPtr)
    )
  {
    resetRawMatrixPointers = false;
  }
  if (resetRawMatrixPointers)
  {
    externData_.dQdxMatrixPtr = tmpdQdxMatrixPtr;
    externData_.dFdxMatrixPtr = tmpdFdxMatrixPtr;
  }

  externData_.nextStaVectorPtr = tmpStateVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStateDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoreVectorPtr;

  // setup the relevant RAW vector pointers:
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  externData_.nextStaDerivVectorRawPtr = &((*externData_.nextStaDerivVectorPtr)[0]);
  externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);

//#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // setup the relevant RAW matrix pointers (down in the devices that need them):
  if (resetRawMatrixPointers || solState_.blockAnalysisFlag)
  {
    this->setupRawMatrixPointers_();
  }
//#endif

  InstanceVector::iterator iter;
  InstanceVector::iterator begin;
  InstanceVector::iterator end;

  // if this is an "inner problem" phase of a Two-Level Newton
  // simulation, then only load the PDE devices.  Everything else just
  // gets "1" on the diagonal.
  //
  // Note, it is possible to just load 1's using a petra call, so I may
  // get rid of the "trivial" matrix stamp stuff soon.

  if (solState_.twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM)
  {
    begin = nonPdeInstancePtrVec_.begin();
    end   = nonPdeInstancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->loadTrivialDAE_FMatrixStamp ();
    }

    begin = pdeInstancePtrVec_.begin();
    end   = pdeInstancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      tmpBool = (*iter)->loadDAEdQdx();bsuccess = bsuccess && tmpBool;
      tmpBool = (*iter)->loadDAEdFdx();bsuccess = bsuccess && tmpBool;
    }
  }
  // Else, do a normal analytical matrix load.
  else
  {
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->loadDAEMatrices 
        (*(externData_.dFdxMatrixPtr) , *(externData_.dQdxMatrixPtr));
    }
  }

  // Run jacobian diagnostic.
  if (devOptions_.testJacobianFlag &&
      (solState_.timeStepNumber >= devOptions_.testJacStartStep &&
       solState_.timeStepNumber <= devOptions_.testJacStopStep)
     )
  {
    if (DEBUG_DEVICE)
    {
      devOptions_.debugLevel -= 2;
    }

    // Test just the specified device(s), if the user specified any.
    if ((devOptions_.testJacDeviceNameGiven))
    {
      begin = testJacDevicePtrVec_.begin();
      end   = testJacDevicePtrVec_.end();
    }
    else // Test all the devices:
    {
      begin = instancePtrVec_.begin();
      end   = instancePtrVec_.end();
    }

    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->testDAEMatrices (nameVec_);
    }

    if (DEBUG_DEVICE)
      devOptions_.debugLevel += 2;
  }

  //Tell Jacobian, fill is complete allowing accumulation if necessary
  externData_.dQdxMatrixPtr->fillComplete();
  externData_.dFdxMatrixPtr->fillComplete();

  Report::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

  if (DEBUG_DEVICE && devOptions_.debugLevel > 1 && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    dout() << section_divider << std::endl;
    dout() <<  "Q-matrix: nonlinear iteration = " << newtonIter << "\n";
    externData_.dQdxMatrixPtr->printPetraObject(dout());
    dout() << std::endl;
    dout() << section_divider << std::endl;
    dout() <<  "F-matrix: nonlinear iteration = " << newtonIter << "\n";
    externData_.dFdxMatrixPtr->printPetraObject(dout());
    dout() << std::endl;
    dout() << section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEVectors
// Purpose       : This function loads the various DAE related vectors to
//                 set up the following expression:
//
//                 f(x) = dQ/dt + F(x) - B(t) = 0
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadRHSVector.
//
//                 Note that this function, unlike the loadRHSVector
//                 function *requires* that vectors be passed in.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load vectors
//                 are zeroed out.  That is, F, Q, B, dFdxdVp and dQdxVp
//                 are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEVectors
  (N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpCurrSolVectorPtr,
   N_LAS_Vector * tmpLastSolVectorPtr,
   N_LAS_Vector * tmpStaVectorPtr,
   N_LAS_Vector * tmpCurrStaVectorPtr,
   N_LAS_Vector * tmpLastStaVectorPtr,
   N_LAS_Vector * tmpStaDerivVectorPtr,
   N_LAS_Vector * tmpStoVectorPtr,
   N_LAS_Vector * tmpCurrStoVectorPtr,
   N_LAS_Vector * tmpLastStoVectorPtr,
   N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
   N_LAS_Vector * tmpQVectorPtr,
   N_LAS_Vector * tmpFVectorPtr,
   N_LAS_Vector * tmpBVectorPtr,
   N_LAS_Vector * tmpdFdxdVpVectorPtr,
   N_LAS_Vector * tmpdQdxdVpVectorPtr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = tmpSolVectorPtr;
  externData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  externData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  externData_.daeQVectorPtr    = tmpQVectorPtr;
  externData_.daeFVectorPtr    = tmpFVectorPtr;
  externData_.daeBVectorPtr    = tmpBVectorPtr;
  externData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  externData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  externData_.nextStaVectorPtr = tmpStaVectorPtr;
  externData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  externData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  externData_.storeLeadCurrQCompPtr = tmpStoLeadCurrQCompVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoVectorPtr;
  externData_.currStoVectorPtr = tmpCurrStoVectorPtr;
  externData_.lastStoVectorPtr = tmpLastStoVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
  externData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // Set up the relevant RAW Pointers:
  setupRawVectorPointers_ ();

  // call all the intermediate vars loads:
  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM)
  {
    begin = pdeDevicePtrVec_.begin ();
    end   = pdeDevicePtrVec_.end ();
  }
  else
  {
    begin = devicePtrVec_.begin ();
    end   = devicePtrVec_.end ();
  }

#ifndef Xyce_EXCLUDE_SECONDARY_STATE
  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->updateSecondaryState 
      (externData_.nextStaDerivVectorRawPtr, externData_.nextStoVectorRawPtr);
    bsuccess = bsuccess && tmpBool;
  }
#endif // Xyce_EXCLUDE_SECONDARY_STATE

#ifdef Xyce_PARALLEL_MPI
  externData_.nextStaVectorPtr->importOverlap();
  externData_.nextStoVectorPtr->importOverlap();
#endif

  if (solState_.twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM)
  {
    int numDevices = pdeDevicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=pdeDevicePtrVec_.at(i)->loadDAEVectors(externData_.nextSolVectorRawPtr,
                                                       externData_.daeFVectorRawPtr,
                                                       externData_.daeQVectorRawPtr,
                                                       externData_.daeBVectorRawPtr,
                                                       externData_.nextStoVectorRawPtr,
                                                       externData_.storeLeadCurrQCompRawPtr);
    }
  }
  else
  {
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->loadDAEVectors(externData_.nextSolVectorRawPtr,
                                                    externData_.daeFVectorRawPtr,
                                                    externData_.daeQVectorRawPtr,
                                                    externData_.daeBVectorRawPtr,
                                                    externData_.nextStoVectorRawPtr,
                                                    externData_.storeLeadCurrQCompRawPtr);
    }
  }

  // dump to the screen:
  if (DEBUG_DEVICE && devOptions_.debugLevel > 1 && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    dout() <<  "Q-vector: nonlinear iteration = " << newtonIter << "\n";
    externData_.daeQVectorPtr->printPetraObject(std::cout);
    dout() << std::endl;
    dout() <<  "F-vector: nonlinear iteration = " << newtonIter << "\n";
    externData_.daeFVectorPtr->printPetraObject(std::cout);
    dout() << std::endl;

    if (devOptions_.voltageLimiterFlag)
    {
      dout() << "\n\n  dFdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dFdxdVpVectorPtr->printPetraObject(std::cout);
      dout() << std::endl;
      dout() << "\n\n  dQdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dQdxdVpVectorPtr->printPetraObject(std::cout);
      dout() << std::endl;
    }

    debugOutput2();
  }

  // Update parallel if necessary
  externData_.daeQVectorPtr->fillComplete();
  externData_.daeFVectorPtr->fillComplete();
  externData_.dFdxdVpVectorPtr->fillComplete();
  externData_.dQdxdVpVectorPtr->fillComplete();

  Report::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

#ifdef Xyce_SIZEOF
  int size = sizeof(*this);
  dout() << "Size of device package after vector load = " << size << std::endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDeviceMask ()
// Purpose       : let devices set elements of a mask telling the time
//                 integrator what equations should be ignored in taking
//                 weighted norms for error control purposes.
// Special Notes : Devices should *only* zero internal variables, and then
//                 only those that absolutely should never be used to
//                 control step size (e.g. excess phase variables in BJTs)
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/18/07
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDeviceMask()
{

  // call all the intermediate vars loads:
  InstanceVector::iterator iter;
  InstanceVector::iterator begin;
  InstanceVector::iterator end;

  begin = instancePtrVec_.begin ();
  end   = instancePtrVec_.end ();

  nonTrivialDeviceMaskFlag = false;

  for (iter=begin; iter!=end;++iter)
  {
    nonTrivialDeviceMaskFlag |= (*iter)->loadDeviceMask();
  }

  externData_.deviceMaskVectorPtr->fillComplete();

  // make sure the system's flag reflects ours:
  externData_.lasSysPtr->setNonTrivialDeviceMaskFlag(nonTrivialDeviceMaskFlag);

  return nonTrivialDeviceMaskFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalPar ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/05
//-----------------------------------------------------------------------------
void DeviceMgr::addGlobalPar (const Util::Param & par)
{
  if (par.getType() == Util::EXPR)
  {
    globals_.global_expressions.push_back(par.getValue<Util::Expression>());
    globals_.global_exp_names.push_back(par.uTag());

    Util::Expression &expression = globals_.global_expressions.back();

    std::vector<std::string> variables;
    expression.get_names(XEXP_VARIABLE, variables);
    std::vector<std::string> specials;
    expression.get_names(XEXP_SPECIAL, specials);
    if (!specials.empty())
    {
      expression.set_sim_time(solState_.currTime);
    }

    for (std::vector<std::string>::const_iterator it = variables.begin(), 
        end = variables.end() ; it != end ; ++it)
    {
      expression.set_var(*it, globals_.global_params[*it]);
    }

    double val;
    expression.evaluateFunction(val);
    globals_.global_params[par.uTag()] = val;
  }
  else
  {
    globals_.global_params[par.uTag()] = par.getImmutableValue<double>();
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getGlobalPar ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schie, SNL, Electrical Systems Modeling
// Creation Date : 01/25/13
//-----------------------------------------------------------------------------
double DeviceMgr::getGlobalPar (const std::string & parName) const
{
  double retVal = 0;
  GlobalParameterMap::const_iterator parLocItr = globals_.global_params.find(parName);
  if (parLocItr != globals_.global_params.end())
  {
    // extract the value for return
    retVal = parLocItr->second;
  }
  else
  {
    Report::UserError() << "Could not find global parameter \"" << parName << "\"";
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::findGlobalPar
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const double *DeviceMgr::findGlobalPar (const std::string & parName) const
{
  GlobalParameterMap::const_iterator it = globals_.global_params.find(parName);
  if (it != globals_.global_params.end())
    return &(*it).second;

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateIntermediateVars_
// Purpose       : This function calls updateIntermediateVars
//                 for the current devices.
// Special Notes :
// Scope         : private
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 12/14/2006
//-----------------------------------------------------------------------------
bool DeviceMgr::updateIntermediateVars_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateIntermediateVars ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updatePrimaryState_
// Purpose       : This function updates primary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updatePrimaryState_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updatePrimaryState ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSecondaryState_
// Purpose       : This function function updates secondary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSecondaryState_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateSecondaryState ();
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : DeviceMgr::updateDependentParameters_
// Purpose        : This function updates all dependent parameters for
//                  the current time step.
// Special Notes  : This was evolved from updateTimeDependentParameters_
// Scope          : private
// Creator        : Dave Shirley
// Creation Date  : 08/17/06
//----------------------------------------------------------------------------
bool DeviceMgr::updateDependentParameters_()
{
  bool bsuccess = true;
  N_LAS_Vector * solVectorPtr = externData_.nextSolVectorPtr;

  GlobalParameterMap & gp = globals_.global_params;
  std::vector<Util::Expression> & ge = globals_.global_expressions;

  if (timeParamsProcessed_ != solState_.currTime)
    parameterChanged_ = true;

  // Update global params for new time and other global params
  int pos = 0;
  for (std::vector<Util::Expression>::iterator g_i = ge.begin(), g_end = ge.end(); 
      g_i != g_end; ++g_i)
  {
    bool changed = false;
    if (g_i->set_sim_time(solState_.currTime))
      changed = true;

    std::vector<std::string> variables;
    g_i->get_names(XEXP_VARIABLE, variables);
    for (std::vector<std::string>::iterator vs_i = variables.begin(), vs_end = variables.end(); 
        vs_i != vs_end; ++vs_i)
    {
      if (g_i->set_var(*vs_i, gp[*vs_i]))
        changed = true;
    }

    if (changed)
    {
      double val;

      parameterChanged_ = true;
      g_i->evaluateFunction(val);
      gp[globals_.global_exp_names[pos]] = val;
    }
    ++pos;
  }

  // do the models:
  if (firstDependent_)
  {
    firstDependent_ = false;

    dependentPtrVec_.clear();

    ModelVector::iterator iterM;
    ModelVector::iterator beginM =modelVector_.begin();
    ModelVector::iterator endM =modelVector_.end();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      if (!(*iterM)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM));
        bool tmpBool = (*iterM)->updateGlobalParameters(gp);
        bsuccess = bsuccess && tmpBool;
        tmpBool = (*iterM)->updateDependentParameters (*solVectorPtr);
        bsuccess = bsuccess && tmpBool;
        (*iterM)->processParams();
        (*iterM)->processInstanceParams();
      }
    }

    // do the instances
    InstanceVector::iterator iter;
    InstanceVector::iterator begin =instancePtrVec_.begin();
    InstanceVector::iterator end =instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      if (!(*iter)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iter));
        bool tmpBool = (*iter)->updateGlobalParameters(gp);
        bsuccess = bsuccess && tmpBool;
        tmpBool = (*iter)->updateDependentParameters (*solVectorPtr);
        bsuccess = bsuccess && tmpBool;
        (*iter)->processParams();
      }
    }
  }
  else
  {
    EntityVector::iterator iter;
    EntityVector::iterator begin = dependentPtrVec_.begin();
    EntityVector::iterator end = dependentPtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      bool changed = false;
      if (parameterChanged_)
      {
        bool tmpBool = (*iter)->updateGlobalParameters(gp);
        changed = changed || tmpBool;
        bsuccess = bsuccess && tmpBool;
      }
      bool tmpBool = (*iter)->updateDependentParameters (*solVectorPtr);
      changed = changed || tmpBool;
      bsuccess = bsuccess && tmpBool;
      if (changed)
      {
        (*iter)->processParams();
        (*iter)->processInstanceParams();
      }
    }
  }
  timeParamsProcessed_ = solState_.currTime;
  parameterChanged_ = false;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadBVectorsforAC
// Purpose       : This function loads the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::loadBVectorsforAC(N_LAS_Vector * bVecRealPtr, N_LAS_Vector * bVecImagPtr)
{
  bool bsuccess = true;

// copy over the passed pointers:
  externData_.bVecRealPtr = bVecRealPtr;
  externData_.bVecImagPtr = bVecImagPtr;

#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
#endif
  // Now reset the relevant RAW pointers:
  externData_.bVecRealRawPtr = &((*externData_.bVecRealPtr)[0]);
  externData_.bVecImagRawPtr = &((*externData_.bVecImagPtr)[0]);

  std::vector<SourceInstance*>::iterator vIter;
  std::vector<SourceInstance*>::iterator vBegin =indepSourceInstancePtrVec_.begin();
  std::vector<SourceInstance*>::iterator vEnd =indepSourceInstancePtrVec_.end();
  for (vIter=vBegin; vIter!=vEnd;++vIter)
  {
    (*vIter)->loadBVectorsforAC(externData_.bVecRealRawPtr, externData_.bVecImagRawPtr);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBMatrixEntriesforMOR()
// Purpose       : This function obtains the indices for the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec,
                                        std::vector<int>& bMatPosEntriesVec)
{
  bool bsuccess = true;

  int lpos, lneg, lbra;
  std::vector<SourceInstance *>::iterator vIter;
  std::vector<SourceInstance *>::iterator vBegin = indepSourceInstancePtrVec_.begin();
  std::vector<SourceInstance *>::iterator vEnd = indepSourceInstancePtrVec_.end();

  for (vIter=vBegin; vIter!=vEnd;++vIter)
  {
    Vsrc::Instance* vsrc = dynamic_cast<Vsrc::Instance *>(*vIter);
     if (vsrc != 0)
     {
       vsrc->getLIDs(lpos, lneg, lbra);
       bMatEntriesVec.push_back(lbra);
       bMatPosEntriesVec.push_back(lpos);
     }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSources
// Purpose       : This function function updates sources for the present
//                 time step.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSources()
{
  bool bsuccess = true;

  setupSolverInfo_ ();

  std::vector<SourceInstance*>::iterator vIter;
  std::vector<SourceInstance*>::iterator vBegin = indepSourceInstancePtrVec_.begin();
  std::vector<SourceInstance*>::iterator vEnd = indepSourceInstancePtrVec_.end();
  for (vIter=vBegin; vIter!=vEnd; ++vIter)
  {
    (*vIter)->updateSource();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setICs
// Purpose       : This function function sets initial conditions for devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/13/00
//-----------------------------------------------------------------------------
bool DeviceMgr::setICs(
   N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpCurrSolVectorPtr,
   N_LAS_Vector * tmpLastSolVectorPtr,
   N_LAS_Vector * tmpStaVectorPtr,
   N_LAS_Vector * tmpCurrStaVectorPtr,
   N_LAS_Vector * tmpLastStaVectorPtr,
   N_LAS_Vector * tmpStaDerivVectorPtr,
   N_LAS_Vector * tmpStoVectorPtr,
   N_LAS_Vector * tmpCurrStoVectorPtr,
   N_LAS_Vector * tmpLastStoVectorPtr,
   N_LAS_Vector * tmpQVectorPtr,
   N_LAS_Vector * tmpFVectorPtr,
   N_LAS_Vector * tmpBVectorPtr,
   N_LAS_Vector * tmpdFdxdVpVectorPtr,
   N_LAS_Vector * tmpdQdxdVpVectorPtr)
{
  bool bsuccess = true;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = tmpSolVectorPtr;
  externData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  externData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  externData_.daeQVectorPtr    = tmpQVectorPtr;
  externData_.daeFVectorPtr    = tmpFVectorPtr;
  externData_.daeBVectorPtr    = tmpBVectorPtr;
  externData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  externData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  externData_.nextStaVectorPtr = tmpStaVectorPtr;
  externData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  externData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoVectorPtr;
  externData_.currStoVectorPtr = tmpCurrStoVectorPtr;
  externData_.lastStoVectorPtr = tmpLastStoVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
  externData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // if IC's on devices are set, we need to ensure that the
  // raw pointers are up to date first.
  setupRawVectorPointers_ ();

  for (InstanceVector::iterator iter = instancePtrVec_.begin(); 
      iter != instancePtrVec_.end(); ++iter)
  {
    (*iter)->setIC();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool DeviceMgr::output()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin = plotFileInstancePtrVec_.begin ();
  InstanceVector::iterator end = plotFileInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    bool tmpBool = (*iter)->outputPlotFiles ();
    bsuccess = bsuccess && tmpBool;
  }

  // only output .OP information once.
  if (!dotOpOutputFlag && solState_.OPspecified)
  {
    dotOpOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::finishOutput
// Purpose       : Same as output, only this one forces the output.
//
// Special Notes : This function was neccessary with the implementation of
//                 outputInterval.  The final output, which needs to
//                 happen after the transient is over, won't neccessarily
//                 happen with outputInterval enabled.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/19/04
//-----------------------------------------------------------------------------
bool DeviceMgr::finishOutput ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  solState_.forceFinalOutput = true;
  tmpBool = output ();
  bsuccess = bsuccess && tmpBool;
  solState_.forceFinalOutput = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dotOpOutput
// Purpose       :
// Special Notes : This is a quick-and-dirty implementation, to get something
//                 working quickly.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/3/12
//-----------------------------------------------------------------------------
void DeviceMgr::dotOpOutput ()
{
  dotOpOutputFlag = true;

  std::map<std::string, Device *> device_map;
  for (EntityTypeIdDeviceMap::const_iterator it = deviceMap_.begin(); 
      it != deviceMap_.end(); ++it)
  {
    device_map[(*it).second->getName()] = (*it).second;
  }

  lout() << section_divider << "\n"
         << "Operating point information:";

  for (std::map<std::string, Device *>::const_iterator it = device_map.begin(); 
      it != device_map.end(); ++it)
  {
    Xyce::Device::printDotOpOutput(lout(), *(*it).second);
  }

  lout() << section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setGlobalFlags
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/2/12
//-----------------------------------------------------------------------------
void DeviceMgr::setGlobalFlags()
{
  int pde_flag = solState_.PDESystemFlag ? 1 : 0;

  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_LOR, &pde_flag, 1);

  solState_.PDESystemFlag = pde_flag != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::getBreakPoints (std::vector<Util::BreakPoint> & breakPointTimes)
{
  InstanceVector::iterator iterI;
  ModelVector::iterator iterM;
  bool bsuccess = true;
  bool tmpBool = true;

  tmpBool = setupSolverInfo_();
  bsuccess = bsuccess && tmpBool;
  setupRawVectorPointers_ ();

  // For some devices we only need to set breakpoints caused by discontinuities
  // in their parameters:

  for (iterM = modelVector_.begin() ; iterM != modelVector_.end() ; ++iterM)
  {
    if (!(*iterM)->getDependentParams().empty())
    {
      tmpBool = (*iterM)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  for (iterI = instancePtrVec_.begin() ; iterI != instancePtrVec_.end() ; ++iterI)
  {
    if (!(*iterI)->getDependentParams().empty())
    {
      tmpBool = (*iterI)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  // Breakpoints for global params:
  std::vector<Util::Expression>::iterator globalExp_i =
    globals_.global_expressions.begin();
  std::vector<Util::Expression>::iterator globalExp_end =
    globals_.global_expressions.end();
  for (; globalExp_i != globalExp_end; ++globalExp_i)
  {
    double bTime = globalExp_i->get_break_time();
    if (bTime > solState_.currTime)
      breakPointTimes.push_back(bTime);
  }

  if (!breakPointInstancesInitialized)
  {
    InstanceVector::iterator beginI = instancePtrVec_.begin ();
    InstanceVector::iterator endI = instancePtrVec_.end ();

    for (iterI=beginI;iterI!=endI;++iterI)
    {
      // this function returns false if it is the base class, and true otherwise.
      bool functionSetup = (*iterI)->getInstanceBreakPoints (breakPointTimes);
      if (functionSetup)
      {
        bpInstancePtrVec_.push_back(*iterI);
      }
    }
    breakPointInstancesInitialized = true;
  }
  else
  {
    InstanceVector::iterator beginI = bpInstancePtrVec_.begin ();
    InstanceVector::iterator endI = bpInstancePtrVec_.end ();
    for (iterI=beginI;iterI!=endI;++iterI)
    {
      bool functionSetup = (*iterI)->getInstanceBreakPoints (breakPointTimes);
    }
  }

  ModelTypeInstanceVectorMap::const_iterator model_type_it = 
    modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) 
  {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); 
        it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.getBreakPoints(breakPointTimes);
    }
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupSolverInfo_
//
// Purpose       :
//
// Special Notes : This function gets called a lot, and this can be kind of
//                 confusing.  Probably, this function is being used to handle
//                 too many different types of data.
//
//                 For example, it gets called at the beginning
//                 of the "updateSources" function.  Why?  Because the sources
//                 need to know which step we are at, and/or what the current
//                 time is, to do their update properly.
//
//                 However, at the moment, updateSources also provides
//                 information that setupSolverInfo needs.  Only after the
//                 sources have been updated do we know if a sweep source has
//                 been reset.  And, the sweepSourceResetFlag  is used by
//                 setupSolverInfo, to set up the initJctFlag boolean.
//                 So it will need to be called at least one more
//                 time, at the beginning of the RHS load.  (which it is).
//
//                 Anyway, any of the functions that are called from the
//                 outside, such as:  updateSources, loadRHSVector,
//                 loadJacobianMatrix, etc.... have no way of knowing when,
//                 w.r.t. the solvers they are being called.  The only way
//                 to do this properly is to have each of them request
//                 the current solver state, before they go to do their
//                 work.  Hence, this function is called a lot.
//
//                 Unfortunately, this has led to a somewhat sloppy and
//                 confusing interface between the solvers and the
//                 device package.  I wanted to avoid having a lot of
//                 function arguments being passed around for each of
//                 these functions, in part because the calling code
//                 (NLS) doesn't know everything. - NLS knows about the newton
//                 step, but it doesn't know the time step, for example.
//
//                 At some point I hope to refactor this.
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/30/01
//-----------------------------------------------------------------------------
bool DeviceMgr::setupSolverInfo_()
{
  bool bsuccess = true;

  N_TIA_TimeIntInfo & tiInfo = solState_.tiInfo;
  N_NLS_NonLinInfo  & nlInfo = solState_.nlInfo;

  // Time integration info:
  anaIntPtr_->getTimeIntInfo(tiInfo);
  solState_.pdt                  = tiInfo.pdt;
  solState_.currTimeStep         = tiInfo.nextTimeStep;
  solState_.lastTimeStep         = tiInfo.currTimeStep;
  solState_.currTime             = tiInfo.nextTime;
  solState_.finalTime            = tiInfo.finalTime;
  solState_.startingTimeStep     = tiInfo.startingTimeStep;
  solState_.bpTol                = tiInfo.bpTol;
  solState_.currentOrder         = tiInfo.currentOrder; // BNB
  solState_.usedOrder            = tiInfo.usedOrder;  // BNB

  if (solState_.mpdeOnFlag == 1)
  {
    solState_.dcopFlag = 0;
    solState_.initTranFlag =  1;
    solState_.beginIntegrationFlag =  1;
  }
  else
  {
    solState_.dcopFlag             = tiInfo.dcopFlag;
    solState_.initTranFlag         = tiInfo.initTranFlag;
    solState_.beginIntegrationFlag = tiInfo.beginIntegrationFlag;
  }

  solState_.inputOPFlag          = tiInfo.inputOPFlag;
  solState_.acopFlag             = tiInfo.acopFlag;
  solState_.tranopFlag           = tiInfo.tranopFlag;
  solState_.transientFlag        = tiInfo.transientFlag;
  solState_.dcsweepFlag          = tiInfo.dcsweepFlag;
  solState_.sweepSourceResetFlag = tiInfo.sweepSourceResetFlag;

  solState_.timeStepNumber       = tiInfo.timeStepNumber;
//  solState_.initTranFlag         = tiInfo.initTranFlag;
//  solState_.beginIntegrationFlag = tiInfo.beginIntegrationFlag;

  solState_.doubleDCOPStep       = tiInfo.doubleDCOPStep;
  solState_.doubleDCOPEnabled    = tiInfo.doubleDCOPEnabled;
  solState_.stepLoopIter         = tiInfo.stepLoopIter;

  // Nonlinear solver info:
  nlsMgrPtr_->getNonLinInfo(nlInfo);
  solState_.newtonIter           = nlInfo.newtonIter;
  solState_.twoLevelNewtonCouplingMode         = nlInfo.twoLevelNewtonCouplingMode;

  // Get LOCA-specific information.  Note - in general, LOCA is only used for
  // steady state calculations - there isn't much point in using it
  // for transient.  A common case is one where LOCA is used for the
  // tranOP, but not for the subsequent transient phase.  The
  // locaEnabledFlag should switch from true to false under
  // that scenario, once the transient phase starts.
  solState_.locaEnabledFlag      = nlInfo.locaFlag;
  if (solState_.locaEnabledFlag) // if no LOCA, these are 0, true, respectively.
  {
    solState_.continuationStepNumber = nlInfo.continuationStep;
    solState_.firstContinuationParam = nlInfo.firstContinuationParam;
    solState_.firstSolveComplete     = nlInfo.firstSolveComplete;
  }
  // Done with LOCA information.

  // Setup the initialize junctions flag.
  // The initJct flag should only be true if we are at the first Newton step of
  // an initial point in the calculation.  Examples include:
  //     - 1st Newton step of the DCOP initialization for transient (tranOp)
  //     - 1st Newton step of first DC sweep step.
  //     - 1st Newton step of a sweep that has been reset.  That typically
  //         happens if the sweep is multi-dimensional, and the inner loop has
  //         cycled back to the beginning again.
  //
  bool resetFlag =  (solState_.timeStepNumber==0) || (solState_.sweepSourceResetFlag);

  // Do this if using LOCA for a DC or tranop calculation.
  if (solState_.dcopFlag && solState_.locaEnabledFlag)
  {
    resetFlag = resetFlag && (solState_.continuationStepNumber==0);
  }

  solState_.initJctFlag = ((solState_.dcopFlag) &&
                            (solState_.newtonIter==0) &&
                             solState_.firstContinuationParam &&
                             !solState_.firstSolveComplete && resetFlag);

  // initFixFlag: try to mimic "MODEINITFIX" of SPICE.  This is set if:
  //   DCOP or TranOP
  //   Not first iteration
  //   Any device not converged

  solState_.initFixFlag = ((solState_.dcopFlag) &&
                            !(allDevsConverged()) &&
                            (solState_.newtonIter!=0) &&
                             solState_.firstContinuationParam &&
                             !solState_.firstSolveComplete && resetFlag);

  if (solState_.dcopFlag)
  {
    solState_.ltraTimeIndex = 0;
    solState_.ltraTimeHistorySize = 10;
    solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
  }

  // One final check.  See if the "external state" has been set.  If so,
  // check to see if it has the initJctFlag set.  If not, then we probably
  // shouldn't either.  The external state comes from a higher up level
  // in a multi-level newton solve.
  //
  // This should be made more detailed later.
  if (externalStateFlag_)
  {
    if (solState_.newtonIter==0 && solState_.dcopFlag)
    {
      solState_.initJctFlag = solStateExternal_.initJctFlag;
    }
  }

  // The first DCOP step of a "double DCOP" simulation is a special case,
  // in which the nonlinear poisson is solved in place of drift-diffusion
  // equations for the PDE devices.  For this initialization problem, the
  // circuit subproblem should not be included.
  if ((solState_.doubleDCOPEnabled) &&
      (solState_.dcopFlag) &&
      (solState_.doubleDCOPStep == 0))
  {
    solState_.twoLevelNewtonCouplingMode       = Nonlinear::INNER_PROBLEM;
  }

  // if necessary, set up the names vector
  if (devOptions_.testJacobianFlag)
  {
    const NodeNamePairMap & nodeNames = outputMgrPtr_->getAllNodes();
    int nodeNameSize = nodeNames.size();
    nameVec_.resize(nodeNameSize+1,"gnd");
    NodeNamePairMap::const_iterator mapI, mapEnd;
    mapEnd = nodeNames.end();
    mapI =  nodeNames.begin();
    for (; mapI != mapEnd ; ++mapI)
    {
      nameVec_[(*mapI).second.first] = (*mapI).first;
    }
  }

  if (DEBUG_DEVICE) {
    if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
    {
      dout() << solState_;
    }

    solState_.debugTimeFlag =
      (solState_.currTime       >= devOptions_.debugMinTime &&
       solState_.currTime       <= devOptions_.debugMaxTime) &&
      (solState_.timeStepNumber >= devOptions_.debugMinTimestep &&
       solState_.timeStepNumber <= devOptions_.debugMaxTimestep);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawVectorPointers_
// Purpose       : set up raw pointers
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawVectorPointers_ ()
{
  if (externData_.daeQVectorPtr != 0)
  {
    externData_.daeQVectorRawPtr    = &((*externData_.daeQVectorPtr)[0]);
  }

  if (externData_.daeFVectorPtr != 0)
  {
    externData_.daeFVectorRawPtr    = &((*externData_.daeFVectorPtr)[0]);
  }

  if (externData_.daeBVectorPtr != 0)
  {
    externData_.daeBVectorRawPtr    = &((*externData_.daeBVectorPtr)[0]);
  }

  if (externData_.dFdxdVpVectorPtr != 0)
  {
    externData_.dFdxdVpVectorRawPtr = &((*externData_.dFdxdVpVectorPtr)[0]);
  }

  if (externData_.dQdxdVpVectorPtr != 0)
  {
    externData_.dQdxdVpVectorRawPtr = &((*externData_.dQdxdVpVectorPtr)[0]);
  }

  if (externData_.nextSolVectorPtr != 0)
  {
    externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  }

  if (externData_.currSolVectorPtr != 0)
  {
    externData_.currSolVectorRawPtr = &((*externData_.currSolVectorPtr)[0]);
  }

  if (externData_.lastSolVectorPtr != 0)
  {
    externData_.lastSolVectorRawPtr = &((*externData_.lastSolVectorPtr)[0]);
  }

  if (externData_.nextStaVectorPtr != 0)
  {
    externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  }

  if (externData_.currStaVectorPtr != 0)
  {
    externData_.currStaVectorRawPtr = &((*externData_.currStaVectorPtr)[0]);
  }

  if (externData_.lastStaVectorPtr != 0)
  {
    externData_.lastStaVectorRawPtr = &((*externData_.lastStaVectorPtr)[0]);
  }

  if (externData_.nextStaDerivVectorPtr != 0)
  {
    externData_.nextStaDerivVectorRawPtr = &((*externData_.nextStaDerivVectorPtr)[0]);
  }

  if (externData_.nextStoVectorPtr != 0)
  {
    externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);
  }

  if (externData_.currStoVectorPtr != 0)
  {
    externData_.currStoVectorRawPtr = &((*externData_.currStoVectorPtr)[0]);
  }

  if (externData_.lastStoVectorPtr != 0)
  {
    externData_.lastStoVectorRawPtr = &((*externData_.lastStoVectorPtr)[0]);
  }

  if (externData_.storeLeadCurrQCompPtr != 0)
  {
    externData_.storeLeadCurrQCompRawPtr = &((*externData_.storeLeadCurrQCompPtr)[0]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawMatrixPointers_
// Purpose       : set up raw pointers for matrices
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/03/09
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawMatrixPointers_ ()
{
    InstanceVector::iterator iter;
    InstanceVector::iterator begin;
    InstanceVector::iterator end;
    begin = instancePtrVec_.begin();
    end   = instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->setupPointers();
    }
    return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
double DeviceMgr::getMaxTimeStepSize ()
{
  double maxStep = devOptions_.defaultMaxTimeStep;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    double step = (*iter)->getMaxTimeStepSize ();
    SourceInstance * srcInst = dynamic_cast<SourceInstance*>(*iter);
    if (!srcInst || !srcInst->getFastSourceFlag())
      maxStep = Xycemin(step, maxStep);
  }

  return maxStep;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::declareCurrentStepAsBreakpoint
// Purpose       : If during a device load, a device must act in a discontinuous
//                 fashion, let the analysis manager know that this step should
//                 be treated as such.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/5/2010
//-----------------------------------------------------------------------------
void DeviceMgr::declareCurrentStepAsBreakpoint()
{
  anaIntPtr_->setBeginningIntegrationFlag(true);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::enablePDEContinuation
// Purpose       : This function turns on the Continuation flag, which lets
//                 the devices which are Continuation enabled know that they
//                 need to set up their variable parameters.
//
// Special Notes : Currently, only PDE devices can take advantage of this
//                 capability, and it is only used in the context of two-level
//                 Newton.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int DeviceMgr::enablePDEContinuation ()
{
  if (DEBUG_DEVICE)
    dout() << "DeviceMgr::enablePDEContinuation" << std::endl;

  bool bsuccess = true;
  solState_.PDEcontinuationFlag  = true;
  solState_.maxPDEContinuationSteps = 1;
  solState_.currPDEContinuationStep = 0;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();

  for (iter=begin; iter!=end;++iter)
  {
    bool tmpSuccess = (*iter)->enablePDEContinuation();
    bsuccess = bsuccess && tmpSuccess;
  }

  // if any of the devices feels that it needs more than the specified
  // number of continuation steps, re-do, with the new step number.
  if (solState_.maxPDEContinuationSteps != 1)
  {
    for (iter=begin; iter!=end;++iter)
    {
      bool tmpSuccess = (*iter)->enablePDEContinuation();
      bsuccess = bsuccess && tmpSuccess;
    }
  }

  int returnedSteps;

  if (!bsuccess)  returnedSteps = -1;
  else            returnedSteps = solState_.maxPDEContinuationSteps;

  return returnedSteps;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool DeviceMgr::disablePDEContinuation ()
{
  bool bsuccess = true;
  solState_.PDEcontinuationFlag = false;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    bool tmpSuccess = (*iter)->disablePDEContinuation();
    bsuccess = bsuccess && tmpSuccess;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcPDESubProblemInfo
//
// Purpose       : Determines the number of PDE sub-problems.,
//
//                 This is mainly used/needed for 2-level problems.
//
//                 Also determines the number of interface nodes per
//                 sub-problem (ie number of electrodes on each device).
//
// Special Notes : Need to modify to work correctly in parallel, probably.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcPDESubProblemInfo()
{
  // now set up numInterfaceNodes_;
  numInterfaceNodes_.reserve(pdeInstancePtrVec_.size());

  for (InstanceVector::const_iterator it = pdeInstancePtrVec_.begin(), end = pdeInstancePtrVec_.end(); it != end; ++it)
  {
    numInterfaceNodes_.push_back((*it)->getNumExtVars());
  }

  calledBeforeCSPI = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumInterfaceNodes
// Purpose       : returns the vector calculaed in calcPDESubProblemInfo.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
void DeviceMgr::getNumInterfaceNodes (std::vector<int> & numINodes)
{
  if (!calledBeforeCSPI)
  {
    calcPDESubProblemInfo ();
  }

  int size = numINodes.size ();
  int size2 = numInterfaceNodes_.size();

  if (size < size2) numINodes.resize(size2);

  for (int i=0;i<size2;++i)
  {
    numINodes[i] = numInterfaceNodes_[i];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->loadDFDV(iElectrode,dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcCouplingTerms (int iPDEDevice, int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->calcConductance(iElectrode, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::raiseDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/23/03
//-----------------------------------------------------------------------------
bool DeviceMgr::raiseDebugLevel(int increment)
{
  if (DEBUG_DEVICE)
    devOptions_.debugLevel += increment;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getHomotopyBlockSize
// Purpose       : Returns the number of mosfet gainscale blocks (a value
//                 in the device options).
// Special Notes : Needed for block homotopy on gainscale.
// Scope         : public
// Creator       : Roger P. Pawlowski, SNL, Computational Sciences
// Creation Date : 01/26/2005
//-----------------------------------------------------------------------------
int DeviceMgr::getHomotopyBlockSize() const
{
  return devOptions_.numGainScaleBlocks;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/04
//-----------------------------------------------------------------------------
bool DeviceMgr::updateTemperature (double val)
{
  bool bsuccess = true, success = true;

  // convert to kelvin:
  double Ctemp = val;
  double Ktemp = val + CONSTCtoK;

  if (DEBUG_DEVICE && devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
    dout() << "In DeviceMgr::updateTemperature.  new C temp = " << Ctemp << " K temp = " << Ktemp << std::endl;

  // First set the global temp.  This is used in each device if the "tempGiven"
  // variable is false.  This should be in Kelvin.
  devOptions_.temp.setVal(Ktemp);

  {
    // loop over the bsim3 models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B3::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  {
    // loop over the bsim4 models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B4::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  {
    // loop over the b3soi models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B3SOI::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  // Loop over all models, call processParams, with CTemp.
  // "XYCEADMS*TEMP is there to force Verilog devices, which might have
  // temperature dependence through "$temperature" instead of a "TEMP"
  // parameter, to work properly.  If so, they need the temperature set in
  // Kelvin.
  std::string tname("TEMP");
  std::string tname2("XYCEADMSMODTEMP");
  std::string tname3("XYCEADMSINSTTEMP");
  for (ModelVector::const_iterator it = modelVector_.begin(); it != modelVector_.end(); ++it)
  {
    success = (*it)->setParam (tname, Ctemp);
    success = (*it)->setParam (tname2, Ktemp) || success;
    success = success && (*it)->processParams ();
  }

  // Loop over device instances, and set the temperature.  This should be
  // in C, if going through processParams, and K if going through
  // the updateTemperature function.
  InstanceVector::const_iterator iter;
  InstanceVector::const_iterator begin;
  InstanceVector::const_iterator end;

  begin = instancePtrVec_.begin();
  end   = instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    success = (*iter)->setParam (tname, Ctemp);
    success = (*iter)->setParam (tname3, Ktemp) || success;
    success = success && (*iter)->processParams ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
bool DeviceMgr::allDevsConverged()
{
  int allDevsConv = true;

  // if two-level, and just the inner problem, only check the PDE devices,
  // coz those are all we loaded.
  InstanceVector::iterator iter;
  InstanceVector::iterator begin;
  InstanceVector::iterator end;
  if (solState_.twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM)
  {
    begin = pdeInstancePtrVec_.begin ();
    end = pdeInstancePtrVec_.end ();
  }
  else
  {
    begin = instancePtrVec_.begin ();
    end = instancePtrVec_.end ();
  }

  for (iter=begin; iter!=end;++iter)
  {
    bool tmpBool =  (*iter)->isConverged();
    allDevsConv = allDevsConv && tmpBool;
  }

  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_LAND, &allDevsConv, 1);

  if (DEBUG_DEVICE && devOptions_.debugLevel > 0)
  {
    if (allDevsConv)
    {
      dout() << "All devices converged!" << std::endl;
    }
    else
    {
      dout() << "At least one device NOT converged!" << std::endl;
    }
  }

  return allDevsConv != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/22/06
//-----------------------------------------------------------------------------
bool DeviceMgr::innerDevsConverged()
{
  int innerDevsConv = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool tmpFlag = extern_device.isInnerSolveConverged();
      innerDevsConv = innerDevsConv && tmpFlag;
    }
  }

  Parallel::AllReduce(pdsMgrPtr_->getPDSComm()->comm(), MPI_LAND, &innerDevsConv, 1);

  return innerDevsConv != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupExternalDevices
// Purpose       : In parallel, we need to setup all external devices
//                 and appropriately setup the list of instances
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool DeviceMgr::setupExternalDevices()
{
  N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();

#ifdef Xyce_PARALLEL_MPI
  InstanceVector &extern_device_vector = modelTypeInstanceVector_[ExternDevice::Traits::modelType()];

  int procID = pdsCommPtr->procID();
  int numProc = pdsCommPtr->numProc();

  int numExt = extern_device_vector.size();
  int numExtTotal = 0;
  pdsCommPtr->sumAll(&numExt, &numExtTotal, 1);

  InstanceVector orig_extern_device_vector = extern_device_vector;

  //resetup the instance vector to have a size totally the global
  //number of ext devices
  if (numExtTotal > 0)
  {
    extern_device_vector.resize(numExtTotal);

    int loc = 0;
    for (int proc = 0; proc < numProc; ++proc)
    {
      int cnt = 0;
      if (proc == procID) cnt = numExt;
      pdsCommPtr->bcast(&cnt, 1, proc);

      for(int i = 0; i < cnt; ++i)
      {
        if (proc == procID)
        {
          ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*extern_device_vector[loc]);
          int size = extern_device.getInstanceBlock().packedByteCount();
          int bufSize = size+100;
          char *buf=new char[bufSize];
          pdsCommPtr->bcast(&size, 1, proc);
          int pos = 0;
          extern_device.getInstanceBlock().pack(buf, bufSize, pos, pdsCommPtr);
          pdsCommPtr->bcast(buf, size, proc);
          extern_device_vector[loc] = orig_extern_device_vector[i];
          delete [] buf;
        }
        else
        {
          int size = 0;
          pdsCommPtr->bcast(&size, 1, proc);
          int bufSize = size+100;
          char *buf=new char[bufSize];
          pdsCommPtr->bcast(buf, size, proc);
          int pos = 0;
          InstanceBlock instance_block;
          instance_block.unpack(buf, bufSize, pos, pdsCommPtr);
          extern_device_vector[loc] = addExtDeviceInstance_(instance_block);
          delete [] buf;
        }
        static_cast<ExternDevice::Instance *>(extern_device_vector[loc])->setOwningProc(proc);
        static_cast<ExternDevice::Instance *>(extern_device_vector[loc])->setComm(pdsCommPtr);
        ++loc;
      }
    }

    assert(loc == numExtTotal);
  }

#else // not Xyce_PARALLEL_MPI
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.setComm(pdsCommPtr);
    }
  }
#endif // Xyce_PARALLEL_MPI

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateExternalDevices_
// Purpose       : Do the actual solve of the external devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DeviceMgr::updateExternalDevices_()
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.runExternalDevice();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addExtDeviceInstance_
// Purpose       : adds an external device instance on the processors
//               : that don't actually own it.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/06
//-----------------------------------------------------------------------------
ExternDevice::Instance *
DeviceMgr::addExtDeviceInstance_(const InstanceBlock &instance_block)
{
  ExternDevice::Instance *external_instance = 0;

  ModelTypeId model_type;

  if (instance_block.getModelName().empty())
  {
    model_type = getModelGroup(instance_block.getInstanceName().getDeviceType());
  }
  else
  {
    model_type = modelTypeMap_[instance_block.getModelName()];
  }

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }
  }

  // Add an instance of this type.
  Device &device = getDeviceByModelType(model_type);
  DeviceInstance *instance = device.addInstance(instance_block, FactoryBlock(devOptions_, solState_, matrixLoadData_, externData_, commandLine_));

  external_instance = static_cast<ExternDevice::Instance *>(instance);

  return external_instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::homotopyStepSuccess
// Purpose       :
// Special Notes : Needed for 2-level
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceMgr::homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals)
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.homotopyStepSuccess (paramNames, paramVals);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::homotopyStepFailure
// Purpose       :
// Special Notes : Needed for 2-level
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void DeviceMgr::homotopyStepFailure ()
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.homotopyStepFailure ();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceMgr::stepSuccess(Analysis::CurrentMode analysis)
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.stepSuccess(analysis);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceMgr::stepFailure(Analysis::CurrentMode analysis)
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.stepFailure(analysis);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::acceptStep
// Purpose       : Communicate to devices that the current step is accepted
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void DeviceMgr::acceptStep()
{
  // The time history for the LTRA device(s) has to be tracked
  // separately because it can be truncated if the user specifies that
  // option. This has to be called before
  // LTRAInstance::acceptStep(). Note that the DCOP is stored at
  // index zero and the first time step is stored at index 1 for the
  // time history vectors associated with the LTRA
  if (solState_.dcopFlag)
  {
    solState_.ltraTimeIndex = 0;
    solState_.ltraTimeHistorySize = 10;
    solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
  }
  else
  {
    solState_.ltraTimeIndex++;
    if (solState_.ltraTimeIndex >= solState_.ltraTimeHistorySize)
    {
      solState_.ltraTimeHistorySize += 10;
      solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
    }
    solState_.ltraTimePoints[solState_.ltraTimeIndex] = solState_.currTime;
  }

  // make sure we have the right solver state.  Putting this here
  // makes it essential that acceptStep be called BEFORE any changes have
  // been made to the solution or state vectors (e.g. rotating them) and
  // while "nextTime" in the solver still references the time for which the
  // solution vector is valid.

  bool tmpBool = setupSolverInfo_();
  solState_.acceptedTime = solState_.currTime;

  for (InstanceVector::iterator iter = instancePtrVec_.begin(); iter != instancePtrVec_.end(); ++iter)
  {
    (*iter)->acceptStep();
  }

  // If the TRYCOMPACT option is set then the LTRA model will try to
  // compact the amount of data stored and speed up the convolutions. If
  // any of the LTRA instances request this, done in their acceptStep()
  // member function, then perform the time-step compaction here.
  if (solState_.ltraDoCompact)
  {
    solState_.ltraTimePoints[solState_.ltraTimeIndex-1] =
      solState_.ltraTimePoints[solState_.ltraTimeIndex];

    solState_.ltraTimeIndex--;

    // reset the flag for the next time step
    solState_.ltraDoCompact = false;
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/07
//-----------------------------------------------------------------------------
bool DeviceMgr::getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec)
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    int numExt = (*model_type_it).second.size();

    tleVec.resize(numExt);
    int iext = 0;
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it, ++iext)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.getInitialQnorm(tleVec[iext]);
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInnerLoopErrorSums
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::getInnerLoopErrorSums (
  std::vector<N_TIA_TwoLevelError> & tleVec)
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    int numExt = (*model_type_it).second.size();

    tleVec.resize(numExt);
    int iext = 0;
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it, ++iext)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.getInnerLoopErrorSum(tleVec[iext]);
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateStateArrays()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::updateStateArrays()
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.updateStateArrays();
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool DeviceMgr::startTimeStep ()
{
  bool bsuccess = true;
  bool tmpBool = setupSolverInfo_();

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end()) {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.startTimeStep ();
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setExternalSolverState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceMgr::setExternalSolverState (const SolverState & ss)
{
  externalStateFlag_ = true;
  solStateExternal_ = ss;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
int DeviceMgr::restartDataSize(bool pack)
{
  int numdoubles = solState_.ltraTimePoints.size();
  int numSize_t = 3;
  int count = sizeof(double) * (numdoubles);
  count += sizeof(size_t) * numSize_t;

  // bump up size for unpacked data.  This is an empirical multiplier.
  if (!pack)
  {
    count *= 3;
  }

  return count;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dumpRestartData
// Purpose       : Output restart data.
// Special Notes : This function is called by the restart manager to output
//                 persistent data for the device package.  It should NOT
//                 include any data from individual devices, as that restart
//                 data is collected elsewhere.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::dumpRestartData
(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  bool retval=true;

  if (pack)
  {
    size_t size=solState_.ltraTimePoints.size();
    comm->pack(&(solState_.ltraTimeIndex), 1, buf, bsize, pos);
    comm->pack(&(solState_.ltraTimeHistorySize), 1, buf, bsize, pos);
    comm->pack(&(size), 1, buf, bsize, pos);
    comm->pack(&(solState_.ltraTimePoints[0]), size, buf, bsize, pos);
  }
  else
  {
    int count = restartDataSize(false);
    int startIndex = pos;
    for(int i = startIndex; i < (startIndex+count); ++i) buf[i] = ' ';

    size_t size=solState_.ltraTimePoints.size();
    std::ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(std::ios::scientific);
    ost << solState_.ltraTimeIndex << " ";
    ost << solState_.ltraTimeHistorySize << " ";
#ifdef Xyce_DEBUG_RESTART
    dout() <<
      "DeviceMgr::getRestartData:  ltraTimeIndex = " << solState_.ltraTimeIndex <<std::endl;
    dout() <<
      "DeviceMgr::getRestartData:  ltraTimeHistorySize = " << solState_.ltraTimeHistorySize <<std::endl;
#endif
    ost << size << " ";
    for (int i=0;i<size;i++)
    {
      ost << solState_.ltraTimePoints[i] << " ";
#ifdef Xyce_DEBUG_RESTART
    dout() <<
      "DeviceMgr::dumpRestartData:  ltraTimePoints["<<i<<"] ="
      << solState_.ltraTimePoints[i]<<std::endl;
#endif
    }

    std::string data(ost.str());
    for(unsigned int i = 0; i < data.length(); ++i) buf[startIndex+i] = data[i];

    // The line above copies the characters of the data string into buf,
    // but doesn't null-terminate buf.
    // it is essential to terminate the buffer with a null, or attempts
    // to construct a string object from it will get memory access problems.
    buf[startIndex+data.length()] = '\0';
    pos += data.length();
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restoreRestartData
// Purpose       : Load restart data.
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  bool retval=true;

  if (pack)
  {
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeIndex), 1);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeHistorySize), 1);
    size_t size=0;
    comm->unpack(buf, bsize, pos, &(size), 1);
    solState_.ltraTimePoints.resize(size);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimePoints[0]), size);
  }
  else
  {
    std::string str1(buf);
    int length = str1.size() - pos;
    std::string str2(str1,pos,length);

    std::istringstream ist(str2);

    ist >> solState_.ltraTimeIndex;
    ist >> solState_.ltraTimeHistorySize;
#ifdef Xyce_DEBUG_RESTART
    dout() <<
      "DeviceMgr::restoreRestartData:  ltraTimeIndex = " << solState_.ltraTimeIndex <<std::endl;
    dout() <<
      "DeviceMgr::restoreRestartData:  ltraTimeHistorySize = " << solState_.ltraTimeHistorySize <<std::endl;
#endif
    size_t size=0;
    ist >> size;
    solState_.ltraTimePoints.resize(size);
    for (int i=0;i<size;++i)
    {
      ist >> solState_.ltraTimePoints[i];
#ifdef Xyce_DEBUG_RESTART
    dout() << "DeviceMgr::restoreRestartData:  ltraTimePoints["<<i<<"] = " << solState_.ltraTimePoints[i] << std::endl;
#endif
    }

    pos += ist.tellg();
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceEntity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
DeviceEntity * DeviceMgr::getDeviceEntity(const std::string & full_param_name) const
{
  std::string entity_name = Xyce::Util::entityNameFromFullParamName(full_param_name).getEncodedName();

  DeviceEntityMap::iterator it = parameterDeviceCache_.find(full_param_name);
  if (it == parameterDeviceCache_.end())
  {
    DeviceEntity *device_entity = findDeviceEntity(deviceMap_.begin(), deviceMap_.end(), entity_name);
    parameterDeviceCache_[full_param_name] = device_entity;
    return device_entity;
  }
  else if (!(*it).second)
  {
    DeviceEntity *device_entity = findDeviceEntity(deviceMap_.begin(), deviceMap_.end(), entity_name);
    (*it).second = device_entity;
    return device_entity;
  }
  else
    return (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDevicesToCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/06/14
//-----------------------------------------------------------------------------
void DeviceMgr::addDevicesToCount(const std::map<std::string,int> & device_map)
{
  // Add in devices to map for counting.
  std::map<std::string,int>::const_iterator dm_begin = device_map.begin();
  std::map<std::string,int>::const_iterator dm_end = device_map.end();
  std::map<std::string,int>::const_iterator dm_iter = dm_begin;
  for ( ; dm_iter != dm_end; ++dm_iter)
  {
    if ( localDeviceCountMap_[dm_iter->first] )
    {
       localDeviceCountMap_[dm_iter->first] += dm_iter->second;
    }
    else
    {
       localDeviceCountMap_[dm_iter->first] = dm_iter->second;
    }
  }
}

} // namespace Device
} // namespace Xyce
