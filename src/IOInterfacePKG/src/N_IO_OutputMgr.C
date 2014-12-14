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
// Filename       : $RCSfile: N_IO_OutputMgr.C,v $
//
// Purpose        : Output Manager
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
// Revision Number: $Revision: 1.457.2.6 $
//
// Revision Date  : $Date: 2014/08/29 21:26:34 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Misc.h>

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_Objective.h>
#include <N_IO_Op.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_OutputFileBase.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterHomotopy.h>
#include <N_IO_OutputterDC.h>
#include <N_IO_OutputterAC.h>
#include <N_IO_OutputterHB.h>
#include <N_IO_OutputterMPDE.h>
#include <N_IO_OutputterTransient.h>
#include <N_IO_OutputterSensitivity.h>
#include <N_IO_OutputterRawOverride.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_mmio.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {

namespace SensitivityOptions {

enum {
  DIRECT   = 0x01,
  ADJOINT  = 0x02,
  SCALED   = 0x04,
  UNSCALED = 0x08
};

}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
NodeNamePairMap::const_iterator findNode(
    const std::string &name,
    const NodeNamePairMap &node_map,
    const AliasNodeMap &alias_map)
{
  NodeNamePairMap::const_iterator node_it = node_map.find(name);
  if (node_it == node_map.end()) {
    AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
    if (alias_node_it != alias_map.end())
      node_it = node_map.find((*alias_node_it).second);
  }

  return node_it;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::OutputMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::OutputMgr(
  CmdParse &                    command_line,
  Util::Op::BuilderManager &    op_builder_manager)
  : title_(command_line.getArgumentValue("netlist")),
    netListFilename_(command_line.getArgumentValue("netlist")),
    filenameSuffix_(),
    opBuilderManager_(op_builder_manager),
    topology_(0),
    deviceInterface_(0),
    analysisManager_(0),
    enableHomotopyFlag_(false),
    enableSensitivityFlag_(false),
    sensitivityOptions_(0),
    PRINTdcstart_(0.0),
    PRINTdcstop_(0.0),
    PRINTdcvalue_(0.0),
    PRINTdcname_(""),
    loadFlag_(false),
    outputOnceAlreadyFlag_(false),
    initialOutputInterval_(0.0),
    STEPEnabledFlag_(false),
    printEndOfSimulationLine_(true),
    outputVersionInRawFile_(false),
    tempSweepFlag_(false),
    outputCalledBefore_(false),
    dcLoopNumber_(0),
    maxDCSteps_(0),
    hdf5FileNameGiven_(false),
    hdf5HeaderWritten_(false),
    hdf5IndexValue_(0)
{
  if (command_line.getArgumentValue("-delim") == "TAB")
    defaultPrintParameters_.delimiter_ = "\t";
  else if (command_line.getArgumentValue("-delim") == "COMMA")
    defaultPrintParameters_.delimiter_ = ",";
  else
    defaultPrintParameters_.delimiter_ = command_line.getArgumentValue("-delim");

  if (command_line.argExists("-a"))
    defaultPrintParameters_.asciiRaw_ = true;

  if (command_line.argExists("-r")) {
    defaultPrintParameters_.overrideRaw_ = true;
    defaultPrintParameters_.filename_ = command_line.getArgumentValue("-r");
    defaultPrintParameters_.format_ = defaultPrintParameters_.asciiRaw_ ? Format::RAW_ASCII : Format::RAW;
  }

  // If the output file is specified on the command line it takes precedence
  else if (command_line.argExists("-o"))
    defaultPrintParameters_.filename_ = command_line.getArgumentValue("-o");
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::~OutputMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::~OutputMgr()
{
  for (OutputterMap::iterator it = outputterMap_.begin(); it != outputterMap_.end(); ++it)
    for (OutputterMap::mapped_type::iterator it2 = (*it).second.begin(); it2 != (*it).second.end(); ++it2)
      delete (*it2);

  for (OpenPathStreamMap::iterator it = openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it)
    delete (*it).second.second;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerAnalysisManager
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool OutputMgr::registerAnalysisManager(
  Analysis::AnalysisManager *   analysis_manager)
{
  analysisManager_ = analysis_manager;

  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *this);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::notify
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgr::notify(
  const Analysis::StepEvent &   step_event)
{
  switch (step_event.state_) {
    case Analysis::StepEvent::INITIALIZE:
      outputState_.stepMaxCount_ = step_event.count_;
      break;

    case Analysis::StepEvent::STEP_STARTED:
      outputState_.stepLoopNumber_ = step_event.count_;
      startStep(outputState_.stepLoopNumber_, outputState_.stepMaxCount_);
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      break;

    case Analysis::StepEvent::FINISH:
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file in given mode, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *
OutputMgr::openFile(
  const std::string &           path,
  std::ios_base::openmode       mode)
{
  OpenPathStreamMap::iterator it= openPathStreamMap_.find(path);

  if (path == "CONSOLE")
    return &Xyce::dout();
  else if (it != openPathStreamMap_.end()) {
    ++(*it).second.first;
    return (*it).second.second;
  }
  else {
    std::ostream *os = new std::ofstream(path.c_str(), mode);
    openPathStreamMap_[path] = std::pair<int, std::ostream *>(1, os);

    if (!os->good())
    {
      Report::UserFatal0() << "Failure opening " << path;
    }

    return os;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file for output only, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *
OutputMgr::openFile(
  const std::string &           path)
{
  return openFile(path, std::ios_base::out);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openBinaryFile
// Purpose       : open named file in binary mode for output only, create
//                 stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *
OutputMgr::openBinaryFile(
  const std::string &           path)
{
  return openFile(path, std::ios_base::out | std::ios_base::binary);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::closeFile
// Purpose       : Close given stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
int
OutputMgr::closeFile(
  std::ostream *                os)
{
  if (os == &Xyce::dout())
    return 1;

  int open_count= 0;

  for (OpenPathStreamMap::iterator it= openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it) {
    if ((*it).second.second == os) {
      open_count= --(*it).second.first;
      if (open_count == 0) {
        delete os;
        openPathStreamMap_.erase(it);
        break;
      }
    }
  }

  return open_count;
}


namespace {

//-----------------------------------------------------------------------------
// Function      : testAndSet
// Purpose       : test if set contains element and then add the element
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jun 24 09:38:58 2014
//-----------------------------------------------------------------------------
///
/// Returns true of set contained element on entry, adds element if it was not found
///
/// @invariant set contains new element on completion
///
/// @param s Set to find element
/// @param t Element to find
///
/// @return true if set contained element before insertion
///
template<class S, class T>
inline bool testAndSet(S &s, const T &t)
{
  bool found = s.find(t) != s.end();
  if (!found)
    s.insert(t);
  return found;
}

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::prepareOutput
// Purpose       : checks to make sure requested values can be output
// Special Notes : Primary purpose is for error checking BEFORE calculation.
//                 Also allows any name translation for output, if needed.
//
//                 This function also now serves as an initialization for
//                 allNodes_.  allNodes_ needs to be set up after topology
//                 is done setting up.
//
// erkeite  8/10/2007:
// Note that this is currently(august 2007) one of 3 different functions
// that check the .print line for errors.  That is a bit confusing,
// unfortunately.  For commentary and some explanation of this function,
// as well as the other two, see the comments for the delayedPrintLineDiagnostics
// function.  This whole functionality really should be refactored, or
// at least made more transparent to other developers.
//
// Function renamed from check_output to prepareOutput and heavily reorganized
// by David Baur on 6/28/2013.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
void OutputMgr::prepareOutput(
  Parallel::Machine                     comm,
  Analysis::Analysis_Mode               analysis_mode)
{
  // Setup rawfile if requested
  if (defaultPrintParameters_.overrideRaw_)
  {
    Outputter::enableRawOverrideOutput(comm, *this, analysis_mode);

    if (activeOutputterStack_.empty())
      activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());

    addActiveOutputter(PrintType::RAW_OVERRIDE, analysis_mode);
  }
  else
  {
    switch (analysis_mode)
    {
      case Analysis::ANP_MODE_INVALID:
      case Analysis::ANP_MODE_DC_OP:
      case Analysis::ANP_MODE_MOR:
        break;

      case Analysis::ANP_MODE_DC_SWEEP:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_DC_SWEEP))
          Outputter::enableDCOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::TRAN, analysis_mode);
        break;

      case Analysis::ANP_MODE_TRANSIENT:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_TRANSIENT))
          Outputter::enableTransientOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::TRAN, analysis_mode);
        break;

      case Analysis::ANP_MODE_MPDE:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_MPDE))
          Outputter::enableMPDEOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::MPDE, analysis_mode);
        addActiveOutputter(PrintType::MPDE_IC, analysis_mode);
        break;

      case Analysis::ANP_MODE_HB:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_HB))
          Outputter::enableHBOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::HB, analysis_mode);
        addActiveOutputter(PrintType::HB_IC, analysis_mode);
        addActiveOutputter(PrintType::HB_STARTUP, analysis_mode);
        break;

      case Analysis::ANP_MODE_AC:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_AC))
          Outputter::enableACOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::AC, analysis_mode);
        addActiveOutputter(PrintType::AC_IC, analysis_mode);
        break;
    }

    if (enableHomotopyFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Analysis_Mode) (Analysis::ANP_MODE_INVALID + 100)))
        Outputter::enableHomotopyOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::HOMOTOPY, analysis_mode);
    }

    if (enableSensitivityFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Analysis_Mode) (Analysis::ANP_MODE_INVALID + 101)))
        Outputter::enableSensitivityOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::SENS, analysis_mode);
    }
  }

  if (outputterMap_.empty())
  {
    Outputter::TimePrn *outputter = new Outputter::TimePrn(comm, *this, defaultPrintParameters_);
    outputterMap_[PrintType::TRAN].push_back(outputter);
  }

  measureManager_->makeMeasureOps(comm, opBuilderManager_);
  fourierManager_->fixupFourierParameters(comm, opBuilderManager_);
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Class         : STEPOptionsReg
// Purpose       : functor for registering STEP options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct STEPOptionsReg : public PkgOptionsReg
{
  STEPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSTEPOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : DCOptionsReg
// Purpose       : functor for registering DC options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct DCOptionsReg : public PkgOptionsReg
{
  DCOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return outputManager_.registerDCOptions(options);
  }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : TranOptionsReg
// Purpose       : functor for registering transient options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct TranOptionsReg : public PkgOptionsReg
{
  TranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerTranOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : TranOptionsReg
// Purpose       : functor for registering transient options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct MPDETranOptionsReg : public PkgOptionsReg
{
  MPDETranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerMPDETranOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : HBOptionsReg
// Purpose       : functor for registering HB options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct HBOptionsReg : public PkgOptionsReg
{
  HBOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerHBOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OptionsReg
// Purpose       : functor for registering Output options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OptionsReg : public PkgOptionsReg
{
  OptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerOutputOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : DeviceOptionsReg
// Purpose       : functor for registering Device options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct DeviceOptionsReg : public PkgOptionsReg
{
  DeviceOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerDeviceOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : PrintOptionsReg
// Purpose       : functor for registering Print options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct PrintOptionsReg : public PkgOptionsReg
{
  PrintOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()(const Util::OptionBlock &print_block)
  {
    return outputManager_.parsePRINTBlock(print_block);
  }

  OutputMgr &   outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : ObjectiveOptionsReg
// Purpose       : functor for registering Objective options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct ObjectiveOptionsReg : public PkgOptionsReg
{
  ObjectiveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.setOBJECTIVEParams( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : LoadOptionsReg
// Purpose       : functor for registering Load options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct LoadOptionsReg : public PkgOptionsReg
{
  LoadOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerLoad( options ); }

  OutputMgr &outputManager_;
};

//-----------------------------------------------------------------------------
// Class         : SensReg
// Purpose       : functor for registering sensitivity options
// Special Notes : Used by package manager submitRegistration method
// Creator       : Eric Keiter
// Creation Date : 02/10/2014
//-----------------------------------------------------------------------------
struct SensReg : public PkgOptionsReg
{
  SensReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return outputManager_.registerSens(options);
  }

  OutputMgr &   outputManager_;
};

//-----------------------------------------------------------------------------
// Class         : SensOptionsReg
// Purpose       : functor for registering sensitivity options
// Special Notes : Used by package manager submitRegistration method
// Creator       : Eric Keiter
// Creation Date : 02/10/2014
//-----------------------------------------------------------------------------
struct SensOptionsReg : public PkgOptionsReg
{
  SensOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSensOptions( options ); }

  OutputMgr &outputManager_;
};

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerPkgOptionsMgr(PkgOptionsMgr &pkgOpt)
{
  pkgOpt.submitRegistration(
    "DC", netListFilename_, new DCOptionsReg(*this));

  pkgOpt.submitRegistration(
    "TRAN", netListFilename_, new TranOptionsReg(*this));

  pkgOpt.submitRegistration(
    "MPDE", netListFilename_, new MPDETranOptionsReg(*this));

  pkgOpt.submitRegistration(
    "HB", netListFilename_, new HBOptionsReg(*this));

  pkgOpt.submitRegistration(
    "STEP", netListFilename_, new STEPOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OUTPUT", netListFilename_, new OptionsReg(*this));

  pkgOpt.submitRegistration(
    "PRINT", netListFilename_, new PrintOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OBJECTIVE", netListFilename_, new ObjectiveOptionsReg(*this));

  pkgOpt.submitRegistration(
    "LOAD", netListFilename_, new LoadOptionsReg(*this));

  pkgOpt.submitRegistration(
    "DEVICE", netListFilename_, new DeviceOptionsReg(*this));

  pkgOpt.submitRegistration(
    "SENS", netListFilename_, new SensReg(*this));

  pkgOpt.submitRegistration(
    "SENSITIVITY", netListFilename_, new SensOptionsReg(*this));

  registerPkgOptionsMgrInitialConditions(pkgOpt);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDCOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerDCOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator
       iterPL= option_block.getParams().begin();
       iterPL != option_block.getParams().end(); ++iterPL)
  {
    if (iterPL->tag() == "PARAM")
    {
      dcParams_.push_back(iterPL->stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerTranOptions(const Util::OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerMPDETranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerMPDETranOptions(const Util::OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerHBOptions(const Util::OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSTEPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerSTEPOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator
      iterPL= option_block.getParams().begin();
      iterPL != option_block.getParams().end();
      ++iterPL)
  {
    if (iterPL->tag() == "PARAM")
    {
      stepParams_.push_back(iterPL->stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDeviceOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/10/10
//-----------------------------------------------------------------------------
bool OutputMgr::registerDeviceOptions(const Util::OptionBlock & option_block)
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerOutputOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/30/01
//-----------------------------------------------------------------------------
bool OutputMgr::registerOutputOptions(const Util::OptionBlock & OB)
{

  // Need to capture a variable length of output intervals.  The use case is:
  //
  // .OPTIONS OUTPUT INITIAL_INTERVAL= <interval> [<t0> <i0> [<t1> <i1>...]]
  //
  // additional options can appear before or after the INITIAL_INTERVAL as in
  //
  //.OPTIONS OUTPUT HDF5FILE= xxxx INITIAL_INTERVAL=<interval> [<t0> <i0> [<t1> <i1>...]]
  // or
  // .OPTIONS OUTPUT INITIAL_INTERVAL= <interval> [<t0> <i0> [<t1> <i1>...]] HDF5FILE=xxxx
  //

  ParameterList::const_iterator iterPL= OB.getParams().begin();
  while (iterPL != OB.getParams().end())
  {
    if (iterPL->tag() == "INITIAL_INTERVAL")
    {
      // start handling a list of intervals
      // set interval value
      initialOutputInterval_= iterPL->getImmutableValue<double>();
      // look for optional time pairs.
      outputIntervalPairs_.clear();
      bool doneWithTimePairs= false;

      // need to point at next parameter
      ++iterPL;

      while ((iterPL != OB.getParams().end()) && !doneWithTimePairs)
      {
        if (iterPL->tag() == "TIME")
        {
          double t= iterPL->getImmutableValue<double>();
          ++iterPL;
          double iv= iterPL->getImmutableValue<double>();
          ++iterPL;

          outputIntervalPairs_.push_back(std::pair<double, double>(t, iv));
        }
        else
        {
          // didn't find a time pair so bail out of this loop.
          doneWithTimePairs= true;
        }
      }
    }
    else if (iterPL->tag()=="HDF5FILENAME")
    {
    // look for other option tags
      hdf5FileNameGiven_= true;
      hdf5FileName_= iterPL->stringValue();
      ++iterPL;
    }
    else if (iterPL->tag()=="PRINTENDOFSIMLINE")
    {
      // look for flag to turn off "End of Xyce(TM) Simulation" line
      printEndOfSimulationLine_= iterPL->getImmutableValue<bool>();
      ++iterPL;
    }
    else if (iterPL->tag()=="OUTPUTVERSIONINRAWFILE")
    {
      // look for flag to toggle output of version in header of RAW file
     outputVersionInRawFile_ = iterPL->getImmutableValue<bool>();
     ++iterPL;
    }
    else
    {
      // silently ignore?
      ++iterPL;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setExternalNetlistParams
// Purpose       : Called from owning class to set any external parameters
//                 set by the user.  Note, these parameter lists may also
//                 request response functions that the output manger must report.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 08/07/2012
//-----------------------------------------------------------------------------
void OutputMgr::setExternalNetlistParams(
  const StringPairVector &      externalNetlistParams)
{
  // externalParams can contain section names from dakota like
  // variables 2, responses 4, derivatives 4.  If we don't find any sections tags
  // then we can assume that all the parameters are just variable names to set.
  // if we find tags, then use the "responses" section to record what response functions
  // need to be reported.

  StringPairVector *section = 0;
  for (StringPairVector::const_iterator it = externalNetlistParams.begin(), end = externalNetlistParams.end(); it != end; it++)
  {
    if ((*it).first == "variables")
      section = &variablesUsedInSimulation_;
    else if ((*it).first == "functions")
      section = &responseFunctionsRequested_;
    else if ((*it).first == "derivative_variables")
      section = &derivativeVariablesRequested_;
    else if ((*it).first == "analysis_components")
      section = &analysisComponentsRequested_;
    else if (section)
      section->push_back(*it);
  }
}

ParameterList
OutputMgr::getVariableList() const
{
  ParameterList parameter_list;

  for (OutputParameterMap::const_iterator it1 = outputParameterMap_.begin(), end1 = outputParameterMap_.end(); it1 != end1; ++it1)
  {
    const OutputParameterMap::mapped_type &parameter_vector = (*it1).second;

    for (OutputParameterMap::mapped_type::const_iterator it2 = parameter_vector.begin(), end2 = parameter_vector.end(); it2 != end2; ++it2)
    {
      const PrintParameters &print_parameters = (*it2);

      // Populate stringStat, nodeStat and instanceState with variable names from print line
      std::map<std::string, bool> stringStat;
      std::map<std::string, bool> nodeStat;
      std::map<std::string, bool> instanceStat;

      for (ParameterList::const_iterator it3 = print_parameters.variableList_.begin(), end3 = print_parameters.variableList_.end(); it3 != end3; ++it3)
      {
        const Util::Param &parameter = (*it3);
        parameter_list.push_back(parameter);
      }
    }
  }

  return parameter_list;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getVariableNames
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey
// Creation Date : 07/28/09
//-----------------------------------------------------------------------------
std::vector<std::string>
OutputMgr::getVariableNames()
{
  std::vector<std::string> name_list;

  for (OutputParameterMap::const_iterator it1 = outputParameterMap_.begin(), end1 = outputParameterMap_.end(); it1 != end1; ++it1)
  {
    const OutputParameterMap::mapped_type &parameter_vector = (*it1).second;

    for (OutputParameterMap::mapped_type::const_iterator it2 = parameter_vector.begin(), end2 = parameter_vector.end(); it2 != end2; ++it2)
    {
      const PrintParameters &print_parameters = (*it2);

      // Populate stringStat, nodeStat and instanceState with variable names from print line
      std::map<std::string, bool> stringStat;
      std::map<std::string, bool> nodeStat;
      std::map<std::string, bool> instanceStat;

      for (ParameterList::const_iterator it3 = print_parameters.variableList_.begin(), end3 = print_parameters.variableList_.end(); it3 != end3; ++it3)
      {
        const Util::Param &parameter = (*it3);
        name_list.push_back(parameter.tag());
      }
    }
  }

  return name_list;
}

namespace {

std::set<std::string>::const_iterator findNode(const std::string &name, const std::set<std::string> &node_set, const AliasNodeMap &alias_map)
{
  std::set<std::string>::const_iterator node_it = node_set.find(name);
  if (node_it == node_set.end()) {
    AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
    if (alias_node_it != alias_map.end())
      node_it = node_set.find((*alias_node_it).second);
  }

  return node_it;
}

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::printLineDiagnostics
// Purpose       :
// Special Notes : erkeite:  8/10/2007
//
// For commentary and some explanation of this function, see the comments
// for the delayedPrintLineDiagnostics function.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/18/06
//-----------------------------------------------------------------------------
void
printLineDiagnostics(
  Parallel::Machine                                                     comm,
  const OutputMgr &                                                     output_manager,
  const Measure::Manager &                                              measure_manager,
  const std::set<std::string> &                                         node_names,
  const std::map<std::string, RefCountPtr<Device::InstanceBlock> > &    device_map,
  IO::AliasNodeMap &                                                    alias_node_map)
{
  for (OutputParameterMap::const_iterator it1 = output_manager.outputParameterMap_.begin(), end1 = output_manager.outputParameterMap_.end(); it1 != end1; ++it1)
  {
    const OutputParameterMap::mapped_type &parameter_vector = (*it1).second;

    for (OutputParameterMap::mapped_type::const_iterator it2 = parameter_vector.begin(), end2 = parameter_vector.end(); it2 != end2; ++it2)
    {
      const PrintParameters &print_parameters = (*it2);

      // Populate stringStat, nodeStat and instanceState with variable names from print line
      std::map<std::string, bool> stringStat;
      std::map<std::string, bool> nodeStat;
      std::map<std::string, bool> instanceStat;

      for (ParameterList::const_iterator it3 = print_parameters.variableList_.begin(), end3 = print_parameters.variableList_.end(); it3 != end3; ++it3)
      {
        const Util::Param &parameter = (*it3);

        std::vector<std::string> nodes;
        std::vector<std::string> instances;
        std::vector<std::string> leads;
        std::vector<std::string> strings;
        std::vector<std::string> special;
        if (Util::hasExpressionTag(parameter) )
        {
          // parameter starts with "{" but may not have been been parsed into an expression.
          // check if there is an underlying expression object with the parameter
          if (parameter.getType() == Util::EXPR)
          {
            parameter.getValue<Util::Expression>().get_names(XEXP_NODE, nodes);
            parameter.getValue<Util::Expression>().get_names(XEXP_INSTANCE, instances);
            parameter.getValue<Util::Expression>().get_names(XEXP_LEAD, leads);
            parameter.getValue<Util::Expression>().get_names(XEXP_STRING, strings);
            parameter.getValue<Util::Expression>().get_names(XEXP_SPECIAL, special);  // special returns vars like TIME
          }
          else if ( ((parameter.getType()) == Util::DBLE) || ((parameter.getType()) == Util::INT) )
          {
          }
          instances.insert(instances.end(), leads.begin(), leads.end());
          strings.insert(strings.end(), special.begin(), special.end());  // add specials to strings
        }
        else
        {
          std::string varType= parameter.tag();
          if ((varType == "I" || ((varType.size() == 2 || varType.size() == 3) && varType[0] == 'I')) && parameter.getImmutableValue<int>() > 0)
          {
            // any devices found in this I(xxx) structure need to be communicated to the device manager
            // so that the lead currents can be calculated
            if (parameter.getImmutableValue<int>() != 1)
            {
              Report::UserError() << "Only one device argument allowed in I() in .PRINT command";
            }
            else
            {
              ++it3;
              if (varType.size() == 2)
              {
                instances.push_back((*it3).tag() + "{" + varType[1] + "}");
              }
              else
              {
                instances.push_back((*it3).tag());
              }
            }
          }
          else if ((varType == "V" || ((varType.size() == 2 || varType.size() == 3) && varType[0] == 'V')) && parameter.getImmutableValue<int>() > 0)
          {
            int numIndices = parameter.getImmutableValue<int>();

          if (numIndices < 1 || numIndices > 2)
          {
            Report::UserError0() << "Only one or two node arguments allowed in V() in .PRINT command";
          }
          else
          {
            for (; numIndices > 0 ; --numIndices)
            {
              ++it3;
              nodes.push_back((*it3).tag());
            }
          }
          }
        else if (varType == "N" && parameter.getImmutableValue<int>() > 0)
        {
          // Don't attempt to check this type of variable.
          //
          // N-variables can only be fully resolved once the devices are allocated,
          // and variables which are internal to the devices have been assigned names.
          //(N-variables are the same as V-variables, except that they can include
          // internal vars).
          if (parameter.getImmutableValue<int>() != 1)
            Report::UserError() << "Only one device argument allowed in N() in .PRINT command";
          else
            ++it3;
        }
        else
        {
          strings.push_back(parameter.tag());
        }
        }

      for (std::vector<std::string>::const_iterator iter_s= strings.begin(); iter_s != strings.end(); ++iter_s)
      {
        const std::string &name= *iter_s;

        bool done = false;
        if (name == "TEMP" || name == "TIME" || name == "FREQUENCY" || name == "INDEX" || name == "OBJFUNC" || name == "SENS" || name == "sweep")
          done = true;

        if (!done)
        {
          // check if this is a step sweep value.
          for (std::vector<std::string>::const_iterator it= output_manager.dcParams_.begin(), end = output_manager.dcParams_.end(); it != end ; ++it)
          {
            if (*it == name)
            {
              done = true;
              break;
            }
          }
        }

        if (!done)
        {
          // check if this is a step sweep value.
          for (std::vector<std::string>::const_iterator it= output_manager.stepParams_.begin(), end = output_manager.stepParams_.end(); it != end ; ++it)
          {
            if (*it == name)
            {
              done = true;
              break;
            }
          }
        }

        double mValue;
        bool mFound;
        measure_manager.getMeasureValue(name, mValue, mFound);
        if (mFound)
          done = true;

        if (!done)
        {
          done = output_manager.deviceInterface_->findParam(name);

          // Note:  when this is called, the devices haven't been allocated
          // yet.  So, this particular check isn't reliable.  If the device and/or
          // parameter is not found by the device package via getParam, then
          // simply save it for a later diagnostic.
          if (!done)
          {
            output_manager.deferredParameterCheck_.insert(name);
          }
          done = true;
        }
        stringStat[name]= done;
      }

      for (std::vector<std::string>::const_iterator iter_s = nodes.begin(); iter_s != nodes.end(); ++iter_s)
      {
        const std::string &name= *iter_s;

        // allow * to pass
        nodeStat[name] = (name == "*") || (findNode(name, node_names, alias_node_map) != node_names.end());
      }

      for (std::vector<std::string>::iterator iter_s= instances.begin(); iter_s != instances.end(); ++iter_s)
      {
        const std::string &name= *iter_s;

        if (name.substr(name.size() - 1, 1) == "}")
        {
          // allow * to pass
          instanceStat[name]= (name == "*") || (device_map.find(name.substr(0, name.size() - 3)) != device_map.end());
        }
        else
        {
          // allow * to pass
          instanceStat[name]= (name == "*") || (device_map.find(name) != device_map.end());
        }
      }
      }

    // Sum uses of names in nodeStat and instanceStat
    std::vector<int> stat;
    for (std::map<std::string, bool>::const_iterator it = nodeStat.begin(), end = nodeStat.end(); it != end; ++it)
    {
      stat.push_back((*it).second ? 1 : 0);
    }

    for (std::map<std::string, bool>::const_iterator it = instanceStat.begin(), end = instanceStat.end(); it != end; ++it)
    {
      stat.push_back((*it).second ? 1 : 0);
    }

    Parallel::AllReduce(comm, MPI_SUM, stat);

// Generate message
    std::ostringstream oss;

    int count = 0;
    int stat_index = 0;
    for (std::map<std::string, bool>::const_iterator it= nodeStat.begin(), end = nodeStat.end(); it != end; ++it)
    {
      if (stat[stat_index++] == 0)
      {
        if (count != 0)
        {
          oss << ", ";
        }
        oss << "node " << (*it).first;
        ++count;
      }
    }

    for (std::map<std::string, bool>::const_iterator it= instanceStat.begin(), end = instanceStat.end(); it != end; ++it)
    {
      if (stat[stat_index++] == 0)
      {
        if (count != 0)
        {
          oss << ", ";
        }
        if (((*it).first)[((*it).first).size() - 1] == '}')
        {
          oss << "I" << ((*it).first).substr(((*it).first).size() - 2, 1)
              << "("
              << ((*it).first).substr(0, ((*it).first).size() - 3)
              << ")";
          ++count;
        }
        else
        {
          oss << "I(" << (*it).first << ")";
          ++count;
        }
      }
    }

    for (std::map<std::string, bool>::const_iterator it= stringStat.begin(); it != stringStat.end(); ++it)
    {
      if (!(*it).second)
      {
        if (count != 0)
        {
          oss << ", ";
        }
        oss << (*it).first;
        ++count;
      }
    }

    if (count != 0)
    {
      Report::UserError0().at(print_parameters.netlistLocation_)
        << "There " << (count == 1 ? "was " : "were ") << count << " undefined symbol"
        << (count == 1 ? "" : "s") << " in .PRINT command: " << oss.str();
    }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::delayedPrintLineDiagnostics
// Purpose       : Yet another function that checks the .print line.
//
// Special Notes : erkeite.  8/10/2007
//
//                 This function has to be called after the devices in
//                 the device package have all been allocated.  The function
//                 printLineDiagnostics is called earlier than that, so
//                 it doesn't handle device-specific parameter outputs.
//
//                 We have 3 completely different functions that do different
//                 aspects of checking the print line:
//
//                (1) printLineDiagnostics:  happens really early, as soon as
//                 topology is constructed, but before devivces are allocated.
//
//                (2) delayedPrintLineDiagnostics:  happens just a little later, to
//                 check .print line fields that are specific to device entities(ie
//                 aren't just I and V solution vars).
//
//                (3) check_output:  Happens at the beginning of the simulation,
//                 after everything is completely set up, but before any solves
//                 happen.  The advantage of this one is that it actually uses
//                 all the same function calls as the regular output call that
//                 will happen later.  So, it is probably the most thorough test
//                 of all.  However, for really big simulations, it happens relatively
//                 late in the setup.  For really large simulations, setup can take
//                 a long time, and it is often useful to get diagnostics earlier
//                 than that.
//
//                 Note from TVR 10/22/2013:  "check_output", referred to
//                 above, has been broken apart and renamed "prepareOutput".
//                 The notes from ERK above are therefore somewhat out of date.
//                 The refactor, though, is not complete, and the replacement
//                 code is not well documented.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
void deferredPrintLineDiagnostics(
  Parallel::Machine     comm,
  OutputMgr &           output_manager)
{
  std::ostringstream oss;
  int count = 0;
  for (std::set<std::string>::const_iterator it = output_manager.deferredParameterCheck_.begin(), end = output_manager.deferredParameterCheck_.end(); it != end; ++it)
  {
    double result = 0.0;
    if (!output_manager.deviceInterface_->getParamAndReduce(*it, result))
    {
      if (count != 0)
        oss << ", ";
      oss << *it;
      ++count;
    }
  }

  if (count != 0)
  {
    Report::UserError0() << "There " << (count == 1 ? "was " : "were ")
                         << count << " undefined symbol" << (count == 1 ? "" : "s")
                         << " in .PRINT command: " << oss.str();
  }
}

template <class T>
class null_back_insert_iterator: public std::back_insert_iterator<T>
{
public:
  null_back_insert_iterator(T& t)
  {}
  null_back_insert_iterator& operator=(const T &t) {
    return *this;
  }
  null_back_insert_iterator& operator*() {
    return *this;
  }
  null_back_insert_iterator& operator++() {
    return *this;
  }
  null_back_insert_iterator operator++ (int) {
    return *this;
  }
};

void OutputMgr::checkPrintParameters(Parallel::Machine comm)
{
  topology_->returnVarTypeVec(typeRefs_);
  topology_->getNodeNames(allNodes_);
  topology_->getStateNodeNames(stateNodes_);
  topology_->getStoreNodeNames(storeNodes_);
  topology_->getExternNodeNames(externNodes_);

  Util::Op::OpList tempOpList;

  for (OutputParameterMap::const_iterator it = outputParameterMap_.begin(), end = outputParameterMap_.end(); it != end; ++it)
  {
    for (std::vector<PrintParameters>::const_iterator it2 = (*it).second.begin(), end2 = (*it).second.end(); it2 != end2; ++it2)
    {
      PrintParameters print_parameters = (*it2);
      fixupPrintParameters(comm, print_parameters);
      makeOps(comm, opBuilderManager_, print_parameters.netlistLocation_,
              print_parameters.variableList_.begin(), print_parameters.variableList_.end(), std::back_inserter<Util::Op::OpList>(tempOpList));
    }
  }

  for (Util::Op::OpList::iterator it = tempOpList.begin(), end = tempOpList.end(); it != end; ++it)
    delete *it;

  if (hdf5FileNameGiven_)
  {
    prepareHDF5Output(comm);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::getOutputIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/30/01
//-----------------------------------------------------------------------------
bool
OutputMgr::getOutputIntervals(
    double &                                    initialInterval,
    std::vector< std::pair<double, double> > &  intervalPairs) const
{
  initialInterval = initialOutputInterval_;
  intervalPairs = outputIntervalPairs_;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_IO OutputMgr::setOBJECTIVEParams
// Purpose       : Sets the OBJECTIVE calculation parameters.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/07/05
//-----------------------------------------------------------------------------
bool OutputMgr::setOBJECTIVEParams(const Util::OptionBlock & option_block)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "In N_IO OutputMgr::setoption_OBJECTIVEParams" << std::endl;
#endif

  std::string name;
  for (ParameterList::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end(); it != end; ++it)
  {
    if ((*it).uTag() == "NAME")
    {
      name = (*it).usVal();
      if (objectiveMap_.find(name) != objectiveMap_.end())
      {
        Report::UserFatal0() << "Duplicate objective name " << name;
      }
      break;
    }
  }

  std::pair<ObjectiveMap::iterator, bool> result = objectiveMap_.insert(ObjectiveMap::value_type(name, Objective()));

  (*result.first).second.initialize(option_block);

  return true;
}

namespace {
struct IsFallback
{
  bool operator()(const PrintParameters &print_parameters)
  {
    return print_parameters.fallback_;
  }
};
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addOutputParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jul 23 07:21:53 2014
//-----------------------------------------------------------------------------
///
/// Adds the print parameters to the output type as follows:
///
/// If the new print parameters is a fallback and no print parameters
/// have been provided then the print parameters are added.
///
/// If the new print is not a fallback (is specific), then the fallback
/// print parameters are removed and this one is added to the list.
///
/// @invariant
///
/// @param output_type          output type to add the print parameters too
/// @param print_parameters     print parameters to add
///
///
void OutputMgr::addOutputPrintParameters(
  OutputType::OutputType        output_type,
  const PrintParameters &       print_parameters)
{
  OutputParameterMap::mapped_type &print_parameters_vector = outputParameterMap_[output_type];

  if (print_parameters.fallback_ && outputParameterMap_[output_type].empty())
    print_parameters_vector.push_back(print_parameters);
  else if (!print_parameters.fallback_) {
    print_parameters_vector.erase(std::remove_if(print_parameters_vector.begin(), print_parameters_vector.end(), IsFallback()), print_parameters_vector.end());
    print_parameters_vector.push_back(print_parameters);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::parsePRINTBlock
// Purpose       : registers set of variables to print for .PRINT
//
// Special Notes : ERK.  This function used to also output the header.
//                 Now, that functionality is down in outputPRINTHeader_.
//
//                 Also, this function allocates expression handlers for
//                 any expressions specified.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
bool OutputMgr::parsePRINTBlock(
  const Util::OptionBlock &     print_block)
{
  PrintParameters       print_parameters = defaultPrintParameters_;
  PrintType::PrintType  print_type;
  Format::Format        format = Format::STD;
  bool                  no_index = false;

  print_parameters.netlistLocation_ = print_block.getNetlistLocation();

  ParameterList::const_iterator iterParam = print_block.getParams().begin();
  for (; iterParam != print_block.getParams().end(); ++iterParam)
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
#endif

    if (iterParam->tag() == "WIDTH")
    {
      print_parameters.streamWidth_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "TYPE")
    {
      std::string s = iterParam->stringValue();
      if (s == "AC")
        print_type = PrintType::AC;
      else if (s == "AC_IC")
        print_type = PrintType::AC_IC;
      else if (s == "DC")
        print_type = PrintType::DC;
      else if (s == "HB")
        print_type = PrintType::HB;
      else if (s == "HB_TD")
        print_type = PrintType::HB_TD;
      else if (s == "HB_FD")
        print_type = PrintType::HB_FD;
      else if (s == "HB_IC")
        print_type = PrintType::HB_IC;
      else if (s == "HB_STARTUP")
        print_type = PrintType::HB_STARTUP;
      else if (s == "HOMOTOPY")
        print_type = PrintType::HOMOTOPY;
      else if (s == "MPDE")
        print_type = PrintType::MPDE;
      else if (s == "MPDE_IC")
        print_type = PrintType::MPDE_IC;
      else if (s == "SENS")
        print_type = PrintType::SENS;
      else if (s == "TRAN")
        print_type = PrintType::TRAN;
      else
      {
        Report::DevelFatal0() << "Unrecognized analysis type " << s;
      }
    }
    else if (iterParam->tag() == "PRECISION")
    {
      print_parameters.streamPrecision_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "FILTER")
    {
      print_parameters.filter_ = iterParam->getImmutableValue<double>();
    }
    else if (iterParam->tag() == "FORMAT")
    {
      std::string s = iterParam->stringValue();
      if (s == "STD")
        format = Format::STD;
      else if (s == "TECPLOT")
        format = Format::TECPLOT;
      else if (s == "DAKOTA")
        format = Format::DAKOTA;
      else if (s == "PROBE")
        format = Format::PROBE;
      else if (s == "NOINDEX")
      {
        format = Format::STD;
        no_index = true;
        print_parameters.printIndexColumn_ = false;
      }
      else if (s == "CSV")
      {
        format = Format::CSV;
        no_index = true;
        print_parameters.printIndexColumn_ = false;
        print_parameters.delimiter_ = ",";
      }
      else if (s == "RAW")
      {
        format = defaultPrintParameters_.asciiRaw_ ? Format::RAW_ASCII : Format::RAW;
      }
      else
      {
        Report::DevelFatal0() << "Unrecognized print format " << s;
      }
      print_parameters.format_ = format;
    }
    else if (iterParam->tag() == "TIMEWIDTH")
    {
      print_parameters.timeWidth_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "TIMESCALEFACTOR")
    {
      print_parameters.outputTimeScaleFactor_ = iterParam->getImmutableValue<double>();
    }
    else if (iterParam->tag() == "FILE")
    {
      // netListFilename_ should be the default unless FILE was set to
      // something other than "" which is its default value.
      if (iterParam->stringValue() != "")
      {
        print_parameters.filename_ = iterParam->stringValue();
      }
    }
    else if (iterParam->tag() == "DELIMITER")
    {
      if (iterParam->stringValue() == "TAB")
      {
        print_parameters.delimiter_ = "\t";
      }
      else if (iterParam->stringValue() == "COMMA")
      {
        print_parameters.delimiter_ = ",";
      }
      else if (iterParam->stringValue() != "")
      {
        Report::UserWarning0() << "Invalid value of DELIMITER in .PRINT statment, ignoring";
      }
    }
    else
    {
      // This must be the first print variable.
      break;
    }
  }

  // Remaining parameters are variables
  print_parameters.variableList_.assign(iterParam, print_block.getParams().end());

  // Increase width to make sure that the columns do not run together
  if (print_parameters.streamWidth_ - print_parameters.streamPrecision_  < 9)
    print_parameters.streamWidth_ = print_parameters.streamPrecision_ + 9;

  // Indicate if an index column should be added
  print_parameters.printIndexColumn_ = !no_index
                            && print_parameters.format_ != Format::PROBE
                            && print_parameters.format_ != Format::TECPLOT
                            && print_parameters.format_ != Format::DAKOTA
                            && print_parameters.format_ != Format::RAW
                            && print_parameters.format_ != Format::RAW_ASCII;

  // Assemble the apropriate flavors of output variable lists based on the PRINT type and format
  if (print_type == PrintType::AC)
  {
    PrintParameters freq_print_parameters = print_parameters;
    if (freq_print_parameters.format_ != Format::PROBE)
    {
      freq_print_parameters.variableList_.push_front(Util::Param("FREQUENCY", 0.0));
    }
    if (freq_print_parameters.printIndexColumn_)
    {
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    }
    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;
    addOutputPrintParameters(OutputType::AC, freq_print_parameters);

    PrintParameters ac_ic_print_parameters = print_parameters;
    ac_ic_print_parameters.fallback_ = true;
    if (ac_ic_print_parameters.format_ != Format::PROBE)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    }
    if (ac_ic_print_parameters.printIndexColumn_)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    }
    addOutputPrintParameters(OutputType::AC_IC, ac_ic_print_parameters);
  }
  else if (print_type == PrintType::AC_IC)
  {
    PrintParameters ac_ic_print_parameters = print_parameters;
    if (ac_ic_print_parameters.format_ != Format::PROBE)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    }
    if (ac_ic_print_parameters.printIndexColumn_)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    }
    addOutputPrintParameters(OutputType::AC_IC, ac_ic_print_parameters);
  }
  else if (print_type == PrintType::HOMOTOPY)
  {
    PrintParameters homotopy_print_parameters = print_parameters;
    if (homotopy_print_parameters.format_ == Format::STD)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.prn";
    }
    else if (homotopy_print_parameters.format_ == Format::TECPLOT)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.dat";
    }
    else if (homotopy_print_parameters.format_ == Format::PROBE)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csd";
    }
    addOutputPrintParameters(OutputType::HOMOTOPY, homotopy_print_parameters);
  }
  else if (print_type == PrintType::SENS)
  {
    PrintParameters sensitivity_print_parameters = print_parameters;

    if (sensitivity_print_parameters.format_ == Format::STD)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
    }
    else if (sensitivity_print_parameters.format_ == Format::TECPLOT)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.dat";
    }
    else if (sensitivity_print_parameters.format_ == Format::PROBE)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.csd";
    }
    else if (sensitivity_print_parameters.format_ == Format::DAKOTA)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.txt";
    }

    std::copy(sensitivityVariableList_.begin(), sensitivityVariableList_.end(), std::back_inserter(sensitivity_print_parameters.variableList_));

    addOutputPrintParameters(OutputType::SENS, sensitivity_print_parameters);
  }
  else if (print_type == PrintType::HB)
  {
    PrintParameters freq_print_parameters = print_parameters;
    freq_print_parameters.fallback_ = true;
    if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(Util::Param("FREQUENCY", 0.0));
    if (freq_print_parameters.printIndexColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;
    addOutputPrintParameters(OutputType::HB_FD, freq_print_parameters);

    PrintParameters time_print_parameters = print_parameters;
    time_print_parameters.fallback_ = true;
    if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (time_print_parameters.printIndexColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_TD, time_print_parameters);

    PrintParameters hb_ic_print_parameters = print_parameters;
    hb_ic_print_parameters.fallback_ = true;
    if (hb_ic_print_parameters.format_ == Format::STD)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.prn";
    else if (hb_ic_print_parameters.format_ == Format::CSV)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.csv";
    else if (hb_ic_print_parameters.format_ == Format::TECPLOT)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.dat";
    if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_ic_print_parameters.printIndexColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_IC, hb_ic_print_parameters);

    PrintParameters hb_startup_print_parameters = print_parameters;
    hb_startup_print_parameters.fallback_ = true;
    if (hb_startup_print_parameters.format_ == Format::STD)
      hb_startup_print_parameters.defaultExtension_ = ".startup.prn";
    else if (freq_print_parameters.format_ == Format::CSV)
      hb_startup_print_parameters.defaultExtension_ = ".startup.csv";
    else if (freq_print_parameters.format_ == Format::TECPLOT)
      hb_startup_print_parameters.defaultExtension_ = ".startup.dat";
    if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_startup_print_parameters.printIndexColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_STARTUP, hb_startup_print_parameters);
  }
  else if (print_type == PrintType::HB_TD)
  {
    PrintParameters time_print_parameters = print_parameters;
    if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (time_print_parameters.printIndexColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_TD, time_print_parameters);
  }
  else if (print_type == PrintType::HB_FD)
  {
    PrintParameters freq_print_parameters = print_parameters;
    if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(Util::Param("FREQUENCY", 0.0));
    if (freq_print_parameters.printIndexColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;
    addOutputPrintParameters(OutputType::HB_FD, freq_print_parameters);
  }
  else if (print_type == PrintType::HB_IC)
  {
    PrintParameters hb_ic_print_parameters = print_parameters;
    if (hb_ic_print_parameters.format_ == Format::STD)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.prn";
    else if (hb_ic_print_parameters.format_ == Format::CSV)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.csv";
    else if (hb_ic_print_parameters.format_ == Format::TECPLOT)
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.dat";
    if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_ic_print_parameters.printIndexColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_IC, hb_ic_print_parameters);
  }
  else if (print_type == PrintType::HB_STARTUP)
  {
    PrintParameters hb_startup_print_parameters = print_parameters;
    if (hb_startup_print_parameters.format_ == Format::STD)
      hb_startup_print_parameters.defaultExtension_ = ".startup.prn";
    else if (hb_startup_print_parameters.format_ == Format::CSV)
      hb_startup_print_parameters.defaultExtension_ = ".startup.csv";
    else if (hb_startup_print_parameters.format_ == Format::TECPLOT)
      hb_startup_print_parameters.defaultExtension_ = ".startup.dat";
    if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_startup_print_parameters.printIndexColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::HB_STARTUP, hb_startup_print_parameters);
  }
  else if (print_type == PrintType::MPDE)
  {
    PrintParameters mpde_print_parameters = print_parameters;
    addOutputPrintParameters(OutputType::MPDE, mpde_print_parameters);

    PrintParameters mpde_ic_print_parameters = print_parameters;
    mpde_ic_print_parameters.fallback_ = true;
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";

    if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.printIndexColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::MPDE_IC, mpde_ic_print_parameters);
  }
  else if (print_type == PrintType::MPDE_IC)
  {
    PrintParameters mpde_ic_print_parameters = print_parameters;
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";

    if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.printIndexColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::MPDE_IC, mpde_ic_print_parameters);
  }
  else if (print_type == PrintType::TRAN)
  {
    addOutputPrintParameters(OutputType::TRAN, print_parameters);
    PrintParameters homotopy_print_parameters = print_parameters;
    homotopy_print_parameters.fallback_ = true;
    if (homotopy_print_parameters.format_ == Format::STD)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.prn";
    }
    else if (homotopy_print_parameters.format_ == Format::TECPLOT)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.dat";
    }
    else if (homotopy_print_parameters.format_ == Format::PROBE)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csd";
    }
    addOutputPrintParameters(OutputType::HOMOTOPY, homotopy_print_parameters);

    PrintParameters mpde_print_parameters = print_parameters;
    mpde_print_parameters.fallback_ = true;
    addOutputPrintParameters(OutputType::MPDE, mpde_print_parameters);

    PrintParameters mpde_ic_print_parameters = print_parameters;
    mpde_ic_print_parameters.fallback_ = true;
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";

    if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.printIndexColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    addOutputPrintParameters(OutputType::MPDE_IC, mpde_ic_print_parameters);
  }
  else if (print_type == PrintType::DC)
  {
    PrintParameters dc_print_parameters = print_parameters;
    if (dc_print_parameters.format_ == Format::STD)
      dc_print_parameters.defaultExtension_ = ".prn";
    addOutputPrintParameters(OutputType::DC, dc_print_parameters);
  }
  else
  {
    Report::UserError0() << "Unrecognized .PRINT type";
  }

  activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerLoad
// Purpose       : registers set of variables to set for .LOAD
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/11/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerLoad(const Util::OptionBlock & OB)
{
  loadFlag_ = true;
  Report::UserWarning0()
    << ".LOAD not supported yet.  Use .INCLUDE instead";
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSens
// Purpose       : registers set of variables to set for .SENS.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
bool OutputMgr::registerSens(const Util::OptionBlock &option_block)
{
  bool bsuccess = true;

  std::vector<std::string> functions;
  std::vector<std::string> parameters;

  for (std::list<Util::Param>::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end(); it != end; ++it)
  {
    std::string tag = (*it).uTag();
    if ((*it).uTag() == "OBJFUNC")
    {
      functions.push_back((*it).stringValue());
    }
    else if (std::string((*it).uTag(), 0, 5) == "PARAM")
    {
      parameters.push_back((*it).stringValue());
    }
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  for (std::vector<std::string>::const_iterator it1 = functions.begin(), end1 = functions.end(); it1 != end1; ++it1)
  {
    const std::string &function = (*it1);

    {
      Util::Marshal mout;
      mout << function << std::string("FUNCTION") << Util::Op::identifier<SensitivityObjFunctionOp>() << int(0);
      sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
    }

    int index = 0;
    for (std::vector<std::string>::const_iterator it2 = parameters.begin(), end2 = parameters.end(); it2 != end2; ++it2, ++index)
    {
      const std::string &parameter = (*it2);

      if (sensitivityOptions_ & SensitivityOptions::DIRECT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::SCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
      }

      if (sensitivityOptions_ & SensitivityOptions::ADJOINT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::DIRECT)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
      }
    }
  }

  enableSensitivityFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSensOptions
// Purpose       : registers set of variables to set for .OPTIONS SENSITIVITY
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
bool OutputMgr::registerSensOptions(const Util::OptionBlock &option_block)
{
  sensitivityOptions_ = 0;

  bool adjointGiven=false;
  bool outputUnscaledGiven=false;

  for (std::list<Util::Param>::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end() ; it != end; ++it)
  {
    if ((*it).uTag() == "ADJOINT")
    {
      adjointGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::ADJOINT;
      }
    }
    else if ((*it).uTag() == "DIRECT" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::DIRECT;
    }
    else if ((*it).uTag() == "OUTPUTSCALED" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::SCALED;
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::UNSCALED;
      }
    }
  }

  // default behavior is for the code to assume adjoint sensitivites are wanted
  // even if ADJOINT=1 wasn't specified.
  if (!adjointGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::ADJOINT;
  }

  // default behavior is for the code to assume unscaled sensitivities are wanted
  // for output.
  if (!outputUnscaledGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::UNSCALED;
  }

  return true;
}

//-----------------------------------------------------------------------------
template<class It, class Out>
void makeParamList(const std::string name, It first, It last, Out out)
{
  for (; first != last; ++first)
  {
    *out++ = Util::Param(name, 1.0);
    *out++ = Util::Param(*first, 0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::removeStarVariables
// Purpose       : Process V(*) and I(*) on .print line
// Special Notes : Replaces v(*) and i(*) that appears on the .print line
//                 with a list of all v() or i() variables as appropriate.
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void removeStarVariables(Parallel::Machine comm, ParameterList &variable_list,
    const NodeNamePairMap &all_nodes, const NodeNamePairMap &external_nodes)
{
  bool vStarFound = false;
  bool iStarFound = false;

  ParameterList::iterator vStarPosition = variable_list.end();
  ParameterList::iterator iStarPosition = variable_list.end();

  // remove the v(*) and i(*) from the .print line
  for (ParameterList::iterator it = variable_list.begin(); it != variable_list.end(); )
  {
    // process the * entries
    if ((*it).tag() == "*")
    {
      // move to the type
      --it;

      // remember type of replacement and location of *
      if ((*it).tag() == "V")
      {
        vStarFound = true;
        vStarPosition = it;
        --vStarPosition;
      }
      else if ((*it).tag() == "I")
      {
        iStarFound = true;
        iStarPosition = it;
        --iStarPosition;
      }

      // remove the v or i
      it = variable_list.erase(it);

      // remove the *
      it = variable_list.erase(it);
    }

    // move to next item on .print line
    else
    {
      ++it;
    }
  }

  std::vector<std::string> v_list;
  std::vector<std::string> i_list;

  if (vStarFound)
  {
    for (NodeNamePairMap::const_iterator it = external_nodes.begin(); it != external_nodes.end() ; ++it)
    {
      ExtendedString tmpStr((*it).first);
      tmpStr.toUpper();

      if (tmpStr.rfind("BRANCH") == std::string::npos)
        v_list.push_back(tmpStr);
    }
  }

  if (iStarFound)
  {
    for (NodeNamePairMap::const_iterator iter_a = all_nodes.begin(); iter_a != all_nodes.end() ; ++iter_a)
    {
      ExtendedString tmpStr((*iter_a).first);
      tmpStr.toUpper();

      size_t pos = tmpStr.rfind("BRANCH");
      if (pos != std::string::npos && tmpStr[0] != 'Y')
      {
        tmpStr = tmpStr.substr(0, pos - 1);

        i_list.push_back(tmpStr);
      }
    }
  }

  Util::Marshal mout;
  mout << v_list << i_list;

  std::vector<std::string> dest;
  Parallel::AllGatherV(comm, mout.str(), dest);

  ParameterList vStarList;
  ParameterList iStarList;

  for (int p = 0; p < Parallel::size(comm); ++p)
  {
    Util::Marshal min(dest[p]);

    std::vector<std::string> x;
    std::vector<std::string> y;
    min >> x >> y;
    makeParamList("V", x.begin(), x.end(), std::back_inserter(vStarList));
    makeParamList("I", y.begin(), y.end(), std::back_inserter(iStarList));
  }

  // advance iterators
  ++vStarPosition;
  ++iStarPosition;

  // append temporary lists to print block, erasing temporary lists
  variable_list.splice(vStarPosition, vStarList);
  variable_list.splice(iStarPosition, iStarList);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setSweepParameters
// Purpose       : Copy .DC or .STEP sweep parameters, and set up flags if
//                 sweeping the TEMP variable.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void OutputMgr::setStepSweep(const Analysis::SweepVector &step_sweep_parameters)
{
  if (!step_sweep_parameters.empty())
  {
    outputState_.stepSweepVector_ =  step_sweep_parameters;
  }

  // check if one of the sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!step_sweep_parameters.empty())
  {
    for (Analysis::SweepVector::const_iterator iterP = outputState_.stepSweepVector_.begin(), end = outputState_.stepSweepVector_.end(); iterP != end; ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        tempSweepFlag_ = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setSweepParameters
// Purpose       : Copy .DC or .STEP sweep parameters, and set up flags if
//                 sweeping the TEMP variable.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void OutputMgr::setDCSweep(const Analysis::SweepVector &dc_sweep_parameters)
{
  if (!dc_sweep_parameters.empty())
  {
    outputState_.dcSweepVector_ = dc_sweep_parameters;
  }

  // check if one of the DC sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!dc_sweep_parameters.empty() && !outputCalledBefore_)
  {
    for (Analysis::SweepVector::const_iterator iterP = outputState_.dcSweepVector_.begin();iterP != outputState_.dcSweepVector_.end(); ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        tempSweepFlag_ = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::fixupPrintParameters
// Purpose       : Perform some .print line checks and munging, primarily
//                 dealing with use of V(*) and I(*)
// Special Notes : Mostly a re-wrapping of work that used to be done in the
//                 check_output function, which was too sprawling to remain
//                 as a single function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/26/2013
//-----------------------------------------------------------------------------
void OutputMgr::fixupPrintParameters(Parallel::Machine comm, PrintParameters &print_parameters)
{
  removeStarVariables(comm, print_parameters.variableList_, allNodes_, externNodes_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::output
// Purpose       : Runs specified output commands
//
// Special Notes : ERK.  I've refactored this so that it receives STEP
//                 and DC sweep parameter information(mainly name and value).
//                 This function is called from the time integrator, which has
//                 all this info.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::output(
  Parallel::Machine             comm,
  const double                  time,
  const int                     stepNumber,
  const int                     maxStep,
  const Analysis::SweepVector & step_sweep_vector,
  const int                     dcNumber,
  const int                     maxDC,
  const Analysis::SweepVector & dc_sweep_vector,
  const N_LAS_Vector &          solnVecPtr,
  const N_LAS_Vector &          stateVecPtr,
  const N_LAS_Vector &          storeVecPtr,
  const std::vector<double> &   objectiveVec,
  const std::vector<double> &   dOdpVec,
  const std::vector<double> &   dOdpAdjVec,
  const std::vector<double> &   scaled_dOdpVec,
  const std::vector<double> &   scaled_dOdpAdjVec,
  bool                          skipPrintLineOutput)
{
  // copy over time:
  outputState_.circuitTime_ = time;

  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;
  if (outputState_.stepMaxCount_ > 0)
  {
    STEPEnabledFlag_ = true;
  }

  // copy the new values into the locally owned vector:
  if (!step_sweep_vector.empty())
  {
    outputState_.stepSweepVector_ = step_sweep_vector;
  }

  // copy over the dc sweep information:
  dcLoopNumber_ = dcNumber;
  maxDCSteps_ = maxDC;

  if (!analysisManager_->getBlockAnalysisFlag())
  {
    // This error test should not be used in the MPDE case, as at least
    // one of the initial conditions that can be
    // used by MPDE is a DC sweep.
    if (maxDCSteps_ > 0 && analysisManager_->getTransientFlag())
    {
      Report::DevelFatal0() << "Print type is inconsistent, maxDCSteps = " << maxDCSteps_;
    }
  }

  if (!dc_sweep_vector.empty())
  {
    outputState_.dcSweepVector_ = dc_sweep_vector;

    // For now just have PRINTdcvalue, etc just be the first parameter.
    const Analysis::SweepParam &firstParam = dc_sweep_vector.front();
    PRINTdcname_   = firstParam.name;
    PRINTdcvalue_  = firstParam.currentVal;
    if (firstParam.type == "LIST")
    {
      PRINTdcstart_  = firstParam.valList.front(); // [0];
      PRINTdcstop_   = firstParam.valList.back(); // [size1-1];
      // int size1 = firstParam.valList.size();
      // PRINTdcstop_   = firstParam.valList[size1-1];
    }
    else
    {
      PRINTdcstart_  = firstParam.startVal;
      PRINTdcstop_   = firstParam.stopVal;
    }
  }

  // Check for temperature:
  outputState_.circuitTemp_ = deviceInterface_->getParamAndReduce("TEMP");

  // call on the measure manager to update any active measures
  // update .measure functions before outputting the .print line
  if (analysisManager_->getTransientFlag())
  {
    measureManager_->updateTranMeasures(comm, outputState_.circuitTime_, &solnVecPtr, &stateVecPtr, &storeVecPtr);
    fourierManager_->updateFourierData(comm, outputState_.circuitTime_, &solnVecPtr);
  }
  else
  {
    measureManager_->updateDCMeasures(comm, outputState_.dcSweepVector_, &solnVecPtr, &stateVecPtr, &storeVecPtr);
  }

  // Needs to pass skipPrintLineOutput
  if (!skipPrintLineOutput)
  {
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator
          it = activeOutputterStack_.back().begin();
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->output(comm, solnVecPtr, stateVecPtr, storeVecPtr);
        (*it)->outputSensitivity(comm, objectiveVec,
                                 dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec,
                                 solnVecPtr, stateVecPtr, storeVecPtr);
      }
    }
  }

  // This will become an outputter
  if (hdf5FileNameGiven_)
  {
    updateHDF5Output(comm, solnVecPtr);
  }

  // transient or dc values for the objective function call.
  double arg1 = 0.0;
  double arg2 = 0.0;
  if (analysisManager_->getTransientFlag())
  {
    arg1 = outputState_.circuitTime_;
  }
  else
  {
    if (outputState_.dcSweepVector_.size() > 0)
    {
      arg1 = outputState_.dcSweepVector_[0].currentVal;
    }
    if (outputState_.dcSweepVector_.size() > 1)
    {
      arg2 = outputState_.dcSweepVector_[1].currentVal;
    }
  }

  // if there are any simulation level functions(objective objects)
  // that need data from this time step, send it to them as well
  for (ObjectiveMap::iterator ob = objectiveMap_.begin(); ob != objectiveMap_.end(); ++ob)
  {
    (*ob).second.setup(comm, *this);
    (*ob).second.save(comm, outputState_.circuitTime_, arg1, arg2, &solnVecPtr, &stateVecPtr, &storeVecPtr);
  }

  outputCalledBefore_ = true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputAC(
  Parallel::Machine     comm,
  double                frequency,
  const N_LAS_Vector &  real_solution_vector,
  const N_LAS_Vector &  imaginary_solution_vector)
{
  outputState_.circuitFrequency_ = frequency;

  measureManager_->updateACMeasures(comm, frequency, &real_solution_vector, &imaginary_solution_vector);
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator
        it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputAC(comm, frequency, real_solution_vector, imaginary_solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHomotopy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const N_LAS_Vector &                  solution_vector)
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator
        it = activeOutputterStack_.back().begin();
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHomotopy(comm, parameter_names, param_values, solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::finishOutput
// Purpose       : Runs specified finish output commands
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 12/13/00
//-----------------------------------------------------------------------------
void OutputMgr::finishOutput()
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator
        it = activeOutputterStack_.back().begin();
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->finishOutput();
    }
  }

  if (hdf5FileNameGiven_)
  {
    closeHDF5Output();
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHB
// Purpose       : .print output for Harmonic Balance runs
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/5/2013
//-----------------------------------------------------------------------------
void OutputMgr::outputHB(
  Parallel::Machine                     comm,
  const int                             stepNumber,
  const int                             maxStep,
  const Analysis::SweepVector &         step_sweep_parameters,
  const std::vector<double> &           timePoints,
  const std::vector<double> &           freqPoints,
  const N_LAS_BlockVector &             timeDomainSolutionVec,
  const N_LAS_BlockVector &             freqDomainSolutionVecReal,
  const N_LAS_BlockVector &             freqDomainSolutionVecImaginary,
  const N_LAS_BlockVector &             timeDomainStoreVec,
  const N_LAS_BlockVector &             freqDomainStoreVecReal,
  const N_LAS_BlockVector &             freqDomainStoreVecImaginary)
{
  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;
  if (outputState_.stepMaxCount_ > 0)
  {
    STEPEnabledFlag_ = true;
  }

  // copy the new values into the locally owned vector:
  if (!(step_sweep_parameters.empty()))
  {
    outputState_.stepSweepVector_ = step_sweep_parameters;
  }

  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin();
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHB(comm, timePoints, freqPoints,
                      timeDomainSolutionVec, freqDomainSolutionVecReal, freqDomainSolutionVecImaginary,
                      timeDomainStoreVec, freqDomainStoreVecReal, freqDomainStoreVecImaginary);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputMPDE
// Purpose       : .PRINT output for mpde runs
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/22/03
//-----------------------------------------------------------------------------
void OutputMgr::outputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points,
  const N_LAS_Vector &          solution_vector)
{
  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputMPDE(comm, time, fast_time_points, solution_vector);
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::steppingComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::startStep(
    int                           step,
    int                           max_step)
{
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin();
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->startStep(step, max_step);
      }
    }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::steppingComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::steppingComplete()
{
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin();
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->steppingComplete();
      }
    }
}

// //-----------------------------------------------------------------------------
// // Function      : getMeasureValue
// // Purpose       : get a .measure value
// // Special Notes :
// // Creator       : Dave Shirley, PSSI
// // Creation Date : 4/29/10
// //-----------------------------------------------------------------------------
// void OutputMgr::getMeasureValue(const std::string &name, double &result, bool &found) const
// {
//   found = false;
//   measureManager_.getMeasureValue(name, result, found);
// }

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputDCOP
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgr::outputDCOP(
  Parallel::Machine     comm,
  const N_LAS_Vector &  solution)
{
  if (outputOnceAlreadyFlag_)
  {
    return;
  }

  if (icData_.output_op_)
  {
    std::ofstream opOut;

    for (NodeNamePairMap::iterator it = allNodes_.begin(), end = allNodes_.end(); it != end ; ++it)
    {
      int ind = (*it).second.first;
      (*it).second.second = solution[ind];
    }

    if (Parallel::rank(comm) == 0)
    {
      std::string outputFileName;

      if (icData_.output_op_file_ == "")
        outputFileName = netListFilename_ + ".op";
      else
        outputFileName = icData_.output_op_file_;
      opOut.open(outputFileName.c_str());
    }

    writeDCOPRestart(comm, opOut, allNodes_);

    if (Parallel::rank(comm) == 0)
      opOut.close();
  }
  else if (icData_.saveFlag_)
  {
    std::ofstream saveOutputStream;

    if (Parallel::rank(comm) == 0)
    {
      std::string outputFileName;
      if (icData_.saveOutputFile_  == "")
      {
        outputFileName = netListFilename_ + ".ic";
      }
      else
      {
        outputFileName = icData_.saveOutputFile_;
      }

      saveOutputStream.open(outputFileName.c_str());
    }

    // check for .IC
    outputIC_or_NODESET(comm, saveOutputStream, icData_.saveFileType_, allNodes_, solution);
    if (Parallel::rank(comm) == 0)
      saveOutputStream.close();
  }

  outputOnceAlreadyFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setupInitialConditions
// Purpose       : This function is an umbrella function, under which the
//                 various types of initial conditions are potentially set up.
//                 These include, but aren't limited to .IC, .NODESET and
//                 .DCOP restart.
//
//                 In addition to setting these things up, it does some
//                 nominal checking to make sure that more than on IC hasn't
//                 been specified(ie DCOP restart and .IC can't both be set).
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool
OutputMgr::setupInitialConditions(
  Parallel::Machine     comm,
  N_LAS_Vector &        solnVec,
  N_LAS_Vector &        flagVec)
{
    bool dcopRestart = false;
    bool dotIC = false;
    bool NODESET = false;
    icData_.icType_ = InitialConditionsData::IC_TYPE_UNDEFINED;

    if (icData_.input_op_ && icData_.ICflag_)
    {
      Report::UserFatal0()
        << "Cannot set both DCOP restart and .IC simultaneously.";
    }

    if (icData_.input_op_ && icData_.nodesetflag_)
    {
      Report::UserFatal0()
        << "Cannot set both DCOP restart and .NODESET simultaneously.";
    }

    if (icData_.ICflag_ && icData_.nodesetflag_)
    {
      Report::UserFatal0()
        << "Cannot set both .IC and .NODESET simultaneously.";
    }

    if (icData_.input_op_)
    {
      std::string inputFileName;
      std::ifstream opIn;

      if (icData_.input_op_file_ == "")
      {
        inputFileName = netListFilename_ + ".op";
      }
      else
      {
        inputFileName = icData_.input_op_file_;
      }
      opIn.open(inputFileName.c_str(), std::ios::in);

      if (opIn.good())
      {
        Xyce::dout() << "Reading in operating point initial estimate from: "
          << inputFileName << std::endl;

        // check for dcop restart(the original)
        dcopRestart = readDCOPRestart(comm, opIn, allNodes_, solnVec, flagVec, icData_.opData_, icData_.op_found_, icData_.total_soln_);
        if (dcopRestart)
          icData_.icType_ = InitialConditionsData::IC_TYPE_DCOP_RESTART;
      }
    }
    else if (icData_.ICflag_)
    {
      // check for .IC
      if (icData_.ICflag_)
        dotIC = setupIC_or_NODESET(comm, allNodes_, solnVec, flagVec, InitialConditionsData::IC_TYPE_IC, ICblockVec_, icData_.opData_, icData_.op_found_, icData_.total_soln_);
      if (dotIC)
        icData_.icType_ = InitialConditionsData::IC_TYPE_IC;
    }
    else if (icData_.nodesetflag_)
    {
      // check for .NODESET
      if (icData_.nodesetflag_)
        NODESET = setupIC_or_NODESET(comm, allNodes_, solnVec, flagVec, InitialConditionsData::IC_TYPE_NODESET, nodesetblockVec_, icData_.opData_, icData_.op_found_, icData_.total_soln_);
      if (NODESET)
        icData_.icType_ = InitialConditionsData::IC_TYPE_NODESET;
    }

    return dcopRestart || dotIC || NODESET;
}

} // namespace IO
} // namespace Xyce
