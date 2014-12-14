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
// Filename      : NetlistImportTool.C
//
// Purpose       : Implement the interface to to read and parse a netlist for
//                 an electrical circuit.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.99.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <N_UTL_fwd.h>

#include <N_IO_NetlistImportTool.h>

#include <N_IO_DistributionTool.h>

#include <N_PDS_ParallelMachine.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

#include <N_IO_DeviceBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OutputMgr.h>

#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInterface.h>

#include <N_IO_CmdParse.h>

#include <N_ERH_ErrorMgr.h>

#include <N_TOP_NodeBlock.h>
#include <N_TOP_NodeDevBlock.h>
#include <N_TOP_InsertionTool.h>
#include <N_TOP_TopologyMgr.h>
#include <N_DEV_SourceData.h>
#include <N_UTL_Xyce.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace IO {

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::factory
// Purpose       : This function returns a pointer to the only
//                 instance of this class. (for the current library
//                 instance, hopefully)
//
// Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
//                 static pointer was returned) but had to be changed
//                 so that the library version of Xyce would work
//                 correctly.
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool* NetlistImportTool::factory(IO::CmdParse & cp,
    Xyce::Topo::Manager & tm)
{
  NetlistImportTool * NIT_ptr = new NetlistImportTool(cp,tm);
  return NIT_ptr;
}



//-------------------------------------------------------------------------
// Function      : NetlistImportTool::~NetlistImportTool
// Purpose       : Destructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::~NetlistImportTool()
{
    delete circuitBlock_;
    delete distToolPtr_;
}

//-----------------------------------------------------------------------------
// Function      : NetlistImportTool::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  return true;
}

typedef std::list<Util::Param> ParameterList;

namespace { // <unnamed>

void setLeadCurrentDevices(const ParameterList &variable_list, Device::DeviceInterface &device_manager)
{
  std::set<std::string> devicesNeedingLeadCurrents;

  for (ParameterList::const_iterator iterParam = variable_list.begin() ; iterParam != variable_list.end(); ++iterParam)
  {
    std::string varType(iterParam->tag());

    if (Util::hasExpressionTag(*iterParam))
    {
      std::vector<std::string> leads;

      Util::Expression exp;
      exp.set(iterParam->tag());
      exp.get_names(XEXP_LEAD, leads);

      // any lead currents found in this expression need to be communicated to the device manager.
      // Multi terminal devices have an extra designator on the name as in name{lead_name}
      // need to remove any {} in the name.
      for (std::vector<std::string>::const_iterator currLeadItr = leads.begin(); currLeadItr != leads.end(); ++currLeadItr)
      {
        size_t leadDesignator = currLeadItr->find_first_of("{");
        devicesNeedingLeadCurrents.insert( currLeadItr->substr(0, leadDesignator));
      }
    }
    else
    {
      if ((varType == "I" || (varType.size() == 2 && varType[0] == 'I')) &&  (iterParam->getImmutableValue<int>() > 0))
      {
        // any devices found in this I(xxx) structure need to be communicated to the device manager
        // so that the lead currents can be calculated
        if (iterParam->getImmutableValue<int>() != 1)
        {
          Report::UserError0() << "Only one device argument allowed in I() in .print";
        }
        else
        {
          ++iterParam;
          if (varType.size() == 2)
          {
            devicesNeedingLeadCurrents.insert(iterParam->tag());
          }
          else
          {
            devicesNeedingLeadCurrents.insert(iterParam->tag());
          }
        }
      }
    }
  }

  // the list of devices that need lead currents.
  if ( Xyce::DEBUG_IO && !devicesNeedingLeadCurrents.empty())
  {
    std::set<std::string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrents.begin();
    std::set<std::string>::iterator endDeviceNameItr = devicesNeedingLeadCurrents.end();
    Xyce::dout() << "Devices for which lead currents were requested: ";
    while ( currentDeviceNameItr != endDeviceNameItr)
    {
      Xyce::dout() << *currentDeviceNameItr << "  ";
      currentDeviceNameItr++;
    }
  }

  // This is the last call before devices are constructed
  // So it's the last time I can isolate lead currents. However it won't be sufficient
  // as we haven't parsed expressions in devices yet.  I'll need to rethink
  // when this call is made to the device manager.  Sometime just after device instance
  // construction but before the store vector is allocated.
  device_manager.setLeadCurrentRequests(devicesNeedingLeadCurrents);
}

} // namespace <unnamed>

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::constructCircuitFromNetlist
// Purpose       : Construct a circuit from a netlist.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
int NetlistImportTool::constructCircuitFromNetlist(const std::string &netlistFile,
                                                   const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
                                                   OutputMgr &output_manager, Measure::Manager &measure_manager)
{
  std::map<std::string,RefCountPtr<Device::InstanceBlock> > dNames;
  std::set<std::string> nNames;

  // build metadata
  metadata_.buildMetadata();

  // Create a circuitBlock instance to hold netlist circuit data.
  circuitBlock_ = new IO::CircuitBlock(
    netlistFile,
    commandLine_,
    metadata_,
    modelNames_ ,
    ssfMap_,
    circuitContext_,
    useCount_,
    dNames,
    nNames,
    aliasNodeMap_, 
    externalNetlistParams );

  // Parse the netlist file.
  if (Xyce::DEBUG_IO)
    Xyce::dout() << "Starting netlist parsing." << std::endl;

  // Get the device insertion tool. Note: the argument determines
  // the insertion tool returned, since only the device insertion
  // is available at this point (08/27/2003) only it is available
  // and the argument is does nothing. This may need to be changed
  // later.

  Xyce::Topo::InsertionTool* insertToolPtr = topMgr_.getInsertionTool("");
  circuitBlock_->registerInsertionTool(insertToolPtr);

  // Regsiter the device interface with the circuit block.
  circuitBlock_->registerDeviceInterface(devIntPtr_);

  distToolPtr_ = new IO::DistributionTool(circuitBlock_, commandLine_, pdsCommPtr_);

  distToolPtr_->registerPkgOptionsMgr(pkgOptMgrPtr_);

  circuitBlock_->registerDistributionTool(distToolPtr_);

  // Read in the circuit context and circuit options on the root processor.
  // NOTE:  The circuit context and circuit options are registered with the distTool.
  if (Parallel::rank(comm_) == 0)
  {
    circuitBlock_->parseNetlistFilePass1();
  }

  // Just in case an error was reported in parsing the netlist. 
  N_ERH_ErrorMgr::safeBarrier(comm_);
  
  // Distribute the circuit context and circuit options to all other processors.
  bool ok = distToolPtr_->broadcastGlobalData();

  if (Parallel::rank(comm_) == 0)
  {
    if (ok)
      ok = circuitBlock_->parseNetlistFilePass2();

  }
  else
  {
    // Wait to receive the circuit data from processor 0.
    if (ok)
      distToolPtr_->receiveCircuitData();
  }

  // Just in case an error was reported in parsing the netlist. 
  N_ERH_ErrorMgr::safeBarrier(comm_);

  // Check for name collisions between devices
  Xyce::IO::checkNodeDevConflicts(circuitBlock_, pdsCommPtr_);

  if (ok) {
    if (Xyce::DEBUG_IO)
      Xyce::dout() << "Completed netlist parsing. ";

    // Write out preprocessed netlist, if requested.
    if (Parallel::rank(comm_) == 0)
    {
      circuitBlock_->writeOutNetlist();
    }

    output_manager.setTitle(circuitBlock_->getTitle());

    setLeadCurrentDevices(output_manager.getVariableList(), *devIntPtr_);
    printLineDiagnostics(comm_, output_manager, measure_manager, nNames, dNames, aliasNodeMap_);
  }

  return 1;
}



//-------------------------------------------------------------------------
// Function      : NetlistImportTool::registerDevMgr
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
bool NetlistImportTool::registerDevMgr(Device::DeviceInterface * devPtr)
{
   return ( (devIntPtr_ = devPtr) != NULL );
}

//-----------------------------------------------------------------------------
// Function      : NetlistImportTool::registerParallelServices
// Purpose       : Registers N_PDS_Comm object for parallel communication.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/17/2003
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerParallelServices(N_PDS_Comm * tmp_pds_ptr)
{
  if( !tmp_pds_ptr )
    return false;

  pdsCommPtr_ = tmp_pds_ptr;
  comm_ = pdsCommPtr_->comm();

  if (Xyce::DEBUG_DISTRIBUTION)
    Xyce::dout() << "Number of Processors:   " << Parallel::size(comm_) << std::endl
                 << "Processor ID:           " << Parallel::rank(comm_) << std::endl;

  return true;
}

//-------------------------------------------------------------------------
// Name          : NetlistImportTool::NetlistImportTool
// Purpose       : Default constructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::NetlistImportTool(
    IO::CmdParse & cp,
    Xyce::Topo::Manager & tm)
  : circuitBlock_(NULL),
    commandLine_(cp),
    topMgr_(tm),
    currentContextPtr_(NULL),
    distToolPtr_(NULL),
    circuitContext_( metadata_, contextList_, currentContextPtr_ ),
    useCount_(0),
    pdsCommPtr_(0),
    comm_(0)
{}

} // namespace IO
} // namespace Xyce
