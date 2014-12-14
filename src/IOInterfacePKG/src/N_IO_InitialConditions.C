//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_IO_InitialConditions.C,v $
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
// Revision Number: $Revision: 1.4.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iterator>

#include <N_ERH_Message.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_Marshal.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : OutputMgr_DCOPOptionsReg
// Purpose       : functor for registering DCOP options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_DCOPOptionsReg : public PkgOptionsReg
{
  OutputMgr_DCOPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerDCOPOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_ICOptionsReg
// Purpose       : functor for registering IC options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_ICOptionsReg : public PkgOptionsReg
{
  OutputMgr_ICOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerIC ( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_NodeSetOptionsReg
// Purpose       : functor for registering NodeSet options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_NodeSetOptionsReg : public PkgOptionsReg
{
  OutputMgr_NodeSetOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerNodeSet( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_SaveOptionsReg
// Purpose       : functor for registering Save options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_SaveOptionsReg : public PkgOptionsReg
{
  OutputMgr_SaveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSave( options ); }

  OutputMgr &outputManager_;
};


bool OutputMgr::registerPkgOptionsMgrInitialConditions(PkgOptionsMgr &pkgOpt)
{
  pkgOpt.submitRegistration(
    "OP_IO", netListFilename_, new OutputMgr_DCOPOptionsReg(*this));

  pkgOpt.submitRegistration(
    "IC", netListFilename_, new OutputMgr_ICOptionsReg(*this));

  pkgOpt.submitRegistration(
    "NODESET", netListFilename_, new OutputMgr_NodeSetOptionsReg(*this));

  pkgOpt.submitRegistration(
    "SAVE", netListFilename_, new OutputMgr_SaveOptionsReg(*this));

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSave
// Purpose       : registers set of variables to set for .SAVE.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/11/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerSave(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << " SAVE Registered::" << std::endl;
#endif

  icData_.saveFlag_ = true;

  ExtendedString sval("");

  ParameterList::const_iterator iterPL = OB.getParams().begin();
  ParameterList::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
   Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif
    if (iterPL->tag() == "TYPE")
    {
      sval = iterPL->stringValue();
      sval.toUpper();

      if (sval == "NODESET" || sval == ".NODESET")
      {
        icData_.saveFileType_ = ".NODESET";
      }
      else if (sval == "IC" || sval == ".IC")
      {
        icData_.saveFileType_ = ".IC";
      }
      else
      {
        Report::UserWarning0() << "Unrecognized type specified on .SAVE command.  Defaulting to .NODESET";
      }
    }
    else if (iterPL->tag() == "FILE")
    {
      icData_.saveOutputFile_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else if (iterPL->tag() == "LEVEL")
    {
      sval = iterPL->stringValue();
      sval.toUpper();

      if (sval == "ALL")
      {
        // do nothing
      }
      else if (sval == "NONE")
      {
        // none means don't output anything.  Pretend the .SAVE line isn't in the netlist.
        icData_.saveFlag_ = false;
        icData_.saveFileLevel_ = "NONE";
      }
      else if (sval == "TOP")
      {
        Report::UserWarning0() << "LEVEL=TOP in .SAVE line not supported.  Ignoring. ";
      }
      else
      {
        Report::UserWarning0() << "Unrecognized LEVEL " << sval << " specified in .SAVE command";
      }
    }
    else
    {
      Report::UserWarning0() << "Parameter " << iterPL->tag() << " not recognized in .SAVE command";
    }


    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerIC
// Purpose       : registers set of variables to set for .IC or .DCVOLT
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerIC(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << " IC Registered::" << std::endl;
#endif

  icData_.ICflag_ = true;

  // save the ICblock_.  This will be needed later.
  // We allow multiple .IC blocks in the netlist, so there needs to be an
  // STL vector of these.
  ICblockVec_.push_back(OB);
  int end = ICblockVec_.size()-1;

#ifdef Xyce_DEBUG_IO
  ParameterList::iterator iterParam = ICblockVec_[end].getParams().begin();
  ParameterList::iterator endParam  = ICblockVec_[end].getParams().end();

  while (iterParam != endParam)
  {
   Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
     ++iterParam;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerNodeSet
// Purpose       : registers set of variables to set for .NODESET
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerNodeSet(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
 Xyce::dout() << " NODESET Registered::" << std::endl;
#endif

 icData_.nodesetflag_ = true;

  // save the nodesetblock_.  This will be needed later.
  // We allow multiple .nodeset blocks in the netlist, so there needs to be an
  // STL vector of these.
  nodesetblockVec_.push_back(OB);
  int end = nodesetblockVec_.size()-1;

#ifdef Xyce_DEBUG_IO
  ParameterList::iterator iterParam = nodesetblockVec_[end].getParams().begin();
  ParameterList::iterator endParam  = nodesetblockVec_[end].getParams().end();

  while (iterParam != endParam)
  {
   Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
     ++iterParam;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDCOPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/08/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerDCOPOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
  {
    if (iterPL->tag() == "INPUT")
    {
      icData_.input_op_ = true;
      icData_.input_op_file_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "OUTPUT")
    {
      icData_.output_op_ = true;
      icData_.output_op_file_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else
    {
      Report::UserWarning0() << "Parameter " << iterPL->tag() << " not recognized in .DCOP command";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getICData
// Purpose       : Provide OP data for LOCA
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/04/06
//-----------------------------------------------------------------------------
NodeNamePairMap & OutputMgr::getICData(int & op_found, int & icType)
{
  op_found = icData_.op_found_;
  icType = icData_.icType_;

  return icData_.opData_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setupIC_or_NODESET
// Purpose       :
// Special Notes : Assumes that allNodes_ is set up.
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool
setupIC_or_NODESET(
  Parallel::Machine                     comm,
  const NodeNamePairMap &               all_nodes,
  N_LAS_Vector &                        nextSolnVec,
  N_LAS_Vector &                        flagSolnVec,
  int                                   icType,
  std::vector<Util::OptionBlock> &      initBlockVec,
  NodeNamePairMap &                     op_data,
  int &                                 op_found_,
  int &                                 total_soln_)
{
  int lid = 0;
  bool success = false;

  int totalNodes = all_nodes.size();

  nextSolnVec.putScalar(0.0);
  flagSolnVec.putScalar(-1.0);

  for (NodeNamePairMap::const_iterator it = all_nodes.begin(), nodes_end = all_nodes.end(); it != nodes_end ; ++it)
  {
    flagSolnVec[(*it).second.first] = 0;
  }

  std::set<std::string> notFoundInCkt;
  std::set<std::string> matched;
  for (int icBlockIndex = 0, icBlockEnd = initBlockVec.size(); icBlockIndex < icBlockEnd; ++icBlockIndex)
  {
    // param loop
    for (ParameterList::const_iterator itPar  = initBlockVec[icBlockIndex].getParams().begin(), endPar = initBlockVec[icBlockIndex].getParams().end(); itPar != endPar; ++itPar)
    {
      std::string node("");
      double value(0.0);

      // the first tag is always "V".  At this point I variables are not supported,
      // and there is an error trap for them in the OptionBlock.C file.
      ++itPar;
      node = itPar->tag();
      ++itPar;

      ExtendedString strValue(itPar->tag());
      if (strValue.isValue())
      {
        value = strValue.Value();
      }
      else
      {
        Report::UserFatal0() << "Problems processing " << icType << " values";
      }

      int found = false;

      NodeNamePairMap::const_iterator iterCI = all_nodes.find(node);
      if (iterCI != all_nodes.end())
      {
        lid = iterCI->second.first;
        nextSolnVec[lid] = value;
        flagSolnVec[lid] = 1;
        found = true;
      }

#ifdef Xyce_DEBUG_IC
      if (found)
      {
        Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
        Xyce::dout()
          << "procID="<<getProcID() << "  "
          << icType + " found, and set: V(" << node << "):  solution["<<lid<<"] = " << value << std::endl;
      }
#endif

      Parallel::AllReduce(comm, MPI_SUM, &found, 1);

      if (found)
      {
        success = true;

        if (matched.find(node) == matched.end())
        {
          matched.insert(node);
        }
      }
      else
      {
        notFoundInCkt.insert(node);
      }
    }
  }

  nextSolnVec.importOverlap();
  flagSolnVec.importOverlap();

  op_found_ = matched.size();
  total_soln_ = totalNodes;

  std::set<std::string>::iterator matched_i = matched.begin();
  std::set<std::string>::iterator matched_end = matched.end();
  for (; matched_i != matched_end ; ++matched_i)
  {
    // only add this to the opData map if it is local to the processor.
    // The allNodes object only contains local stuff.  Otherwise there
    // isn't much point.
    NodeNamePairMap::const_iterator iterCI = all_nodes.find(*matched_i);
    //if (all_nodes.find(*matched_i) != all_nodes.end())
    if (iterCI != all_nodes.end())
    {
      op_data[*matched_i].first = iterCI->second.first;
      op_data[*matched_i].second = nextSolnVec[op_data[*matched_i].first];
    }
  }

#ifdef Xyce_DEBUG_IO
  // Identify nodes in the circuit which were not set.
  std::set<std::string> local;
  for (NodeNamePairMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end ; ++it)
  {
    if (matched.find((*it).first) == matched.end())
      local.insert((*it).first);
  }

  Util::Marshal mout;
  mout << local;

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  std::set<std::string> notSpecified;
  if (Parallel::rank(comm) == 0) {
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);

      std::vector<std::string> x;
      min >> x;
      std::copy(x.begin(), x.end(), std::inserter(notSpecified, notSpecified.begin()));
    }
  }
#endif // debug io

  if (Parallel::rank(comm) == 0)
  {

#ifdef Xyce_DEBUG_IO
    if (totalNodes > matched.size())
    {
      Xyce::dout() << icType << ":  Initialized " << matched.size() << " of a possible " << totalNodes << " nodes" << std::endl;
    }
    else
    {
      Xyce::dout() << icType << ":  All " << totalNodes << " nodes initialized" << std::endl;
    }
#endif

    if (notFoundInCkt.size() > 0)
    {
      lout() << icType << ":  Nodes specified on ." << icType << " line, but not present in circuit. (ignoring):" << std::endl;

      for (std::set<std::string>::iterator nm_i = notFoundInCkt.begin(), end = notFoundInCkt.end(); nm_i != end; ++nm_i)
      {
        lout() << *nm_i << std::endl;
      }
    }

#ifdef Xyce_DEBUG_IO
    // Note: for a typical .IC line, this list of unspecified nodes will be
    // a very long list.  So, don't output except in debug mode.
    if (notSpecified.size() > 0)
    {
      dout() << icType << ":  Nodes present in circuit, but not specified on ." << icType << " line(ignoring):" << std::endl;

      for (std::set<std::string>::iterator nm_i = notSpecified.begin(), end = notSpecified.end(); nm_i != end; ++nm_i)
      {
        dout() << *nm_i << std::endl;
      }
    }
#endif // Xyce_DEBUG_IO
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputIC_or_NODESET
// Purpose       : Outputs either an *.ic file, either in .ic or .nodeset format.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
void
outputIC_or_NODESET(
  Parallel::Machine             comm,
  std::ofstream &               os,
  const std::string &           saveFileType_,
  NodeNamePairMap &             all_nodes,
  const N_LAS_Vector &          solution)
{
  for (NodeNamePairMap::iterator it = all_nodes.begin(), name_end = all_nodes.end(); it != name_end ; ++it)
  {
    int index =(*it).second.first;
    (*it).second.second = solution[index];
  }

  Util::Marshal mout;
  mout << all_nodes;

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  if (Parallel::rank(comm) == 0) {
    NodeNamePairMap global;
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);

      NodeNamePairMap x;
      min >> x;
      for (NodeNamePairMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end ; ++it)
      {
        ExtendedString tmpName((*it).first);
        tmpName.toUpper();
        std::string::size_type pos = tmpName.rfind("BRANCH");
        if (pos == std::string::npos)
        {
          os << saveFileType_ << " V(";
          os <<(*it).first << ") = " << (*it).second.second << std::endl;
        }
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::writeDCOPRestart
// Purpose       : Output DC operating point solution vars for use in subsequent
//                 runs to speed up DCOP
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/10/06
//-----------------------------------------------------------------------------
void writeDCOPRestart(Parallel::Machine comm, std::ofstream &os, const NodeNamePairMap &all_nodes)
{
  Util::Marshal mout;
  mout << all_nodes;
  
  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  if (Parallel::rank(comm) == 0) {
    NodeNamePairMap global;
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);

      NodeNamePairMap x;
      min >> x;
      for (NodeNamePairMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end ; ++it)
      {
        os << (*it).first << "   " << (*it).second.second << std::endl;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::readDCOPRestart
// Purpose       : Attempt to input better starting estimate of DC operating
//                 point from a previous run
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/11/06
//-----------------------------------------------------------------------------
bool
readDCOPRestart(
  Parallel::Machine             comm,
  std::ifstream &               is,
  const NodeNamePairMap &       all_nodes,
  N_LAS_Vector &                nextSolnVec,
  N_LAS_Vector &                flagSolnVec,
  NodeNamePairMap &             op_data,
  int &                         op_found_,
  int &                         total_soln_)
{
  bool success = false;
  std::set<std::string> matched;
  std::set<std::string> read;

  flagSolnVec.putScalar(-1.0);
  nextSolnVec.putScalar(0.0);

  // Clear flags
  for (NodeNamePairMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end ; ++it)
  {
    flagSolnVec[(*it).second.first] = 0;
  }

  // Populate nextSolnVec, flagSolnVec, op_data and matched for all nodes mathcing on this processor
  for (;;)
  {
    std::string node;
    double value;
    
    is >> node;
    if (node == "")
      break;
    is >> value;

    // On root proc, keep track of all node names
    if (Parallel::rank(comm) == 0)
      read.insert(node);

    NodeNamePairMap::const_iterator it = all_nodes.find(node);
    if (it != all_nodes.end())
    {
      int ind = (*it).second.first;
      nextSolnVec[ind] = value;
      flagSolnVec[ind] = 1;
      if (matched.find(node) == matched.end())
      {
        matched.insert(node);
        op_data[node].first = (*it).second.first;
        op_data[node].second = nextSolnVec[(*it).second.first];
      }
    }
  }

  // Overlap the results
  nextSolnVec.importOverlap();
  flagSolnVec.importOverlap();

  int all_node_count = all_nodes.size();

  Parallel::AllReduce(comm, MPI_SUM, &all_node_count, 1);

  Util::Marshal mout;
  mout << matched;

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  int missing_count = 0;
  if (Parallel::rank(comm) == 0) {
    std::set<std::string> x;
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);
      min >> x;
      std::copy(x.begin(), x.end(), std::inserter(matched, matched.begin()));
    }

    std::set<std::string> missing;
    std::set_difference(read.begin(), read.end(), matched.begin(), matched.end(), std::inserter(missing, missing.begin()));

    if (read.size() > matched.size())
      Xyce::dout() << "DCOP restart:  Initialized " << matched.size() << " of a possible " << read.size() << " nodes" << std::endl;
    else
      Xyce::dout() << "DCOP restart:  All " << all_node_count << " nodes initialized" << std::endl;

    missing_count = missing.size();
    if (missing.size() > 0)
    {
      Xyce::dout() << "DCOP restart:  Nodes specified but not present in circuit:" << std::endl;
      for (std::set<std::string>::const_iterator it = missing.begin(), end = missing.end(); it != end; ++it)
        Xyce::dout() << (*it) << std::endl;
    }

    if (all_node_count - read.size() > 0)
      Xyce::dout() << "DCOP restart:  " << (all_node_count - read.size()) << " nodes present in circuit but not specified:" << std::endl;
  }

  Parallel::Broadcast(comm, &missing_count, 1, 0);

  if (missing_count == 0)
    success = true;

  return success;
}

} // namespace IO
} // namespace Xyce

