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
// Filename       : $RCSfile: N_IO_CircuitBlock.h,v $
//
// Purpose        : Declare the circuit level containers for holding netlist
//                  circuit data and the associated circuit level methods.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/06/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.88 $
//
// Revision Date  : $Date: 2014/08/07 23:17:58 $
//
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_CircuitBlock_h
#define Xyce_N_IO_CircuitBlock_h

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>

#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_DeviceBlock.h>
#include <N_UTL_OptionBlock.h>
#include <N_TOP_InsertionTool.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {

typedef std::pair<std::ifstream*,SpiceSeparatedFieldTool*> FileSSFPair;

//-----------------------------------------------------------------------------
// Class         : CircuitBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/06/2001
//-----------------------------------------------------------------------------

class CircuitBlock
{
public:
  // Constructor.
  CircuitBlock(
     const std::string &                                                fileName,
     CmdParse &                                                         cp,
     CircuitMetadata &                                                  md,
     std::map<std::string,int> &                                        mn,
     std::map<std::string,FileSSFPair> &                                ssfm,
     CircuitContext &                                                   cc,
     int &                                                              uc,
     std::map<std::string, RCP<Device::InstanceBlock> > &               dNames,
     std::set<std::string> &                                            nNames,
     AliasNodeMap &                                                     alias_node_map,
     const std::vector< std::pair< std::string, std::string> > &        externalNetlistParams);

  // Constructor.
  CircuitBlock(
     const std::string &                                                fileName,
     const std::vector<SpiceSeparatedFieldTool::StringToken> &          parsedInputLine,
     CmdParse &                                                         cp,
     CircuitMetadata &                                                  md,
     std::map<std::string,int> &                                        mn,
     std::map<std::string,FileSSFPair> &                                ssfm,
     CircuitContext &                                                   cc,
     int &                                                              uc,
     CircuitBlock *                                                     mainCircPtr,
     CircuitBlock *                                                     parentCircPtr,
     DistributionTool*                                                  dtPtr,
     Topo::InsertionTool*                                               itPtr,
     Device::DeviceInterface*                                           diPtr,
     std::map<std::string, RCP<Device::InstanceBlock> > &               dNames,
     std::set<std::string> &                                            nNames,
     AliasNodeMap &                                                     alias_node_map,
     const std::vector< std::pair< std::string, std::string> > &        externalNetlistParams,
     bool                                                               removeCvar,
     bool                                                               removeDvar,
     bool                                                               removeIvar,
     bool                                                               removeLvar,
     bool                                                               removeMvar,
     bool                                                               removeQvar,
     bool                                                               removeRvar,
     bool                                                               removeVvar,
     bool                                                               replgndvar);

  // Destructor.
  ~CircuitBlock();

private:
  // Copy Constructor.
  CircuitBlock(CircuitBlock const& rhsCB);
  CircuitBlock &operator=(const CircuitBlock& rhsCB);

public:

  // Get the name of the subcircuit, if this is a subcircuit  
  const std::string& getName() const {
    return name_;
  }

  // Get the name of the circuit, if this is the root
  const std::string &getTitle() const {
    return title_;
  }

  const AliasNodeMap &getAliasNodeMap() const {
    return aliasNodeMap_;
  }

  CircuitContext* getCircuitContextPtr() {
    return &circuitContext_;
  }

  // This function parses the netlist file and fills in the
  // details of the circuit. This is phase 1 of netlist parsing.
  // The devices cannot be completely handled in this phase.
  bool parseNetlistFilePass1();
  bool parseNetlistFilePass1(const std::string &libSelect, std::string libInside);

  // Perform the second pass over the netlist, this phase primarily
  // handles devices.
  bool parseNetlistFilePass2();

  // Perform special pass for mutual inductances
  bool parseMutualInductances();

  // Print the contents of CircuitBlock.
  void print();

  // Set data_->ssfPtr_ .
  void setSSFPtr( SpiceSeparatedFieldTool* ssfPtr );

  void setStartPosition();
  void setEndPosition();
  void setFilePosition(std::streampos const& position);
  void setLinePosition(int const& position);
  const std::streampos getStartPosition() const;
  const std::streampos getEndPosition() const;
  int getLineStartPosition() const;
  int getLineEndPosition() const;

  // Extract subcircuit data from parsedLine.
  bool extractSubcircuitData(std::string, int);

  // Instatiate all of the devices in the current (sub)circuit. This
  // method will be invoked during as a part of pass 2 operations. It
  // will be invoked to create instances of the devices in the main
  // circuit and in each subcircuit instance.
  bool instantiateDevices(std::string libSelect, std::string libInside);

  int getDeviceNames (std::vector<std::string> &names);
  void checkDeviceNames (const std::vector<std::string> &names);

  // Add a device to the circuit.
  void addTableData(DeviceBlock & device );

  // Add a mutual inductor to the circuit.
  void addMutualInductor( DeviceBlock& device, CircuitContext* context );

  // Add a model to the circuit.
  void addModel(ParameterBlock & model, std::string const& modelPrefix);

  // Add a set of options corresponding to a .OPTIONS netlist line
  // to the circuit.
  void addOptions(const OptionBlock &options );

  // Search the subcircuitInstanceTable of the current circuit block for the
  // subcircuit of the given name. If it is not found, recursively
  // search each parent subcircuit. Return a pointer to the circuit
  // block if it is found, otherwise return NULL.
  CircuitBlock* findSubcircuit( std::string const& subcircuitName );

  void registerDistributionTool(DistributionTool* dtPtr);
  void registerInsertionTool(Xyce::Topo::InsertionTool* insertionToolPtr);
  void registerDeviceInterface(Device::DeviceInterface* devIntPtr);

  // Receive the circuit context (from the distribution tool).
  bool receiveCircuitContext(CircuitContext&  ccIn);

  void registerGlobalParams() const;

  // Change netlist file name in slave procs
  void setFileName ( std::string & );

  // Process a device line on processor zero, or serial.
  bool handleDeviceLine(
     std::vector<SpiceSeparatedFieldTool::StringToken> const& deviceLine,
     const std::string &libSelect="", const std::string &libInside="");

  // resolve expressions in optionBlocks like .print
  bool resolveExpressionsInOptionBlocks();

  // write out a copy of the netlist if it is requested.
  void writeOutNetlist();

public:

  // Lookup table for initial conditions
  std::map< std::string, std::vector< SpiceSeparatedFieldTool::StringToken > >
  initCondIndex;

private:
  std::string netlistFileName_;
  std::string title_;    // For top level circuit, given by first line of netlist.
  std::string name_;     // For subcircuits.
  std::string analysisName_;

  std::set<std::string> & nodeNames_;
  std::map<std::string,int> & modelNames_;
  ModelMap modMap_;
  std::map<std::string, RCP<Device::InstanceBlock> > & deviceNames_;
  std::list<Util::OptionBlock> optionsTable_;
  Util::OptionBlock deviceOptions_;

  int & useCount_;  // Counts copies so deletion of circuitBlockTable_ 
                    // can be done properly in the class destructor.
  std::map<std::string, CircuitBlock*> circuitBlockTable_;

  // keep track of K lines that need extracted
  std::multimap< CircuitContext *, DeviceBlock > rawMIs_;

  CmdParse & commandLine_;
  
  // Circuit Context object to hold context information.
  CircuitContext & circuitContext_;

  CircuitMetadata & metadata_;

  std::vector< std::pair< std::string, std::string> > externalNetlistParams_;

  // Checking the netlist syntax.
  bool netlistSave_;

  int devProcessedNumber_;

  std::vector<std::string> nodeList_;         // External nodes of a subcircuit.
  std::vector<DeviceBlock> deviceList_;

  // This is a map of node aliases.  Interface nodes to a subcircuit are removed as
  // the subcircuit is expanded into the netlist. We'll store the names of the interface
  // nodes in case the user accesses them elsewhere (as in a print statement)
  // The keys are the alias names and the values are the real cicuit node names
  AliasNodeMap &                aliasNodeMap_;
  std::set< std::string >       aliasNodeMapHelper_;

  std::ifstream netlistIn_;

  SpiceSeparatedFieldTool* ssfPtr_;

  std::streampos fileStartPosition_;
  std::streampos fileEndPosition_;
  int lineStartPosition_;
  int lineEndPosition_;

  // Top level circuit pointer
  CircuitBlock* mainCircuitPtr_;

  // For subcircuits, points to the circuitBlock instance that contains
  // this subcircuit. NULL for top level.
  CircuitBlock* parentCircuitPtr_; 

  DistributionTool* distToolPtr_;

  std::map<std::string,FileSSFPair> & ssfMap_;

  DeviceBlock device_;

  // These are added to check allow removal of "redundant" devices (where
  // all device nodes are the same).
  bool remove_redundant_C_;     //capacitors
  bool remove_redundant_D_;     //diodes
  bool remove_redundant_I_;     //independent current sources
  bool remove_redundant_L_;     //inductors
  bool remove_redundant_M_;     //MOSFETS
  bool remove_redundant_Q_;     //BJTs
  bool remove_redundant_R_;     //resistors
  bool remove_redundant_V_;     //independent voltage sources

  // This is added to turn on an option where we replace any occurrence of
  // "GND", "GND!", or "GROUND" with "0" (so that all four terms are
  // synonyms)
  bool replace_ground_;

  std::vector<SpiceSeparatedFieldTool::StringToken> parsedLine_;

  ParameterBlock tmpModel_;

  Xyce::Topo::InsertionTool* insertionToolPtr_;
  Device::DeviceInterface* devIntPtr_;

  //This function preprocesses the netlist file to provide the user the
  //option of removing "redundant" devices (devices where all the nodes are
  //the same.  The info gathered here affects the phase 1 and phase 2 parse.
  bool parsePreprocess(const std::string & netlistFileName);

  //This function will reproduce a copy of the netlist file under the name
  //netlistfilename_copy.cir (to be used, in the end, to produce netlist files
  //which contain large resistors connecting "dangling" nodes to ground.
  void produceUnflattenedNetlist(const std::string & netlistFileName);

  // Handle a netlist line, determine the line type and take the
  // appropriate action.
  bool handleLinePass1( bool & result, std::map<std::string,int> & par,
                        std::map<std::string,int> & fun, ModelMap & modMap,
                        std::map<std::string,int> & sub, 
                        const std::string &libSelect, std::string &libInside );

  bool getLinePass2(std::vector<SpiceSeparatedFieldTool::StringToken>& line,
                    const std::string &libSelect, std::string &libInside);

  bool removeTwoTerminalDevice(const char linetype,
                               const ExtendedString & node1,
                               const ExtendedString & node2);

  bool removeThreeTerminalDevice(const char linetype,
                                 const ExtendedString & node1,
                                 const ExtendedString & node2,
                                 const ExtendedString & node3);

  // Handle a line for Mutual Inductance Pass
  bool getLinePassMI();

  // Handle a netlist .include or .lib line, return the include file, and lib strings.
  void handleIncludeLine(
      std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
      const ExtendedString &, std::string& includefile, std::string &libSelect, std::string &libInside);

  // Handle a netlist .endl line.
  void handleEndlLine(
      std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
      const std::string &libSelect,
      std::string &libInside);

  // Post process analysis commands.
  bool handleAnalysis();

  // Post process a mutual inductor in the current circuit.
  bool handleMutualInductance( DeviceBlock & device );

  // Parse the given include file adding the contents
  // to the current CircuitBlock.
  bool parseIncludeFile(std::string const& includeFile, std::string const& libSelect,
       std::map<std::string,int> & par, std::map<std::string,int> & fun,
       ModelMap & modMap, std::map<std::string,int> & sub);

  // Parse the given include file for 2nd pass
  bool parseIncludeFile2(std::string const& includeFiles,
          const std::string &libSelect);

  // Retrieve separate IC= data from line or external file and temporarily
  // store in CircuitBlock
  void handleInitCond(
   std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine );

  // Fully parse and instantiate a single device.
  bool instantiateDevice( DeviceBlock & device,
      std::string & prefix, const std::map<std::string,std::string>& nodeMap,
      const std::string &libSelect, const std::string &libInside);

  // Expand a subcircuit instance by adding the devices and
  // device models that compose the subcircuit to the main
  // (top level) circuit. Prepend device names and nodes with
  // subcircuitPrefix.
  bool expandSubcircuitInstance(DeviceBlock & subcircuitInstance,
          const std::string &libSelect, const std::string &libInside);

  //
  bool getFileLocation( std::string const& path, std::string const& file,
                        std::string& location);
  bool getFileLocation( std::string const& path, std::string const& file,
                        char separator, std::string& location);

  // read a line of input from an istream
  void readline( std::istream & in, char * line );

};

// check for name collisions between nodes and devices
void checkNodeDevConflicts(CircuitBlock* cktBlk, N_PDS_Comm* pdsComm);

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitBlock N_IO_CircuitBlock;

#endif // Xyce_N_IO_CircuitBlock_h
