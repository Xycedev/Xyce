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
// Filename       : $RCSfile: N_DEV_DeviceBlock.h,v $
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
// Revision Number: $Revision: 1.33 $
//
// Revision Date  : $Date: 2014/05/19 20:00:59 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceBlock_h
#define Xyce_N_DEV_DeviceBlock_h

#include <string>
#include <vector>
#include <iosfwd>

#include <N_DEV_Param.h>
#include <N_DEV_InstanceName.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Packable.h>


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ModelBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
///
/// ModelBlock represents a .MODEL line from the netlist.
///
class ModelBlock : public Packable
{
  friend std::ostream& operator<<(std::ostream& os, const ModelBlock & mb);

public:
  ModelBlock(const std::string &name = "", const std::string &type = "", int level = 1);

  ModelBlock(const ModelBlock & right);
  ModelBlock &operator=(const ModelBlock & right);

  ~ModelBlock();

  const ModelName &getName() const
  {
    return name_;
  }

  void setName(const ModelName &name)
  {
    name_ = name;
  }

  const std::string &getType() const
  {
    return type_;
  }

  void setType(const std::string &type)
  {
    type_ = type;
  }

  int getLevel() const
  {
    return level_;
  }

  void setLevel(int level)
  {
    level_ = level;
  }

  const NetlistLocation &getNetlistLocation() const
  {
    return netlistLocation_;
  }

  void setNetlistLocation(const NetlistLocation &netlist_location)
  {
    netlistLocation_ = netlist_location;
  }

  bool operator==(const ModelBlock &right) const
  {
    return equal_nocase(name_, right.name_);
  }

  bool operator!=(const ModelBlock &right) const
  {
    return !equal_nocase(name_, right.name_);
  }

  void clear();

  // Packing Utils
  Packable * instance() const;
  int packedByteCount() const;

  void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const;
  void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm );

private:
  ModelName             name_;                  ///< Model name
  std::string           type_;                  ///< Model type
  int                   level_;                 ///< Device level
  NetlistLocation       netlistLocation_;       ///< Path and line number of .MODEL command

public:
  std::vector<Param>    params;                 ///< Parameters from the line
};

//-----------------------------------------------------------------------------
// Class         : InstanceBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
///
/// InstanceBlock represent a device instance line from the netlist.
///
class InstanceBlock : public Packable
{
  friend std::ostream& operator<<(std::ostream& os, const InstanceBlock & ib);

public:
  InstanceBlock(const std::string &name = std::string());

  InstanceBlock(const InstanceBlock & right);
  InstanceBlock & operator=(InstanceBlock & right);

  ~InstanceBlock();

  const InstanceName &getInstanceName() const 
  {
    return name_;
  }

  void setInstanceName(const InstanceName &name) 
  {
    name_ = name;
  }

  const ModelName &getModelName() const 
  {
    return modelName_;
  }

  void setModelName(const ModelName &modelName) 
  {
    modelName_ = modelName;
  }

  const NetlistLocation &getNetlistLocation() const {
    return netlistLocation_;
  }

  void setNetlistLocation(const NetlistLocation &netlist_location)
  {
    netlistLocation_ = netlist_location;
  }

  bool operator==(const InstanceBlock &right) const
  {
    return equal_nocase(name_.getEncodedName(), right.name_.getEncodedName());
  }

  bool operator!=(const InstanceBlock &right) const
  {
    return !equal_nocase(name_.getEncodedName(), right.name_.getEncodedName());
  }

  void clear ();

  //Packing Utils
  Packable * instance() const;
  int packedByteCount() const;

  void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const;
  void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm );

private:
  InstanceName          name_;                  ///< Device instance name
  ModelName             modelName_;             ///< Model name if provided
  NetlistLocation       netlistLocation_;       ///< Path and line number of .MODEL command

public:
  std::vector<Param> params;

  int iNumNodes;
  int numIntVars;
  int numExtVars;
  int numStateVars;

  bool modelFlag;
  bool sourceFlag;
  bool bsourceFlag;
  bool offFlag;
  bool off;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::InstanceBlock N_DEV_InstanceBlock;
typedef Xyce::Device::ModelBlock N_DEV_ModelBlock;

#endif
