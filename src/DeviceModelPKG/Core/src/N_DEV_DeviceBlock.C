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
// Filename       : $RCSfile: N_DEV_DeviceBlock.C,v $
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
// Revision Number: $Revision: 1.38 $
//
// Revision Date  : $Date: 2014/05/19 20:00:59 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <iostream>

#include <N_DEV_DeviceBlock.h>
#include <N_ERH_ErrorMgr.h>

#include <N_PDS_Comm.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const std::string & name, const std::string &type, int level)
  : name_(name),
    level_(level),
    type_(type),
    netlistLocation_()
{}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const ModelBlock &right)
  : name_(right.name_),
    type_(right.type_),
    level_(right.level_),
    netlistLocation_(right.netlistLocation_),
    params(right.params)
{}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::~ModelBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::~ModelBlock()
{
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock & ModelBlock::operator=(const ModelBlock &right)
{
  name_   = right.name_;
  type_   = right.type_;
  level_ =  right.level_;

  netlistLocation_ = right.netlistLocation_;

  params = right.params;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream & os, const ModelBlock & mb)
{
  std::vector<Param>::const_iterator it_pL, end_pL;

  os << "Model Block" << std::endl;
  os << "Model:  " << mb.name_ << std::endl;
  os << " type:  " << mb.type_ << std::endl;
  os << " Level: " << mb.level_ << std::endl;
  os << " netlistLocation_: " << mb.netlistLocation_ << std::endl;
  os << " Tagged Params" << std::endl;
  os << " -------------" << std::endl;

  it_pL=mb.params.begin();
  end_pL=mb.params.end();
  for ( ; it_pL != end_pL; ++it_pL)
  {
    os << it_pL->tag() << "\t" << it_pL->stringValue() << std::endl;
  }

  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void ModelBlock::clear()
{
  name_ = "";
  type_ = "";
  level_ = 0;
  params.clear();

  netlistLocation_ = NetlistLocation();
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * ModelBlock::instance() const
{
  return new ModelBlock();
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int ModelBlock::packedByteCount() const
{

  int byteCount = 0;

  int size, length, i;

  std::vector<Param>::const_iterator it_tpL;

  //----- count name
  length = name_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count type
  length = type_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count level
  byteCount += sizeof(int);

  //----- count params
  size = params.size();
  byteCount += sizeof(int);
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += it_tpL->packedByteCount();
  }

  //----- count netlistFileName_
  length = netlistLocation_.getPath().length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count lineNumber_
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::pack
// Purpose       : Packs ModelBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void ModelBlock::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{

  int size, length;
  int i;
  std::vector<Param>::const_iterator it_tpL;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  length = name_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( name_.c_str(), length, buf, bsize, pos );

  //----- pack type
  length = type_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( type_.c_str(), length, buf, bsize, pos );

  //----- pack level
  comm->pack(&level_, 1, buf, bsize, pos );

  //----- pack params
  size = params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    it_tpL->pack( buf, bsize, pos, comm );
  }

  //----- pack netlistFileName_
  length = netlistLocation_.getPath().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( netlistLocation_.getPath().c_str(), length, buf, bsize, pos );

  //----- packlineNumber_
  int line_number = netlistLocation_.getLineNumber();
  comm->pack(&line_number, 1, buf, bsize, pos );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packed " << pos << " bytes for ModelBlock: " <<
    name_ << std::endl;
#endif

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "ModelBlock::pack") << "Predicted pos does not match actual pos";
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::unpack
// Purpose       : Unpacks ModelBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void ModelBlock::unpack(char * pB, int bsize,int & pos, N_PDS_Comm * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack type
  comm->unpack( pB, bsize, pos, &length, 1 );
  type_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack level
  comm->unpack( pB, bsize, pos, &level_, 1 );

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  params.clear();
  Param dp;
  for( i = 0; i < size; ++i )
  {
    dp.unpack( pB, bsize, pos, comm );
    params.push_back( dp );
  }

  //----- unpack netlistFileName_
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistLocation_.setPath(std::string( (pB+pos), length));
  pos += length;

  //----- unpack lineNumber_
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  netlistLocation_.setLineNumber(line_number);

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacked " << pos << " bytes for ModelBlock: " <<
    name_ << std::endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const std::string &name)
  : name_(name),
    modelName_(),
    iNumNodes(0),
    numIntVars(0),
    numExtVars(0),
    numStateVars(0),
    modelFlag(0),
    sourceFlag(0),
    bsourceFlag(0),
    offFlag(0),
    off(0),
    netlistLocation_()
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const InstanceBlock &right)
  : name_      (right.name_),
    modelName_(right.modelName_),
    iNumNodes (right.iNumNodes),
    numIntVars(right.numIntVars),
    numExtVars(right.numExtVars),
    numStateVars(right.numStateVars),
    modelFlag (right.modelFlag),
    sourceFlag(right.sourceFlag),
    bsourceFlag(right.bsourceFlag),
    offFlag   (right.offFlag),
    off       (right.off),
    netlistLocation_(right.netlistLocation_),
    params    (right.params)
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::~InstanceBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::~InstanceBlock ()
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock & InstanceBlock::operator=(InstanceBlock &right)
{
  name_      = right.name_;
  modelName_ = right.modelName_;
  iNumNodes = right.iNumNodes;
  numIntVars= right.numIntVars;
  numExtVars= right.numExtVars;
  numStateVars= right.numStateVars;
  modelFlag = right.modelFlag;
  sourceFlag= right.sourceFlag;
  bsourceFlag= right.bsourceFlag;
  offFlag   = right.offFlag;
  off       = right.off;
  netlistLocation_  = right.netlistLocation_;
  params    = right.params;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void InstanceBlock::clear ()
{
  name_ = InstanceName();
  modelName_  = "";
  iNumNodes  = 0;
  numIntVars = 0;
  numExtVars = 0;
  numStateVars = 0;
  modelFlag  = 0;
  sourceFlag = 0;
  bsourceFlag = 0;
  offFlag    = 0;
  off        = 0;
  netlistLocation_ = NetlistLocation();

  params.clear();
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream & os, const InstanceBlock & ib)
{
  std::vector<Param>::const_iterator it_tpL, end_tpL;

  os << "Instance Block" << std::endl;
  os << "Name:    " << ib.name_ << std::endl;
  os << " Model:  " << ib.getModelName() << std::endl;
  os << " # Nodes: " << ib.iNumNodes << std::endl;
  os << " # Int Vars: " << ib.numIntVars << std::endl;
  os << " # Ext Vars: " << ib.numExtVars << std::endl;
  os << " # State Vars: " << ib.numStateVars << std::endl;
  os << " modelFlag: " << ib.modelFlag << std::endl;
  os << " sourceFlag: " << ib.sourceFlag << std::endl;
  os << " bsourceFlag: " << ib.bsourceFlag << std::endl;
  os << " offFlag: " << ib.offFlag << std::endl;
  os << " off: " << ib.off << std::endl;
  os << " netlistLocation_: " << ib.netlistLocation_ << std::endl;
  os << " Tagged Params" << std::endl;
  os << " -------------" << std::endl;

  it_tpL=ib.params.begin();
  end_tpL=ib.params.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    os << it_tpL->tag() << "\t" << it_tpL->stringValue() << std::endl;
  }

  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * InstanceBlock::instance() const
{
  return new InstanceBlock();
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::packedByteCount
// Purpose       : count bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int InstanceBlock::packedByteCount() const
{

  int byteCount = 0;

  int size, length, i;
  std::vector<Param>::const_iterator it_tpL;

  //----- count name
  const std::string &name = getInstanceName().getEncodedName();
  length = name.length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count getModelName()
  length = getModelName().length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count params
  size = params.size();
  byteCount += sizeof(int);
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += it_tpL->packedByteCount();
  }

  //----- count iNumNodes
  byteCount += sizeof(int);

  //----- count numIntVars
  byteCount += sizeof(int);

  //----- countnumExtVars
  byteCount += sizeof(int);

  //----- count numStateVars
  byteCount += sizeof(int);

  //----- count modelFlag
  byteCount += sizeof(int);

  //----- count sourceFlag
  byteCount += sizeof(int);

  //----- count bsourceFlag
  byteCount += sizeof(int);

  //----- count offFlag
  byteCount += sizeof(int);

  //----- pack off
  byteCount += sizeof(int);

  //----- count netlistFileName_
  length = netlistLocation_.getPath().length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count lineNumber_
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::pack
// Purpose       : Packs InstanceBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void InstanceBlock::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{

  int size, length;
  int i;
  std::vector<Param>::const_iterator it_tpL;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  const std::string &name = getInstanceName().getEncodedName();
  length = name.length();
  comm->pack(&length, 1, buf, bsize, pos);
  comm->pack(name.c_str(), length, buf, bsize, pos);

  //----- pack getModelName()
  length = getModelName().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( getModelName().c_str(), length, buf, bsize, pos );

  //----- pack params
  size = params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    it_tpL->pack( buf, bsize, pos, comm );
  }

  //----- pack iNumNodes
  comm->pack(&iNumNodes, 1, buf, bsize, pos );

  //----- pack numIntVars
  comm->pack(&numIntVars, 1, buf, bsize, pos );

  //----- pack numExtVars
  comm->pack(&numExtVars, 1, buf, bsize, pos );

  //----- pack numStateVars
  comm->pack(&numStateVars, 1, buf, bsize, pos );

  //----- pack modelFlag
  i = modelFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack sourceFlag
  i = sourceFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack bsourceFlag
  i = bsourceFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack offFlag
  i = offFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack off
  i = off;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack name
  length = netlistLocation_.getPath().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( netlistLocation_.getPath().c_str(), length, buf, bsize, pos );

  //----- pack lineNumber_
  int line_number = netlistLocation_.getLineNumber();
  comm->pack(&line_number, 1, buf, bsize, pos );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packed " << pos << " bytes for InstanceBlock: " <<
    getInstanceName() << std::endl;
#endif
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "InstanceBlock::pack") << "Predicted pos does not match actual pos";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::unpack
// Purpose       : Unpacks InstanceBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void InstanceBlock::unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = InstanceName(std::string( (pB+pos), length));
  pos += length;

  //----- unpack getModelName()
  comm->unpack( pB, bsize, pos, &length, 1 );
  modelName_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  params.clear();
  Param dp;
  for( i = 0; i < size; ++i )
  {
    dp.unpack( pB, bsize, pos, comm );
    params.push_back( dp );
  }

  //----- unpack iNumNodes
  comm->unpack( pB, bsize, pos, &iNumNodes, 1 );

  //----- unpack numIntVars
  comm->unpack( pB, bsize, pos, &numIntVars, 1 );

  //----- unpack numExtVars
  comm->unpack( pB, bsize, pos, &numExtVars, 1 );

  //----- unpack numStateVars
  comm->unpack( pB, bsize, pos, &numStateVars, 1 );

  //----- unpack modelFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  modelFlag = ( i != 0 );

  //----- unpack sourceFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  sourceFlag = ( i != 0 );

  //----- unpack bsourceFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  bsourceFlag = ( i != 0 );

  //----- unpack offFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  offFlag = ( i != 0 );

  //----- unpack off
  comm->unpack( pB, bsize, pos, &i, 1 );
  off = ( i != 0 );

  //----- unpack netlistFileName_
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistLocation_.setPath(std::string( (pB+pos), length));
  pos += length;

  //----- unpack lineNumber_
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  netlistLocation_.setLineNumber(line_number);

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacked " << pos << " bytes for InstanceBlock: " <<
    getInstanceName() << std::endl;
#endif

}

} // namespace Device
} // namespace Xyce
