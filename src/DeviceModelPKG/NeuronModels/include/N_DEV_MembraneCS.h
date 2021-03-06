//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_MembraneCS.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 08/11/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/03/19 17:23:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MembraneCS_h
#define Xyce_N_DEV_MembraneCS_h

#include <N_DEV_MembraneModel.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MembraneCS
// Purpose       : This is class defines a membrane with the Connor Stevens equations
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembraneCS : public MembraneModel
{
public:
  MembraneCS(const SolverState & ss1);
  ~MembraneCS() {}

  void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp );
  void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea);
  void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea);
  void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea);
  void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea);
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembraneCS N_DEV_MembraneCS;

#endif
