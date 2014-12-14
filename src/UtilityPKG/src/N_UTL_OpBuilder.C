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
// Filename       : $RCSfile: N_UTL_OpBuilder.C,v $
//
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Util {
namespace Op {

//-----------------------------------------------------------------------------
// Function      : Builder::createOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:56:45 2014
//-----------------------------------------------------------------------------
Operator *
Builder::createOp(
  ParameterList::const_iterator &       it) const
{
  ParameterList::const_iterator it2 = it;

  Operator *new_op = makeOp(it2);
  if (new_op)
    it = it2;

  return new_op;
}

//-----------------------------------------------------------------------------
// Function      : BuilderManager::createOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:00:36 2014
//-----------------------------------------------------------------------------
Operator *
BuilderManager::createOp(
  ParameterList::const_iterator &       it) const
{
  for (BuilderVector::const_iterator it1 = opBuilderVector_.begin(), end1 = opBuilderVector_.end(); it1 != end1; ++it1) {
    Operator *new_op = (*it1)->createOp(it);
    if (new_op)
      return new_op;
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : BuilderManager::findCreateFunction
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:01:47 2014
//-----------------------------------------------------------------------------
CreateFunction
BuilderManager::findCreateFunction(
  Identifier    id) const
{
  CreateMap::const_iterator it = opCreateMap_.find(id);
  if (it != opCreateMap_.end())
    return (*it).second;

  return 0;
}

} // namespace Op
} // namespace Util
} // namespace Xyce
