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
// Filename       : $RCSfile: N_DEV_ExternCodeInterface.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/09/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_ExternCodeInterface.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ExternCodeInterface::ExternCodeInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
ExternCodeInterface::ExternCodeInterface()
{}

//-----------------------------------------------------------------------------
// Function      : ExternCodeInterface::ExternCodeInterface
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
ExternCodeInterface::ExternCodeInterface
  (const ExternCodeInterface &right)
{}

//-----------------------------------------------------------------------------
// Function      : ExternCodeInterface::~ExternCodeInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
ExternCodeInterface::~ExternCodeInterface()

{}

} // namespace Device
} // namespace Xyce
