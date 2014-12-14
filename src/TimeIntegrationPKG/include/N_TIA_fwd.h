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
// Filename       : $RCSfile: N_TIA_fwd.h,v $
//
// Purpose        : AC analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Ting Mei   
//
// Creation Date  : 01/11
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
// Revision Date  : $Date: 2014/07/16 17:32:41 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_fwd_h
#define Xyce_N_TIA_fwd_h

namespace Xyce {
namespace TimeIntg {

} // namespace TimeIntg
} // namespace Xyce

enum Integration_Method_Mode
{ // Time Integrator Method Mode
  TIAMethod_NONE                        = 0,
  TIAMethod_BACKWARD_EULER              = 1,
  TIAMethod_BACKWARD_DIFFERENTIATION_2  = 2,
  TIAMethod_TRAPEZOIDAL                 = 3,
  TIAMethod_VARIABLE_THETA              = 4,
  TIAMethod_A_CONTRACTIVE_2             = 5,
  TIAMethod_BACKWARD_DIFFERENTIATION_15 = 6,
  TIAMethod_ONESTEP                     = 7,
  TIAMethod_GEAR_12                     = 8
};


class N_TIA_BackwardDifferentiation15;
class N_TIA_DataStore;
class N_TIA_MPDEInterface;
class N_TIA_OneStep;
class N_TIA_StepErrorControl;
class N_TIA_TIAParams;
class N_TIA_TimeIntInfo;
class N_TIA_TimeIntegrationMethods;
class N_TIA_TwoLevelError;
class N_TIA_WorkingIntegrationMethod;

#endif // Xyce_N_TIA_fwd_h
