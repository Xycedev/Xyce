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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_fwd.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_fwd_h
#define Xyce_N_NLS_fwd_h

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Enum          : AnalysisMode
// Description   : These anlysis modes influence the choice of parameters for
//                 the nonlinear solver. The mode is set by the
//                 setAnalysisMode() function in the Manager class.
//-----------------------------------------------------------------------------
enum AnalysisMode
{
  DC_OP,
  DC_SWEEP,
  TRANSIENT,
  HB_MODE,
  NUM_MODES
};

class Manager;
class ParamMgr;
class ReturnCodes;
class NonLinInfo;
class NLParams;
class NonLinearSolver;
class ConductanceExtractor;
class Sensitivity;
class TwoLevelNewton;
class ConstraintBT;
class MatrixFreeEpetraOperator;
class DampedNewton;

} // namespace Nonlinear
} // namespace Xyce

namespace LOCA {
  class Stepper;
  namespace MultiContinuation {
    class AbstractGroup;
  }
  class GlobalData;
  namespace StatusTest {
    class Wrapper;
  }
}

namespace N_NLS_NOX {
  class SharedSystem;
  class Group;
  class Vector;
}

namespace N_NLS_LOCA {
  class Group;
}

typedef Xyce::Nonlinear::DampedNewton N_NLS_DampedNewton;
typedef Xyce::Nonlinear::MatrixFreeEpetraOperator N_NLS_MatrixFreeEpetraOperator;
typedef Xyce::Nonlinear::ConstraintBT N_NLS_ConstraintBT;
typedef Xyce::Nonlinear::ParamMgr N_NLS_ParamMgr;
typedef Xyce::Nonlinear::TwoLevelNewton N_NLS_TwoLevelNewton;
typedef Xyce::Nonlinear::ConductanceExtractor N_NLS_ConductanceExtractor;
typedef Xyce::Nonlinear::Sensitivity N_NLS_Sensitivity;
typedef Xyce::Nonlinear::NonLinearSolver N_NLS_NonLinearSolver;
typedef Xyce::Nonlinear::NLParams N_NLS_NLParams;
typedef Xyce::Nonlinear::NonLinInfo N_NLS_NonLinInfo;
typedef Xyce::Nonlinear::ReturnCodes N_NLS_ReturnCodes;
typedef Xyce::Nonlinear::Manager N_NLS_Manager;

#endif // Xyce_N_UTL_fwd_h
