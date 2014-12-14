//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Class         : N_UTL_ExpressionData
//
// Purpose       : This class manages a single expression.
//
// Special Notes : This class owns a single expression class pointer,
//                 and all the auxilliary data needed to manage its
//                 usage.
//
//                 I felt that this would be a good way to avoid
//                 bloat in functions like outputPRINT.
//
//                 I made all the data public.  Right now, the
//                 evaluate function just returns the value of the
//                 expression.
//
//                 The derivative values are all in
//                 the valDerivs vector after each evaluation.  If
//                 it becomes necessary later to access derivatives
//                 for some reason, a developer only has to access
//                 the public valDerivs vector, right after an
//                 evaluate function call.
//
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------

#ifndef N_UTL_ExpressionData_H
#define N_UTL_ExpressionData_H

#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Xyce.h>

class N_LAS_Vector;

namespace Xyce {
namespace Util {

class ExpressionData
{
public:
  ExpressionData (const std::string &expression);
  ExpressionData (const Expression &expression);

  ~ExpressionData();

  const std::string &getExpression() const {
    return expressionString_;
  }

  int setup(Parallel::Machine comm, const IO::OutputMgr &output_manager);

  double evaluate(Parallel::Machine comm, double current_circuit_time,
                  const N_LAS_Vector *solnVecPtr, const N_LAS_Vector *stateVecPtr, const N_LAS_Vector * stoVecPtr, const N_LAS_Vector *solnVecImagPtr = 0) const;

private:
  Expression *                  expression_;                    ///< Compiled Expression
  std::string                   expressionString_;              ///< Original expression string
  int                           unresolvedSymbolCount_;         ///< Count of unresolved symbols in expression
  Op::OpList                    expressionOps_;                 ///< Ops to compute variables
  mutable std::vector<double>   variableValues_;                ///< Cache of computed variables
};

} // namespace Util
} // namespace Xyce

typedef Xyce::Util::ExpressionData N_UTL_ExpressionData;

#endif  // N_UTL_ExpressionData_H
