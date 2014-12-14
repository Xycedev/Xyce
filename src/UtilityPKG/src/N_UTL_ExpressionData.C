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
// Filename       : $RCSfile: N_UTL_ExpressionData.C,v $
//
// Purpose        : Handle data for one expression object
//
// Special Notes  :
//
// Creator        : Richard Schiek
//
// Creation Date  : 8/24/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_Op.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Vector.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Util {

static const int NOT_SETUP = std::numeric_limits<int>::max();

//-----------------------------------------------------------------------------
// Function      : ExpressionData::ExpressionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::ExpressionData (
  const std::string &   expression)
  : expression_(0),
    expressionString_(expression),
    unresolvedSymbolCount_(NOT_SETUP)
{}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::ExpressionData
// Purpose       : constructor -- used to take an existing expression pointer
//                 not a string that can be converted to an expression
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::ExpressionData (
  const Expression &  expression)
  : expression_(new Expression( expression )),
    expressionString_(expression_->get_expression()),
    unresolvedSymbolCount_(NOT_SETUP)
{}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::~ExpressionData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::~ExpressionData ()
{
  delete expression_;

  for (Util::Op::OpList::iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
double
ExpressionData::evaluate(
  Parallel::Machine     comm,
  double                current_circuit_time,
  const N_LAS_Vector *  solnVecPtr,
  const N_LAS_Vector *  stateVecPtr,
  const N_LAS_Vector *  stoVecPtr,
  const N_LAS_Vector *  solnVecImagPtr) const
{
  if (unresolvedSymbolCount_ == NOT_SETUP)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Must call setup() prior to evaluate()";
  }
  else if (unresolvedSymbolCount_ > 0)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Unresolved symbols in expression";
  }

  double value = 0.0;

  if (solnVecPtr)
  {
    // loop over expressionOps_ to get all the values.
    variableValues_.clear();
    for (Util::Op::OpList::const_iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
      variableValues_.push_back(Util::Op::getValue(comm, *(*it), Util::Op::OpData(0, solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr, 0)).real());

    // STD and DDT are implicitly time dependent.  Check the underlying expression
    // for time dependence, and if it is set the current time in the underlying expression
    if (expression_->isTimeDependent())
    {
      expression_->set_sim_time(current_circuit_time);
    }

    // now get expression value and partial derivatives.
    // Note that for these purposes (.PRINT output) we only
    // really need the value, not the derivs.
    if (expression_)
    {
      expression_->evaluateFunction(value, variableValues_);
    }
  }

  return value;
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes : just the declaration, definition at end of file
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {
void convertNodalComputation(std::string &nodalComputation,
                             std::list<Param> &paramList);
} // namespace (unnamed)

//-----------------------------------------------------------------------------
// Function      : ExpressionData::setup
//
// Purpose       : Manages all the basic setup for this class.
//
// Special Notes : Originally, all this stuff was in the constructor,
//                 but it needs to happen after topology is
//                 completely done setting itself up, and the
//                 constructor call was too early for that.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/10/04
//-----------------------------------------------------------------------------
int ExpressionData::setup(Parallel::Machine comm, const IO::OutputMgr &output_manager)
{
  unresolvedSymbolCount_ = 0;

  // allocate expression pointer if we need to
  if( expression_ == NULL )
  {
    expression_ = new Expression(expressionString_);
  }

  std::vector<std::string> varNames;


  // query the expression object for all of its dependent vars.
  expression_->get_names(XEXP_ALL, varNames);

  // this varNames vec is a list of string representations of all of
  // the vars in the expression.  use this list to make a list of
  // Param objects which we can then call the
  // IO::OutputMgr::setParamContextType() and then use
  // getPrintValue_() to get the results.  This reuses code that is
  // already in place to handle solution, state, store and lead
  // currents.

  std::list<Param> param_list;

  for (std::vector<std::string>::iterator it = varNames.begin(), end = varNames.end(); it != end; ++it)
  {
    // based on the type of variable, create the needed Param
    // objects for setParmContextType_ to work.
    int varType = expression_->get_type( *it );
    switch (varType)
    {
      case XEXP_NODAL_COMPUTATION:
        // deconstruct the string and turn it into params, push back into
        // param_list
        convertNodalComputation(*it,param_list);
        break;
      case XEXP_NODE:
        param_list.push_back( Param( "V" , 1 ) );
        param_list.push_back( Param( *it , 0.0 ) );
        break;
      case XEXP_INSTANCE:
      case XEXP_LEAD:
        {
          char leadDesignator = expression_->get_lead_designator( *it);
          std::string currentName("I");
          if( (leadDesignator != 0) && (leadDesignator!=' ') )
          {
            currentName = currentName + leadDesignator;
          }
          param_list.push_back( Param( currentName , 1 ) );
          param_list.push_back( Param( *it , 0.0 ) );
        }
        break;
      case XEXP_SPECIAL:
      case XEXP_STRING:
        param_list.push_back( Param( *it , 0.0 ) );
        break;
      case XEXP_VARIABLE:
        // this case is a global param that must be resolved at each use because
        // it can change during a simulation
        param_list.push_back( Param( "GLOBAL_PARAMETER" , *it ) );
        break;
      default:
        {
          std::string errorMsg("Can't find context for expression variable ");
          errorMsg += *it;
          errorMsg += " in full expression ";
          errorMsg += expressionString_;
          errorMsg += " Please check to ensure this parameter is correct and set in your netlist.";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,errorMsg);
          // increment the number of unresolved strings reported in the
          // return value.
          ++unresolvedSymbolCount_;
          // Just issue the warning for now.  Trying to resolve the
          // string in param_list will cause problems in parallel
          // as not all expressions on all procs will have the same
          // unresolved vars.  The class owner will have to check that
          // at least one of the processors with this expression has
          // all the strings resolved.
          // param_list.push_back( Param( *it , 0.0 ) );
        }
    }
  }

  Parallel::AllReduce(comm, MPI_MIN, &unresolvedSymbolCount_, 1);

  if (unresolvedSymbolCount_ > 0 )
  {
    Report::UserFatal0() << "Can't resolve all symbols in expression " << expressionString_;
  }

  Util::Op::makeOps(comm, output_manager.getOpBuilderManager(), param_list.begin(), param_list.end(), std::back_inserter(expressionOps_));

  return unresolvedSymbolCount_;
}


//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : convertNodalComputation
// Purpose       : given a nodal expression string (e.g. "VM(A,B)"),
//                 construct the set of Params that makeOps would expect for
//                 it
// Special Notes :
//
// Scope         : file-local
// Creator       : Tom Russo
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
void convertNodalComputation(std::string &nodalComputation,
                             std::list<Param> &paramList)
{
  std::list<Param> tempParamList;

  std::size_t firstParen = nodalComputation.find_first_of("(");
  std::size_t lastParen = nodalComputation.find_first_of("(");
  // the length of the name of the param is actually equal to the position
  // of the first paren
  std::string compName=nodalComputation.substr(0,firstParen);
  std::string args=nodalComputation.substr(firstParen+1,nodalComputation.length()-compName.length()-2);

#ifdef Xyce_DEBUG_EXPRESSION
  Xyce::dout() << "Processing nodalComputation : " << nodalComputation
               << Util::push<< std::endl;
  Xyce::dout() << "name of computation: " << compName << std::endl;
  Xyce::dout() << "args: " << args << Util::push << std::endl;
#endif // Xyce_DEBUG_EXPRESSION

  std::size_t firstComma=args.find_first_of(",");
  while (firstComma != std::string::npos)
  {
    std::string arg = args.substr(0,firstComma);
    std::size_t argsLength = args.length();
    args = args.substr(firstComma+1,argsLength-arg.length()-1);
    firstComma = args.find_first_of(",");
    tempParamList.push_back(Param(arg,0.0));
#ifdef Xyce_DEBUG_EXPRESSION
    Xyce::dout() << "arg " << arg << std::endl;
#endif
  }

  tempParamList.push_back(Param(args,0.0));

#ifdef Xyce_DEBUG_EXPRESSION
  Xyce::dout() << "Remaining arg " << args << std::endl;
  Xyce::dout() << "There were " << tempParamList.size() << " args."
               << Util::pop << std::endl;
  Xyce::dout() << Util::pop << std::endl;
#endif

  paramList.push_back(Param(compName,static_cast<double>(tempParamList.size())));
  std::copy (tempParamList.begin(), tempParamList.end(), std::back_inserter(paramList));

} // namespace (unnammed)

}
} // namespace Util
} // namespace Xyce
