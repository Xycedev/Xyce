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
// Filename       : $RCSfile: N_UTL_Op.C,v $
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
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/08/08 20:02:39 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>
#include <N_ERH_Message.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

namespace Xyce {
namespace Util {
namespace Op {

//-----------------------------------------------------------------------------
// Function      : identifier_base
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 17 17:02:25 2014
//-----------------------------------------------------------------------------
Identifier identifier_base()
{
  static Identifier id = reinterpret_cast<Identifier>(&identifier_base);
  return id;
}

//-----------------------------------------------------------------------------
// Function      : ReduceSum::reduce
// Purpose       : perform an MPI sum reduction operation across all processors.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ReduceSum::reduce(Parallel::Machine comm, complex result)
{
  Parallel::AllReduce(comm, MPI_SUM, &result, 1);

  return result;
}

namespace {

void parameterNameAndArgs(std::string &name, std::vector<std::string> &args, ParameterList::const_iterator &it)
{
  const std::string &param_tag = (*it).tag();

  if (param_tag[0] == 'V' || param_tag[0] == 'I' || param_tag[0] == 'N')
  {
    std::ostringstream oss;
    oss << param_tag << "(";
    int arg_count = (*it).getImmutableValue<int>();
    for (int i = 0; i < arg_count; ++i)
    {
      ++it;
      if (i != 0)
        oss << ",";
      oss << (*it).tag();
      args.push_back((*it).tag());
    }
    oss << ")";
    name = oss.str();
  }
}

}

//-----------------------------------------------------------------------------
// Function      : makeOp
// Purpose       : given a parameter list iterator, construct an Op for
//                 the item so defined
// Special Notes :
// Scope         : global (Xyce::IO namespace, but not in any header)
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
Operator *
makeOp(
  const BuilderManager &      op_builder_manager,
  ParameterList::const_iterator & param_it)
{
  Operator *new_op = op_builder_manager.createOp(param_it);

  if (!new_op)
  {
    std::string param_tag = (*param_it).tag();
    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, param_it);

    new_op = new UndefinedOp(param_tag);
    new_op->addArgs(args.begin(), args.end());
  }

  return new_op;
}

//-----------------------------------------------------------------------------
// Function      : makeOps
// Purpose       : given an output manager, a parameter list
//                (defined by begin and end iterators), and a back_insert
//                iterator for an OpList, create all the the Ops with makeOp,
//                synchronize across processors, validate them,  and if
//                all is well, put them onto the end of the OpList using the
//                back_inserter.
//
// Special Notes :
// Scope         : global, Xyce::IO namespace
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void makeOps(Parallel::Machine comm, const BuilderManager &op_builder_manager, const NetlistLocation &netlist_location, ParameterList::const_iterator begin, ParameterList::const_iterator end, std::back_insert_iterator<OpList> inserter)
{
  if (identifier<ConstantOp>() == identifier<UndefinedOp>())
    Report::DevelFatal0() << "Unique Op identifier function fails, see N_UTL_Op.h, Xyce::identifier<T>";

  if (identifier<UndefinedOp>() != identifier<UndefinedOp>())
    Report::DevelFatal0() << "Consistency of Op identifier function fails, see N_UTL_Op.h, Xyce::identifier<T>";

  OpList ops;

  for (ParameterList::const_iterator it = begin; it != end; ++it)
    ops.push_back(makeOp(op_builder_manager, it));

  // Resolve parallel here
  std::vector<Identifier> op_identifier;
  for (OpList::const_iterator it = ops.begin(); it != ops.end(); ++it)
    op_identifier.push_back((*it)->id());

  Parallel::AllReduce(comm, Parallel::op_identifier_compare_op(), op_identifier);

  // Validate ops and report errors
  std::vector<Identifier>::const_iterator op_identifier_it = op_identifier.begin();
  for (std::vector<Operator *>::iterator it = ops.begin(); it != ops.end(); ++it, ++op_identifier_it)
  {
    if ((*op_identifier_it) == 0) // Defined inconsistently
    {
      Report::UserError().at(netlist_location) << "Function or variable is defined differently on different processors";
    }
    else if ((*it)->id() == identifier<UndefinedOp>()) // not defined on this processor
    {
      std::string name = (*it)->getName();
      const std::vector<std::string> &arg_list = (*it)->getArgs();
      if (!arg_list.empty())
      {
        name += "(";
        for (std::vector<std::string>::const_iterator it = arg_list.begin(); it != arg_list.end(); ++it)
        {
          if (it != arg_list.begin())
            name += ",";
          name += (*it);
        }
        name += ")";
      }
      if ((*op_identifier_it) == identifier<UndefinedOp>()) // not defined anywhere
        Report::UserError().at(netlist_location) << "Function or variable " << name << " is not defined";
      else
      {
        CreateFunction f = op_builder_manager.findCreateFunction(*op_identifier_it);
        delete *it;
        (*it) = f(name);
      }
    }
    else if ((*it)->id() != (*op_identifier_it))
      Report::UserError().at(netlist_location) << "Differing types for " << (*it)->getName() << " discovered across processors";
  }

  std::copy(ops.begin(), ops.end(), inserter);
}

//-----------------------------------------------------------------------------
// Function      : getValue
// Purpose       : Given a communicator, an op, and solution/state/store
//                 vectors, evaluate the op and return its value
// Special Notes : This effectively replaces the old "getPrintValue" method,
//                 and is used throughout the Outputter class.
// Scope         : global, Xyce::IO namespace
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
getValue(
  Parallel::Machine             comm,
  const Util::Op::Operator &    op,
  const OpData &                op_data)
{
  return op(comm, op_data);
}

} // namespace Op
} // namespace Util
} // namespace Xyce
