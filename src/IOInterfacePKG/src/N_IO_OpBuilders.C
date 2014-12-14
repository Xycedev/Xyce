//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_IO_OpBuilders.C,v $
//
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2014/08/28 22:37:43 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

#include <N_IO_Op.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

typedef std::list<Util::Param> ParameterList;

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

} // namespace <unnamed>

struct CircuitTemperatureOpBuilder : public Util::Op::Builder
{
  CircuitTemperatureOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitTemperatureOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrTemperatureOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "TEMP") {
      new_op  = new OutputMgrTemperatureOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct CircuitTimeOpBuilder : public Util::Op::Builder
{
  CircuitTimeOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitTimeOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrTimeOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "TIME") {
      new_op  = new OutputMgrTimeOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct CircuitFrequencyOpBuilder : public Util::Op::Builder
{
  CircuitFrequencyOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitFrequencyOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrFrequencyOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "FREQUENCY") {
      new_op  = new OutputMgrFrequencyOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct StepSweepOpBuilder : public Util::Op::Builder
{
  StepSweepOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~StepSweepOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrStepSweepOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    for (size_t i = 0; i < outputManager_.getStepParamVec().size(); ++i)
    {
      if (param_tag == outputManager_.getStepParamVec()[i].name)
      {
        new_op = new OutputMgrStepSweepOp(param_tag, outputManager_, i);
        break;
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct DCSweepOpBuilder : public Util::Op::Builder
{
  DCSweepOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~DCSweepOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrDCSweepOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    for (size_t i = 0; i < outputManager_.getDCParamVec().size(); ++i)
    {
      if (param_tag == outputManager_.getDCParamVec()[i].name)
      {
        new_op = new OutputMgrDCSweepOp(param_tag, outputManager_, i);
        break;
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct DCSweepCurrentValueOpBuilder : public Util::Op::Builder
{
  DCSweepCurrentValueOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~DCSweepCurrentValueOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrDCSweepCurrentValueOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "sweep") {
      new_op  = new OutputMgrDCSweepCurrentValueOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct InternalVariableOpBuilder : public Util::Op::Builder
{
  InternalVariableOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~InternalVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<StateOp>();
    builder_manager.addCreateFunction<StoreOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'N' && args.size() == 1)
    {
      NodeNamePairMap::const_iterator it = findNode(args[0], outputManager_.getAllNodes(), outputManager_.getAliasNodeMap());
      if (it != outputManager_.getAllNodes().end())
      {
        int index = (*it).second.first;
        if (param_tag == "N" )
        {
          new_op = new SolutionOp(name, index);
        }
        else if (param_tag == "NR" )
        {
          new_op = new SolutionOp(name, index);
        }
        else if (param_tag == "NI" )
        {
          new_op = new SolutionImaginaryOp(name, index);
        }
        else if (param_tag == "NM" )
        {
          new_op = new SolutionMagnitudeOp(name, index);
        }
        else if (param_tag == "NP" )
        {
          new_op = new SolutionPhaseOp(name, index);
        }
        else if (param_tag == "NDB" )
        {
          new_op = new SolutionDecibelsOp(name, index);
        }
      }
      else
      {
        it = findNode(args[0], outputManager_.getStateNodes(), outputManager_.getAliasNodeMap());
        if (it != outputManager_.getStateNodes().end())
        {
          new_op = new StateOp(name, (*it).second.first);
        }
        else
        {
          it = outputManager_.getStoreNodes().find(args[0]);
          if (it != outputManager_.getStoreNodes().end())
          {
            new_op = new StoreOp(name, (*it).second.first);
          }
          else
          {
            new_op = new Util::Op::UndefinedOp(param_tag);
          }
        }
      }

      if (new_op)
        new_op->addArg(args[0]);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct VoltageVariableOpBuilder : public Util::Op::Builder
{
  VoltageVariableOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~VoltageVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionRealOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<VoltageDifferenceOp>();
    builder_manager.addCreateFunction<VoltageDifferenceRealOp>();
    builder_manager.addCreateFunction<VoltageDifferenceImaginaryOp>();
    builder_manager.addCreateFunction<VoltageDifferenceMagnitudeOp>();
    builder_manager.addCreateFunction<VoltageDifferencePhaseOp>();
    builder_manager.addCreateFunction<VoltageDifferenceDecibelsOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'V' && args.size() > 0)
    {

      // Solution variable
      if (args.size() == 1)
      {
        int index = -1;
        NodeNamePairMap::const_iterator it = findNode(args[0], outputManager_.getAllNodes(), outputManager_.getAliasNodeMap());
        if (it != outputManager_.getAllNodes().end())
          index = (*it).second.first;

        if (param_tag == "V" )
        {
          new_op = new SolutionOp(name, index);
        }
        else if (param_tag == "VR" )
        {
          new_op = new SolutionRealOp(name, index);
        }
        else if (param_tag == "VI" )
        {
          new_op = new SolutionImaginaryOp(name, index);
        }
        else if (param_tag == "VM" )
        {
          new_op = new SolutionMagnitudeOp(name, index);
        }
        else if (param_tag == "VP" )
        {
          new_op = new SolutionPhaseOp(name, index);
        }
        else if (param_tag == "VDB" )
        {
          new_op = new SolutionDecibelsOp(name, index);
        }

        if (new_op)
          new_op->addArg(args[0]);
      }

      // Volatage Difference
      if (args.size() == 2)
      {
        int index1 = -1;
        int index2 = -1;

        NodeNamePairMap::const_iterator it = findNode(args[0], outputManager_.getAllNodes(), outputManager_.getAliasNodeMap());
        if (it != outputManager_.getAllNodes().end())
          index1 = (*it).second.first;

        it = findNode(args[1], outputManager_.getAllNodes(), outputManager_.getAliasNodeMap());
        if (it != outputManager_.getAllNodes().end())
          index2 = (*it).second.first;

        if (param_tag == "V" )
        {
          new_op = new VoltageDifferenceOp(name, index1, index2);
        }
        else if (param_tag == "VR" )
        {
          new_op = new VoltageDifferenceRealOp(name, index1, index2);
        }
        else if (param_tag == "VI" )
        {
          new_op = new VoltageDifferenceImaginaryOp(name, index1, index2);
        }
        else if (param_tag == "VM" )
        {
          new_op = new VoltageDifferenceMagnitudeOp(name, index1, index2);
        }
        else if (param_tag == "VP" )
        {
          new_op = new VoltageDifferencePhaseOp(name, index1, index2);
        }
        else if (param_tag == "VDB" )
        {
          new_op = new VoltageDifferenceDecibelsOp(name, index1, index2);
        }

        if (new_op)
          new_op->addArgs(args.begin(), args.end());
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

namespace {
static const char * const func_names[] = {"II", "IR", "IP", "IM", "IDB"};
}

struct CurrentVariableOpBuilder : public Util::Op::Builder
{

  CurrentVariableOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CurrentVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionRealOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<StoreOp>();
    builder_manager.addCreateFunction<StoreRealOp>();
    builder_manager.addCreateFunction<StoreImaginaryOp>();
    builder_manager.addCreateFunction<StoreMagnitudeOp>();
    builder_manager.addCreateFunction<StorePhaseOp>();
    builder_manager.addCreateFunction<StoreDecibelsOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'I')
    {
      // Node name could be circuit_context:DeviceTypeDeviceName while internally it should be DeviceType:circuit_context:DeviceName.
      std::string modifiedName;
      std::string::size_type lastColonInName = args[0].find_last_of(":");
      if ((lastColonInName != std::string::npos) && (lastColonInName + 1 < args[0].length()))
      {
        std::string::iterator deviceName = args[0].begin() + lastColonInName+1;
        std::string::iterator namePrefixEnd = args[0].begin() + lastColonInName;
        modifiedName.append(deviceName, deviceName + 1);
        modifiedName.append(":");
        modifiedName.append(args[0].begin(), namePrefixEnd + 1);
        modifiedName.append(deviceName + 1, args[0].end());
      }
      else
      {
        modifiedName = args[0];
      }

      // could be a device lead current "DEV_I" or a branch current.
      // so we don't have to duplicate solution vars(branch currents) in the
      // store vector, look for each type.
      bool param_func = std::find(func_names, func_names + sizeof(func_names)/sizeof(func_names[0]), param_tag) != func_names + sizeof(func_names)/sizeof(func_names[0]);
      std::string store_name = modifiedName + ":DEV_" + (param_func ? "I" : param_tag);  // if it is in the state/store vec.

      // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
      // we try to find devices based on store_name in the following blocks of code,
      // so do this modification now.
      std::string::size_type space = store_name.find_first_of(" ");
      if (space != std::string::npos)
      {
        if (space == 4 && store_name.substr(0, 4) == "YPDE")
        {
          store_name.replace(4, 1, ":");
        }
      }

      // Search store
      NodeNamePairMap::const_iterator it = outputManager_.getStoreNodes().find(store_name);
      if (it != outputManager_.getStoreNodes().end())
      {
        int index = (*it).second.first;
        if (param_tag == "IR" )
        {
          new_op = new StoreRealOp(name, index);
        }
        else if (param_tag == "II" )
        {
          new_op = new StoreImaginaryOp(name, index);
        }
        else if (param_tag == "IM" )
        {
          new_op = new StoreMagnitudeOp(name, index);
        }
        else if (param_tag == "IP" )
        {
          new_op = new StorePhaseOp(name, index);
        }
        else if (param_tag == "IDB" )
        {
          new_op = new StoreDecibelsOp(name, index);
        }
        else // IC, IE, IB
        {
          new_op = new StoreOp(name, index);
        }
      }

      // Search solution
      if (!new_op)
      {
        std::string solution_name = modifiedName + "_BRANCH";         // if it is in the solution vec.
        it = outputManager_.getAllNodes().find(solution_name);
        if (it != outputManager_.getAllNodes().end())
        {
          int index = (*it).second.first;
          if (param_tag == "I" )
          {
            new_op = new SolutionOp(name, index);
          }
          else if (param_tag == "IR" )
          {
            new_op = new SolutionRealOp(name, index);
          }
          else if (param_tag == "II" )
          {
            new_op = new SolutionImaginaryOp(name, index);
          }
          else if (param_tag == "IM" )
          {
            new_op = new SolutionMagnitudeOp(name, index);
          }
          else if (param_tag == "IP" )
          {
            new_op = new SolutionPhaseOp(name, index);
          }
          else if (param_tag == "IDB" )
          {
            new_op = new SolutionDecibelsOp(name, index);
          }
        }
      }

      if (new_op)
        new_op->addArg(args[0]);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

struct ExpressionOpBuilder : public Util::Op::Builder
{
  ExpressionOpBuilder(Parallel::Machine comm, const OutputMgr & output_manager)
    : comm_(comm),
      outputManager_(output_manager)
  {}

  virtual ~ExpressionOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<ExpressionOp>();
    builder_manager.addCreateFunction<Util::Op::ConstantOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    int param_type = (*it).getType();

    if (Util::hasExpressionTag(*it))
    {
      if (param_type == Util::EXPR)
      {
        Util::Param &param = const_cast<Util::Param &>(*it);

        new_op = new ExpressionOp(param_tag, param.getValue<Util::Expression>(), comm_, outputManager_);
      }
      else if( (param_type == Util::DBLE) || (param_type == Util::INT) )
      {
        new_op = new Util::Op::ConstantOp(param_tag, (*it).getImmutableValue<double>());
      }
      else
      {
        new_op = new ExpressionOp(param_tag, param_tag, comm_, outputManager_);
      }
    }
    
    return new_op;
  }

private:
  const Parallel::Machine       comm_;
  const OutputMgr &             outputManager_;
};

struct ObjectiveOpBuilder : public Util::Op::Builder
{
  ObjectiveOpBuilder(Parallel::Machine comm, const OutputMgr & output_manager)
    : comm_(comm),
      outputManager_(output_manager)
  {}

  virtual ~ObjectiveOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<ObjectiveOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    int param_type = (*it).getType();

    if (outputManager_.getObjectiveMap().find(param_tag) != outputManager_.getObjectiveMap().end())
    {
      Objective &objective = outputManager_.getObjective(param_tag);
      new_op = new ObjectiveOp(param_tag, comm_, outputManager_, objective);
    }

    return new_op;
  }

private:
  const Parallel::Machine       comm_;
  const OutputMgr &             outputManager_;
};


struct MeasurementOpBuilder : public Util::Op::Builder
{
  MeasurementOpBuilder(const Measure::Manager & measure_manager)
    : measureManager_(measure_manager)
  {}

  virtual ~MeasurementOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<MeasureOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    const Measure::Base *measure = measureManager_.find(param_tag);
    if (measure)
    {
      new_op = new MeasureOp(param_tag, *measure);
    }

    return new_op;
  }

private:
  const Measure::Manager &             measureManager_;
};

struct CircuitIndexOpBuilder : public Util::Op::Builder
{
  CircuitIndexOpBuilder()
  {}

  virtual ~CircuitIndexOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<CurrentIndexOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "INDEX") {
      new_op  = new CurrentIndexOp(param_tag);
    }

    return new_op;
  }
};

struct SensitivityOpBuilder : public Util::Op::Builder
{
  SensitivityOpBuilder()
  {}

  virtual ~SensitivityOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SensitivityObjFunctionOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectScaledOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointScaledOp>();
  }

  virtual Util::Op::Operator *makeOp(ParameterList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "SENS")
    {
      std::string function;
      std::string parameter;
      ptrdiff_t type;
      int index;

      Util::Marshal min(param_string);

      min >> function >> parameter >> type >> index;

      std::string name = "d" + function + "/d(" + parameter + ")";
      if (type == Util::Op::identifier<SensitivityObjFunctionOp>())
        new_op = new SensitivityObjFunctionOp(function, 0);
      else if (type == Util::Op::identifier<SensitivitydOdpDirectOp>())
        new_op = new SensitivitydOdpDirectOp(name + "_Dir", index);
      else if (type == Util::Op::identifier<SensitivitydOdpDirectScaledOp>())
        new_op = new SensitivitydOdpDirectScaledOp(name + "_Dir_scaled", index);
      else if (type == Util::Op::identifier<SensitivitydOdpAdjointOp>())
        new_op = new SensitivitydOdpAdjointOp(name + "_Adj", index);
      if (type == Util::Op::identifier<SensitivitydOdpAdjointScaledOp>())
        new_op = new SensitivitydOdpAdjointScaledOp(name + "_Adj_scaled", index);
    }

    return new_op;
  }
};

void registerOpBuilders(Util::Op::BuilderManager &op_builder_manager, Parallel::Machine comm, OutputMgr &output_manager)
{
  op_builder_manager.addBuilder(new CircuitTemperatureOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitTimeOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitFrequencyOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitIndexOpBuilder());
  op_builder_manager.addBuilder(new SensitivityOpBuilder());
  op_builder_manager.addBuilder(new ExpressionOpBuilder(comm, output_manager));
  op_builder_manager.addBuilder(new StepSweepOpBuilder(output_manager));
  op_builder_manager.addBuilder(new DCSweepOpBuilder(output_manager));
  op_builder_manager.addBuilder(new DCSweepCurrentValueOpBuilder(output_manager));
  op_builder_manager.addBuilder(new ObjectiveOpBuilder(comm, output_manager));
  op_builder_manager.addBuilder(new InternalVariableOpBuilder(output_manager));
  op_builder_manager.addBuilder(new VoltageVariableOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CurrentVariableOpBuilder(output_manager));
}

void registerOpBuilders(Util::Op::BuilderManager &op_builder_manager, Parallel::Machine comm, Measure::Manager &measure_manager)
{
  op_builder_manager.addBuilder(new MeasurementOpBuilder(measure_manager));
}

} // namespace IO
} // namespace Xyce
