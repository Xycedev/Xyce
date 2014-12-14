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
// Filename       : $RCSfile: N_IO_Op.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
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
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/08/11 18:03:53 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Op_h
#define Xyce_N_IO_Op_h

#include <iterator>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_Op.h>
#include <N_IO_Measure_fwd.h>
#include <N_UTL_fwd.h>

#include <N_IO_OutputMgr.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_UTL_Op.h>
#include <N_UTL_ExpressionData.h>

namespace Xyce {
namespace IO {

class CurrentIndexOp : public Util::Op::Op<CurrentIndexOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  CurrentIndexOp(const std::string &name)
    : Base(name)
  {}

  virtual ~CurrentIndexOp()
  {}

  static complex get(const CurrentIndexOp &op, const Util::Op::OpData &op_data);
};

class OutputMgrTimeOp : public Util::Op::Op<OutputMgrTimeOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrTimeOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrTimeOp()
  {}

  static complex get(const OutputMgrTimeOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

class OutputMgrFrequencyOp : public Util::Op::Op<OutputMgrFrequencyOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrFrequencyOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrFrequencyOp()
  {}

  static complex get(const OutputMgrFrequencyOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};


class OutputMgrTemperatureOp : public Util::Op::Op<OutputMgrTemperatureOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrTemperatureOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrTemperatureOp()
  {}

  static complex get(const OutputMgrTemperatureOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

class OutputMgrStepSweepOp : public Util::Op::Op<OutputMgrStepSweepOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrStepSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    : Base(name),
      outputMgr_(output_manager),
      index_(index)
  {}

  virtual ~OutputMgrStepSweepOp()
  {}

  static complex get(const OutputMgrStepSweepOp &op, const Util::Op::OpData &op_data);

  const int           index_;
  const OutputMgr &   outputMgr_;
};

class OutputMgrDCSweepOp : public Util::Op::Op<OutputMgrDCSweepOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrDCSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    : Base(name),
      outputMgr_(output_manager),
      index_(index)
  {}

  virtual ~OutputMgrDCSweepOp()
  {}

  static complex get(const OutputMgrDCSweepOp &op, const Util::Op::OpData &op_data);

  const int           index_;
  const OutputMgr &   outputMgr_;
};

class OutputMgrDCSweepCurrentValueOp : public Util::Op::Op<OutputMgrDCSweepCurrentValueOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrDCSweepCurrentValueOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrDCSweepCurrentValueOp()
  {}

  static complex get(const OutputMgrDCSweepCurrentValueOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

class ObjectiveOp : public Util::Op::Op<ObjectiveOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  ObjectiveOp(const std::string &name, Parallel::Machine comm, const OutputMgr &output_manager, Objective &objective)
    : Base(name),
      comm_(comm),
      outputMgr_(output_manager),
      objective_(objective)
  {}

  virtual ~ObjectiveOp()
  {}

  static complex get(const ObjectiveOp &op, const Util::Op::OpData &op_data);

  Objective &           objective_;
  Parallel::Machine     comm_;
  const OutputMgr &     outputMgr_;
};

class SolutionOp : public Util::Op::Op<SolutionOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  SolutionOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~SolutionOp()
  {}

  static complex get(const SolutionOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

class SolutionRealOp : public Util::Op::Op<SolutionRealOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  SolutionRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~SolutionRealOp()
  {}

  static complex get(const SolutionRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class SolutionImaginaryOp : public Util::Op::Op<SolutionImaginaryOp, Util::Op::ReduceSum>
{
public:
  SolutionImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionImaginaryOp()
  {}

  static complex get(const SolutionImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class SolutionMagnitudeOp : public Util::Op::Op<SolutionMagnitudeOp, Util::Op::ReduceSum>
{
public:
  SolutionMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionMagnitudeOp()
  {}

  static complex get(const SolutionMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class SolutionPhaseOp : public Util::Op::Op<SolutionPhaseOp, Util::Op::ReduceSum>
{
public:
  SolutionPhaseOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionPhaseOp()
  {}

  static complex get(const SolutionPhaseOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class SolutionDecibelsOp : public Util::Op::Op<SolutionDecibelsOp, Util::Op::ReduceSum>
{
public:
  SolutionDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionDecibelsOp()
  {}

  static complex get(const SolutionDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class VoltageDifferenceOp : public Util::Op::Op<VoltageDifferenceOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  VoltageDifferenceOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceOp()
  {}

  static complex get(const VoltageDifferenceOp &op, const Util::Op::OpData &op_data);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceRealOp : public Util::Op::Op<VoltageDifferenceRealOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  VoltageDifferenceRealOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceRealOp()
  {}

  static complex get(const VoltageDifferenceRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceImaginaryOp : public Util::Op::Op<VoltageDifferenceImaginaryOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceImaginaryOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceImaginaryOp()
  {}

  static complex get(const VoltageDifferenceImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceMagnitudeOp : public Util::Op::Op<VoltageDifferenceMagnitudeOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceMagnitudeOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceMagnitudeOp()
  {}

  static complex get(const VoltageDifferenceMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferencePhaseOp : public Util::Op::Op<VoltageDifferencePhaseOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferencePhaseOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferencePhaseOp()
  {}

  static complex get(const VoltageDifferencePhaseOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceDecibelsOp : public Util::Op::Op<VoltageDifferenceDecibelsOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceDecibelsOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceDecibelsOp()
  {}

  static complex get(const VoltageDifferenceDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class StateOp : public Util::Op::Op<StateOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StateOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StateOp()
  {}

  static complex get(const StateOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

class StoreOp : public Util::Op::Op<StoreOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StoreOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreOp()
  {}

  static complex get(const StoreOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class StoreRealOp : public Util::Op::Op<StoreRealOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StoreRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~StoreRealOp()
  {}

  static complex get(const StoreRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class StoreImaginaryOp : public Util::Op::Op<StoreImaginaryOp, Util::Op::ReduceSum>
{
public:
  StoreImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreImaginaryOp()
  {}

  static complex get(const StoreImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class StoreMagnitudeOp : public Util::Op::Op<StoreMagnitudeOp, Util::Op::ReduceSum>
{
public:
  StoreMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreMagnitudeOp()
  {}

  static complex get(const StoreMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class StorePhaseOp : public Util::Op::Op<StorePhaseOp, Util::Op::ReduceSum>
{
public:
  StorePhaseOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StorePhaseOp()
  {}

  static complex get(const StorePhaseOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class StoreDecibelsOp : public Util::Op::Op<StoreDecibelsOp, Util::Op::ReduceSum>
{
public:
  StoreDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreDecibelsOp()
  {}

  static complex get(const StoreDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class SensitivityObjFunctionOp : public Util::Op::Op<SensitivityObjFunctionOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivityObjFunctionOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivityObjFunctionOp()
  {}

  static complex get(const SensitivityObjFunctionOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpDirectOp : public Util::Op::Op<SensitivitydOdpDirectOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpDirectOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpDirectOp()
  {}

  static complex get(const SensitivitydOdpDirectOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpDirectScaledOp : public Util::Op::Op<SensitivitydOdpDirectScaledOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpDirectScaledOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpDirectScaledOp()
  {}

  static complex get(const SensitivitydOdpDirectScaledOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpAdjointOp : public Util::Op::Op<SensitivitydOdpAdjointOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpAdjointOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpAdjointOp()
  {}

  static complex get(const SensitivitydOdpAdjointOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpAdjointScaledOp : public Util::Op::Op<SensitivitydOdpAdjointScaledOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpAdjointScaledOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpAdjointScaledOp()
  {}

  static complex get(const SensitivitydOdpAdjointScaledOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class MeasureOp : public Util::Op::Op<MeasureOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  MeasureOp(const std::string &name, const Measure::Base &measure)
    : Base(name),
      measure_(measure)
  {}

  virtual ~MeasureOp()
  {}

  static complex get(const MeasureOp &op, const Util::Op::OpData &op_data);

  const Measure::Base &       measure_;
};

class ExpressionOp : public Util::Op::Op<ExpressionOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  ExpressionOp(const std::string &name, Util::Expression &expression, Parallel::Machine comm, const OutputMgr &output_manager);
  ExpressionOp(const std::string &name, const std::string &expression, Parallel::Machine comm, const OutputMgr &output_manager);

  virtual ~ExpressionOp()
  {}

  void init(Parallel::Machine comm, const OutputMgr &output_manager);

  static complex get(const ExpressionOp &op, const Util::Op::OpData &op_data);

  mutable Util::ExpressionData          expressionData_;
  Parallel::Machine                     comm_;
  const OutputMgr &                     outputMgr_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Op_h
