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
// Filename       : $RCSfile: N_IO_MeasureBase.h,v $
//
// Purpose        : Base class for Measure types
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.32.2.2 $
//
// Revision Date  : $Date: 2014/08/28 22:41:46 $
//
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureBase_h
#define Xyce_N_IO_MeasureBase_h

#include <string>
#include <list>
#include <vector>

#include <N_ANP_SweepParam.h>
#include <N_IO_fwd.h>
#include <N_LAS_Vector.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {
namespace Measure {

typedef std::list<Util::Param> ParameterList;

//-------------------------------------------------------------------------
// Class         : MeasureBase
// Purpose       : Base class for common analysis functions
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Base
{
  public:
    Base( const Util::OptionBlock & measureBlock);

    virtual ~Base();

    virtual void prepareOutputVariables() = 0;
    virtual void reset() {initialized_ = false;}
    virtual void updateTran(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector * stateVec, const N_LAS_Vector *storeVec) {}
    virtual void updateDC(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec ) {}
    virtual void updateAC(Parallel::Machine comm, const double frequency, const N_LAS_Vector *solnVec, const N_LAS_Vector *imaginaryVec) {}

protected:
  // used by individual measure classes to update the output variables 
  // on which they depend 
    void updateOutputVars(Parallel::Machine comm, std::vector<double> & outputVarVec, const double circuitTime,
                          const N_LAS_Vector *solnVec, const N_LAS_Vector * stateVec, const N_LAS_Vector *storeVec, const N_LAS_Vector *imaginaryVec );

public:
  bool finishedCalculation() {
      return calculationDone_;
    }

    void makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);

    // these functions implement measurement window criterial
    // such as TD (delay time) and Rise/Fall/Cross counts
    bool withinTransientWindow( double time );
    bool withinRiseFallCrossWindow( double measureVal, double crossVal  );
    bool withinFromToWindow( double time );
    bool withinMinMaxThreash( double value);

    // used to call the output manager's getPrgetImmutableValue<int>()
    double getOutputValue(Parallel::Machine comm, Util::Op::Operator *op, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec, const N_LAS_Vector *imaginaryVec );

    // used to get the measurement result
    virtual double getMeasureResult() {
      return calculationResult_;
    }

    // used to print the measurement result to an output stream object
    virtual std::ostream& printMeasureResult(std::ostream& os)
    {
      os << name_ << " = " << this->getMeasureResult() << std::endl;
      return os;
    }

    // this is the user defined name for this measurement
    std::string name_;

    // this is the mode under which the measurement is active (DC, AC or TRAN)
    std::string mode_;

    // This is the type of measurement to be done.  A text string is stored here
    // but derrived classes do the work.  Here are the types:
    // TRIG TARG                 -- see MeasureRiseFallDelay
    // AVG, MAX, MIN, PP, RMS    -- see MeasureStatistics
    // FIND, WHEN                -- see MeasureFindWhen
    // DERIVATIVE, DERIV         -- see MeasureDerivativeEvaluation
    // INTEGRAL                  -- see Measure Integral Evaluation
    // ERROR                     -- see MeasureRelativeError
    std::string type_;

    // this bool is set by the constructor of the derived class.  So those that
    // are not supported (i.e. not fully implemented) can warn the user and
    // then the measure manager can based on this flag not add them to the active list.
    bool typeSupported_;
    
    // this is a flag indicating if the measure is set up to start doing it's 
    // calculation.  If it's reset to "false" then the measure should start over 
    bool initialized_;

    // this is the list of output variables needed in the measure.
    // it could be one or more values.  Since each var in the Util::Param list
    // can take up multiple spots on the list, keep a count of the number and
    // a vector of iterators pointing to the start of each var in the list
    int numDepSolVars_;
    Util::ParameterList depSolVarIterVector_;
    Util::Op::OpList outputVars_;
    double outputValueTarget_;
    bool  outputValueTargetGiven_;
    double lastOutputValue_;

    // many controls on how the measure calculation is done are set via keyword=val
    // we'll parse those out and hold them in the base class:
    double td_;
    bool tdGiven_;
    double goal_;
    double weight_;
    double minval_;
    double at_;
    double from_;
    bool fromGiven_;
    double to_;
    bool toGiven_;
    double ymin_;
    double ymax_;
    int rise_;
    bool riseGiven_;
    int fall_;
    bool fallGiven_;
    int cross_;
    bool crossGiven_;
    int actualRise_;
    bool isRising_;
    // max and min threshold are used to set upper and lower bounds on a value
    // prior to computation.
    double maxThresh_;
    bool maxThreshGiven_;
    double minThresh_;
    bool minThreshGiven_;

    int actualFall_;
    bool isFalling_;
    int actualCross_;

    // these are used by the Duty, On and Off measures to give values as to
    // when a variable is "on" or "off"
    double onValue_;
    bool onValueGiven_;
    double offValue_;
    bool offValueGiven_;

    // in the Rise/fall/delay measure the trigger and target can be nodes of the
    // circuit so they are of the type Util::Param
    Util::Param trig_;
    Util::Param targ_;
    double trigOutputValueTarget_;
    bool trigOutputValueTargetGiven_;
    double targOutputValueTarget_;
    bool targOutputValueTargetGiven_;
    double trigFracMax_;  // fraction of the maxima for the trigger value
    bool trigFracMaxGiven_;
    double targFracMax_;  // fraction of the maxima for the target value
    bool targFracMaxGiven_;

    // this is used by max to measure the time when a requested fraction of the
    // maximum is reached.  as in 90% of max value.  Also applies to min as well.
    double fractionToExtrema_;
    bool fractionToExtremaGiven_;

    // this is used by fourier to determine how many harmonics to compute and sample grid size.
    int numFreq_;
    int gridSize_;

    // many measurements will finish before the end of a simulation.  This flag is used
    // to indicate if the measurement is done.
    bool calculationDone_;

    // This is where the results are stored.
    double calculationResult_;

    //  used in the Error measure where least squared error fit to data is measured.
    std::string dataFileName_;           // name of file that may hold data
    std::string comparisonFunctionName_; // what type of comparison will be done as in L1NORM, L2NORM or INFNORM 
    int independentVarColumn_;            // if the data file has more then two columns, this is the indpendent var column (defaults to 0)
    int dependentVarColumn_;              // if there are more then two columns in the data file, this is the dependent var (defaults to 1);
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::Base N_IO_MeasureBase;

#endif // Xyce_N_IO_MeasureBase_h

