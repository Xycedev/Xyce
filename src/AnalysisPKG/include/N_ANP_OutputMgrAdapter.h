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
// Filename       : $RCSfile: N_ANP_OutputMgrAdapter.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.52.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:48 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_OutputMgrAdapter_h
#define Xyce_N_ANP_OutputMgrAdapter_h

#include <Teuchos_SerialDenseMatrix.hpp>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : OutputMgrAdapter
// Purpose       : Inteface class for the output manager
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class OutputMgrAdapter : public Util::ListenerAutoSubscribe<StepEvent>
{
public:
  OutputMgrAdapter(Parallel::Machine comm, Util::Notifier<StepEvent> &step_notifier);

  virtual ~OutputMgrAdapter();

  void notify(const StepEvent &event);

  void registerOutputMgr(IO::OutputMgr * outputMgrPtr)
  {
    outputManager_ = outputMgrPtr;
  }

  void addOutputResults(const Util::OptionBlock & option_block);
  void addOutputResponse(const Util::OptionBlock & option_block);
  
  void setStepParamVec(const std::vector<SweepParam> & paramVec);

  void setDCParamVec(const std::vector<SweepParam> & paramVec);

  const std::vector<SweepParam> &getStepParamVec() const
  {
    return *stepParamVector_;
  }

  const std::vector<SweepParam> &getDCParamVec() const
  {
    return dcParamVector_;
  }

  // accessor methods
  int getStepAnalysisStepNumber() const
  {
    return stepAnalysisStepNumber_;
  }

  int getDCAnalysisStepNumber() const
  {
    return dcAnalysisStepNumber_;
  }

  int getDCAnalysisMaxSteps() const
  {
    return dcAnalysisMaxSteps_;
  }

  IO::OutputMgr &getOutputManager()
  {
    return *outputManager_;
  }

  void setDCAnalysisStepNumber( int num )
  {
    dcAnalysisStepNumber_ = num;
  }

  void setDCAnalysisMaxSteps( int num )
  {
    dcAnalysisMaxSteps_ = num;
  }

  void tranOutput(
    double time, N_LAS_Vector & currSolutionPtr, N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr, 
    std::vector<double> & objectiveVec_, 
    std::vector<double> & dOdpVec_, 
    std::vector<double> & dOdpAdjVec_,
    std::vector<double> & scaled_dOdpVec_, 
    std::vector<double> & scaled_dOdpAdjVec_,
    bool skipPrintLineOutput = false);

  void dcOutput( 
    int dcStepNumber, 
    N_LAS_Vector & currSolutionPtr, N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr,
    std::vector<double> & objectiveVec_, 
    std::vector<double> & dOdpVec_, 
    std::vector<double> & dOdpAdjVec_, 
    std::vector<double> & scaled_dOdpVec_, 
    std::vector<double> & scaled_dOdpAdjVec_);

  void outputResult(const N_LAS_Vector &currSolutionPtr, const N_LAS_Vector &currStatePtr, const N_LAS_Vector &currStorePtr );

  void steppingComplete();

  void finishOutput();

  bool setupInitialConditions ( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec);

  void outputDCOP(const N_LAS_Vector &solution);

  void outputMPDE( double time, const std::vector<double> &fast_time_points, const N_LAS_Vector & solution_vector);

  void outputHB(
    const std::vector< double > & timePoints, const std::vector< double > & freqPoints,
    const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
    const N_LAS_BlockVector & freqDomainSolnVecImaginary, const N_LAS_BlockVector & timeDomainStoreVec,
    const N_LAS_BlockVector & freqDomainStoreVecReal, const N_LAS_BlockVector & freqDomainStoreVecImaginary);

  void outputAC(double freq, const N_LAS_Vector & solnVecRealPtr, const N_LAS_Vector & solnVecImaginaryPtr);

  void outputMORTF(bool origSys, double freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H );

  void resetOutputMORTF();

  void outputROM(
    const Teuchos::SerialDenseMatrix<int, double>& Ghat, const Teuchos::SerialDenseMatrix<int, double>& Chat,
    const Teuchos::SerialDenseMatrix<int, double>& Bhat, const Teuchos::SerialDenseMatrix<int, double>& Lhat );

  void outputROM(
    const N_LAS_Matrix& Ghat, const N_LAS_Matrix& Chat,
    const Teuchos::SerialDenseMatrix<int, double>& Bhat,
    const Teuchos::SerialDenseMatrix<int, double>& Lhat );

  bool getOutputIntervals(double & initialInterval, std::vector<std::pair< double, double > > * intervalPairs) const;

  void outputHomotopy( const std::vector<std::string> & paramNames, const std::vector<double> & paramVals, N_LAS_Vector & solnVecPtr );

  const Xyce::NodeNamePairMap & getAllNodes() const;

private:
  Parallel::Machine             comm_;
  IO::OutputMgr *               outputManager_;
  IO::Measure::Manager *        measureManager_;
  IO::FourierMgr *              fourierManager_;
  IO::OutputMOR *               outputMOR_;
  IO::OutputResults *           outputResults_;
  IO::OutputResponse *          outputResponse_;

  std::vector<SweepParam>                     emptyParamVector_;
  const std::vector<SweepParam> *             stepParamVector_;
  std::vector<SweepParam>                     dcParamVector_;

  int stepAnalysisStepNumber_;
  int stepAnalysisMaxSteps_;
  int dcAnalysisStepNumber_;
  int dcAnalysisMaxSteps_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::OutputMgrAdapter N_ANP_OutputMgrAdapter;

#endif // Xyce_N_ANP_OutputMgrAdapter_h
