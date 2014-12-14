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
// Filename       : $RCSfile: N_ANP_fwd.h,v $
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
// Revision Number: $Revision: 1.7 $
// Revision Date  : $Date: 2014/08/05 21:41:35 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_fwd_h
#define Xyce_N_ANP_fwd_h

#include <vector>

namespace Xyce {
namespace Analysis {

enum Analysis_Mode
{
  ANP_MODE_INVALID,
  ANP_MODE_DC_OP,
  ANP_MODE_DC_SWEEP,
  ANP_MODE_TRANSIENT,
  ANP_MODE_MPDE,
  ANP_MODE_HB,
  ANP_MODE_AC,
  ANP_MODE_MOR
};

enum CurrentMode {CURRENT_MODE_TRANOP = 0, CURRENT_MODE_TRANSIENT = 1, CURRENT_MODE_DC_SWEEP = 2};

class ModelEvaluator;
class ModelEvaluator_Stateless;
class OutputMgrAdapter;
class SweepParam;

class AC;
class AnalysisBase;
class AnalysisManager;
class Transient;
class Step;
class DCSweep;
class MPDE;
class HB;
class Dakota;
class MOR;

typedef std::vector<SweepParam> SweepVector;

} // namespace Analysis
} // namespace Xyce


typedef Xyce::Analysis::ModelEvaluator N_ANP_ModelEvaluator;
typedef Xyce::Analysis::ModelEvaluator_Stateless N_ANP_ModelEvaluator_Stateless;
typedef Xyce::Analysis::OutputMgrAdapter N_ANP_OutputMgrAdapter;
typedef Xyce::Analysis::SweepParam N_ANP_SweepParam;


typedef Xyce::Analysis::AC N_ANP_AC;
typedef Xyce::Analysis::AnalysisBase N_ANP_AnalysisBase;
typedef Xyce::Analysis::AnalysisManager N_ANP_AnalysisManager;
typedef Xyce::Analysis::Transient N_ANP_Transient;
typedef Xyce::Analysis::Step N_ANP_Step;
typedef Xyce::Analysis::DCSweep N_ANP_DCSweep;
typedef Xyce::Analysis::MPDE N_ANP_MPDE;
typedef Xyce::Analysis::HB N_ANP_HB;
typedef Xyce::Analysis::Dakota N_ANP_Dakota;
typedef Xyce::Analysis::MOR N_ANP_MOP;

#endif // Xyce_N_ANP_fwd_h
