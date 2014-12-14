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
// Filename       : $RCSfile: N_DEV_Specie.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Lawrence C Musson, SNL
//
// Creation Date  : 04/16/2014
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/05/01 22:27:13 $
//
// Current Owner  : $Author: lcmusso $
//-----------------------------------------------------------------------------

#include "N_DEV_Specie.h"
#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <string>
#include <N_DEV_Const.h>
#include <vector>
#include <iostream>

namespace Xyce {
  namespace Device {


    //-----------------------------------------------------------------------------
    // Function      : Specie::setBCEnhancedDiffusion
    // Purpose       : sets up bourgoin corbett enhancement
    // Special Notes : 
    // Scope         : public
    // Creator       : Lawrence C Musson, SNL
    // Creation Date : 04/16/2014
    //-----------------------------------------------------------------------------

    void Specie::setBCEnhancedDiffusion(int cI, double s, double tV, double hL)
    {
      enhancedDiffusion = true;
      sigma = s;
      thermalVelocity = tV;
      hopLength = hL;
      carrierIndex = cI;
    }

    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------


    //-----------------------------------------------------------------------------
    // Function      : Specie::getDiffusionCoefficient
    // Purpose       : compute the diffusion coefficient included bourgoin corbett enhancement
    // Special Notes : 
    // Scope         : public
    // Creator       : Lawrence C Musson, SNL
    // Creation Date : 04/16/2014
    //-----------------------------------------------------------------------------

    double Specie::getDiffusionCoefficient(double Temperature, 
                                           std::vector<double> &concs,std::vector<double> &constant_vec)
    {

      double DF = DiffusionPrefactor*exp(-ActivationEnergy/(CONSTboltz*Temperature/CONSTQ));


      double preDF = DF;

      if(enhancedDiffusion)
        DF += sigma*thermalVelocity*constant_vec[carrierIndex]*hopLength*hopLength/6.0;

      return DF;

    }

  }; //namespace Device
}; //namespace Xyce
