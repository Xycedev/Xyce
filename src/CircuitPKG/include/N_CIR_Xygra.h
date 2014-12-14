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
// Filename       : $RCSfile: N_CIR_Xygra.h,v $
//
// Purpose        : Provide a class for Xyce/Alegra coupling
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 8/21/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/04/16 23:47:41 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_XYGRA_H
#define Xyce_N_CIR_XYGRA_H

#include <map>
#include <string>

#include <N_CIR_Xyce.h>

#include <N_DEV_fwd.h>

//-----------------------------------------------------------------------------
// Class         : N_CIR_Xygra
// Purpose       :
// Special Notes :
//
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/21/08
//-----------------------------------------------------------------------------
class N_CIR_Xygra : public Xyce::Circuit::Simulator
{
  public:
    N_CIR_Xygra()
    {}
    
    virtual ~N_CIR_Xygra()
    {}
    
    int xygraGetNumNodes(const std::string & deviceName);
    int xygraGetNumWindings(const std::string & deviceName);
    void xygraGetCoilWindings(const std::string & deviceName, std::vector<int> & cW);
    void xygraGetCoilNames(const std::string & deviceName, std::vector<std::string> & cN);
    bool xygraSetConductances(const std::string & deviceName, const std::vector< std::vector<double> > & cM);
    bool xygraSetK(const std::string & deviceName, const std::vector< std::vector<double> > & kM, const double t=0);
    bool xygraSetSources(const std::string & deviceName, const std::vector< double > & sV, const double t=0);
    bool xygraGetVoltages(const std::string & deviceName, std::vector< double > & vN);

    Xyce::Device::Xygra::Instance *getXygraInstance_(const std::string & deviceName);

  private:
    std::map<std::string, Xyce::Device::Xygra::Instance *> xygraDeviceMap_;
};

#endif
