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
// Filename       : $RCSfile: N_LAS_HBPrecondFactory.h,v $
//
// Purpose        : Preconditioner Factory for HB
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 10/01/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/07/15 19:03:35 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_HBPrecondFactory_h
#define Xyce_N_LAS_HBPrecondFactory_h

// ---------- Standard Includes ----------

#include <string>

// ----------   Xyce Includes   ----------

#include <N_LAS_PrecondFactory.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_OptionBlock.h>

// ----------  Fwd Declares  -------------

class N_LAS_Preconditioner;
class N_LAS_Problem;
class N_LAS_HBBuilder;
class N_LOA_HBLoader;
class N_LAS_Builder;
class N_LAS_System;

//-----------------------------------------------------------------------------
// Class         : N_LAS_HBPrecondFactory
// Purpose       : 
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
class N_LAS_HBPrecondFactory : public N_LAS_PrecondFactory
{
public:

  // Default Constructor, sets null preconditioner.
  N_LAS_HBPrecondFactory();

  // Basic Constructor, sets preconditioner factory options.
  N_LAS_HBPrecondFactory( const N_UTL_OptionBlock & OB );

  // Destructor
  virtual ~N_LAS_HBPrecondFactory() {} 

  // Creates a new preconditioner (matrix based).
  // NOTE:  This type of creation is not supported by this preconditioner factory.
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_Problem> & problem )
  {
    Xyce::Report::DevelFatal0().in("N_LAS_HBPrecondFactory::create()") << " using N_LAS_Problem is not supported!";
    return Teuchos::null;
  }

  // Creates a new preconditioner (matrix free).
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_System> & lasSystem );

  // Set the time step(s) being used in the HB analysis.
  // NOTE:  This is only useful for FD preconditioning techniques.
  void setTimeSteps( const std::vector<double> & timeSteps )
    { timeSteps_ = timeSteps; }

  // Set the fast times being used in the HB analysis.
  void setFastTimes( const std::vector<double> & times )
    { times_ = times; }

  // Register the application system builder
  void registerAppBuilder( N_LAS_Builder *appBuilderPtr ) 
    { appBuilderPtr_ = appBuilderPtr; }

  // Register the application system loader
  void registerHBLoader( const Teuchos::RCP<N_LOA_HBLoader>& hbLoaderPtr ) 
    { hbLoaderPtr_ = hbLoaderPtr; }

  // Register the HB builder 
  void registerHBBuilder( const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder ) 
    { hbBuilderPtr_ = hbBuilder; }

private:

  std::vector<double> times_, timeSteps_;
  std::string precType_;
  N_LAS_Builder *               appBuilderPtr_;
  Teuchos::RCP<N_LAS_System> lasSysPtr_;
  Teuchos::RCP<N_LOA_HBLoader> hbLoaderPtr_;
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr_;
  Teuchos::RCP<const N_UTL_OptionBlock> OB_;

  // Copy constructor.
  N_LAS_HBPrecondFactory( const N_LAS_HBPrecondFactory& pf );
};


#endif

