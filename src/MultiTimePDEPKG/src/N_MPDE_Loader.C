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
// Filename      : $RCSfile: N_MPDE_Loader.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.87 $
// Revision Date  : $Date: 2014/05/22 17:40:32 $
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_MPDE_Loader.h>
#include <N_MPDE_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_Manager.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <Epetra_MultiVector.h>

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::~N_MPDE_Loader
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_MPDE_Loader::~N_MPDE_Loader()
{
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::setFastTimes
// Purpose       : Assign times for fast time scale
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
void N_MPDE_Loader::setFastTimes( const std::vector<double> & times )
{ 
  times_ = times;
  constructPeriodicTimes();
}


void N_MPDE_Loader::setPeriodFlags( const std::vector<bool> & periodicFlags )
{ 
  nonPeriodic_ = periodicFlags;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::constructPeriodicTimes
// Purpose       : Make a copy of times_ vector with padding
//                 at the beginning and end to make the calculation
//                 of derivatives easier around the periodicy condition
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
void N_MPDE_Loader::constructPeriodicTimes()
{
  // we will pad our array of times with 2 times the width, one
  // at the top and one at the bottom of the array
  periodicTimesOffset_ = fastTimeDiscRCPtr_->Width();
  int timesSize = times_.size();
  period_ = times_[timesSize - 1];
  periodicTimes_.resize( timesSize + 2*periodicTimesOffset_ );
  for( int i=0; i< periodicTimesOffset_; i++ )
  {
    periodicTimes_[i] = times_[i + timesSize - periodicTimesOffset_ - 1] - period_;
  }
  for( int i=periodicTimesOffset_; i< (timesSize + periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = times_[i - periodicTimesOffset_];
  }
  for( int i=(timesSize+periodicTimesOffset_); i< (timesSize + 2*periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = period_ - times_[i - (timesSize+periodicTimesOffset_) + 1];
  }
  
#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "periodicTimes_ array is" << std::endl;
  for( int i=0; i<(int)periodicTimes_.size(); i++)
  {
    Xyce::dout() << "  periodicTimes_[ " << i << " ] = " << periodicTimes_[i];
    if( (i >= periodicTimesOffset_) && (i < (timesSize + periodicTimesOffset_) ) )
    {
      Xyce::dout() << "\ttimes_[ " << (i-periodicTimesOffset_) << " ] = " << times_[i-periodicTimesOffset_];
    } 
    Xyce::dout() << std::endl;
  }
#endif
  
}


/*
//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::findNonPeriodicSignals_
// Purpose       : Searches solution vector across periodic boundary
//                 and attempts to identify non-periodic signals
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling 
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::findNonPeriodicSignals_(const N_LAS_BlockVector & solutionBlockVector)
{
  int n2 = times_.size();     // number of points on periodic domain
  double maxAbsDiff = 1.0;
  int numSolutionVars = solutionBlockVector.blockSize();
  
  // now compare soluiton vars at block(0) to block(n2-1).  If
  // these differ by abstol? then mark this as non-periodic
  if( nonPeriodic_.size() != numSolutionVars )
  {
    nonPeriodic_.resize( numSolutionVars );
  }
  
  for (int i=0 ; i<numSolutionVars ; ++i)
  {
    double signalDiff = fabs( solutionBlockVector.block(n2-1)[i] - solutionBlockVector.block(0)[i] );
    if( signalDiff > maxAbsDiff )
    {
      nonPeriodic_[i] = true;
    }
    else
    {
      nonPeriodic_[i] = false;
    }
    
//#ifdef Xyce_DEBUG_MPDE
    Xyce::dout() << "bool N_MPDE_Loader::findNonPeriodicSignals_() For solution variable " << i
      << " The difference over the periodic boundary is " << signalDiff
      << " This will be treated as ";
      
    if( nonPeriodic_[i] )
    {
      Xyce::dout() << " Non-periodic";
    }
    else
    {
      Xyce::dout() << " Periodic";
    }
    
    Xyce::dout() << std::endl;      
//#endif // Xyce_DEBUG_MPDE
  }
  return true;
}
*/

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::loadDAEMatrices( N_LAS_Vector * X,
                                     N_LAS_Vector * S,
                                     N_LAS_Vector * dSdt,
                                     N_LAS_Vector * Store,
                                     N_LAS_Matrix * dQdx,
                                     N_LAS_Matrix * dFdx)
{
#ifdef Xyce_DEBUG_MPDE
  if (mpdeMgrRCPtr_->debugLevel > 0)
  {
    Xyce::dout() << std::endl
           << Xyce::section_divider << std::endl
           << "  N_MPDE_Loader::loadDAEMatrices" << std::endl
           << "  warpMPDE flag = " << warpMPDE_ << std::endl;
  }
#endif // Xyce_DEBUG_MPDE
  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  N_LAS_Vector appdSdt( *appNextStaVecPtr_ );

  N_LAS_BlockMatrix & bdQdx = *dynamic_cast<N_LAS_BlockMatrix*>(dQdx);
  N_LAS_BlockMatrix & bdFdx = *dynamic_cast<N_LAS_BlockMatrix*>(dFdx);
  N_LAS_BlockVector & bX = *dynamic_cast<N_LAS_BlockVector*>(X);

#ifdef Xyce_FLEXIBLE_DAE_LOADS
  N_LAS_BlockVector & bS = *dynamic_cast<N_LAS_BlockVector*>(S);
  N_LAS_BlockVector & bdSdt = *dynamic_cast<N_LAS_BlockVector*>(dSdt);
  N_LAS_BlockVector & bStore = *dynamic_cast<N_LAS_BlockVector*>(Store);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  
  int BlockCount = bX.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
#ifdef Xyce_FLEXIBLE_DAE_LOADS
    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    mpdeDevInterfacePtr_->setFastTime( times_[i] );

    //Update the sources
    appLoaderPtr_->updateSources();

    *appNextVecPtr_ = bX.block(i);
    *appNextStaVecPtr_ = bS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);

    appLoaderPtr_->loadDAEMatrices( appNextVecPtr_, appNextStaVecPtr_, &appdSdt, 
        appNextStoVecPtr_, appdQdxPtr_, appdFdxPtr_);

    bdQdx.block(i,i).add( *appdQdxPtr_ );
    bdFdx.block(i,i).add( *appdFdxPtr_ );
#else
    //For now, the matrices are loaded during the loadDAEVectors method
    //Just copied here
    bdQdx.block(i,i).add( bmdQdxPtr_->block(i,i) );
    bdFdx.block(i,i).add( bmdFdxPtr_->block(i,i) );

#endif // Xyce_FLEXIBLE_DAE_LOADS
  }
  if (warpMPDE_)
  {
    // Add \dot{phi} = omega dQdx contribution here
    int phiGID = warpMPDEPhasePtr_->getPhiGID();
    int phiLID = bX.pmap()->globalToLocalIndex(phiGID);
    if (phiLID >= 0)
    {
      std::vector<int> colIndices;
      std::vector<double> coeffs;
      colIndices.push_back(phiGID);
      coeffs.push_back(1.0);
      bdQdx.replaceAugmentedRow( phiGID, colIndices.size(), &coeffs[0], &colIndices[0] );
    
      // tscoffe/tmei 08/10/05:  Add omega dependency for warped MPDE to every dq/dt2 equation
      // tscoffe 12/12/06:  For WaMPDE, we're augmenting two columns, omega & phi, and
      // the indexing starts at zero, so the index of omega is 0.
      // tscoffe 01/11/07:  Add \dot{phi} = omega dFdx derivative term here:
      (*bOmegadQdt2Ptr_)[phiLID] = -1.0;
    }
    bdFdx.replaceAugmentedColumn(warpMPDEPhasePtr_->getOmegaGID(),*bOmegadQdt2Ptr_);
  }

  // Fast Time scale discretization terms:
  // These are d(dQ/dt1)/dx terms, but go into bdFdx.  
  // For this procedure, need to re-use the app matrix, appdQdx.
  N_LAS_Matrix & dQdxdt = *appdQdxPtr_;

  int Start = fastTimeDiscRCPtr_->Start();
  int Width = fastTimeDiscRCPtr_->Width();
  
  const std::vector<double> & Coeffs = fastTimeDiscRCPtr_->Coeffs();
  
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Processing off diagonal matrix blocks on row " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    int Loc;
    int indexT1 = i + Start + periodicTimesOffset_;
    int indexT2 = indexT1 + Width - 1;
    double invh2 = 1.0 / (periodicTimes_[indexT2] - periodicTimes_[indexT1]);
    
    for( int j = 0; j < Width; ++j )
    {
      Loc = i + (j+Start);
      
      if( Loc < 0 )
      {
        Loc += BlockCount;
      }
      else if( Loc > (BlockCount-1) )
      {
        Loc -= BlockCount;
      }
      
      dQdxdt.put(0.0);
      dQdxdt.add( bdQdx.block(Loc,Loc) );
      dQdxdt.scale( Coeffs[j]*invh2 );
      bdFdx.block(i,Loc).add( dQdxdt );
    }
  }

  // tscoffe/tmei 08/10/05:  Add omega equation for warped MPDE
  if ( warpMPDE_ ) 
  {
    std::vector<int> colIndices;
    std::vector<double> coeffs;
    int omegaGID = warpMPDEPhasePtr_->getOmegaGID();
    warpMPDEPhasePtr_->getPhaseConditionDerivative(&bX,times_,&colIndices,&coeffs);
    bdFdx.replaceAugmentedRow( omegaGID, colIndices.size(), &coeffs[0], &colIndices[0] );
  }

  // Now that the matrix loading is finished, call fillComplete().
  dQdx->fillComplete();
  dFdx->fillComplete();

  // For BlockMatrix objects, synchronize the global copy of the block matrix.
  bdQdx.assembleGlobalMatrix();
  bdFdx.assembleGlobalMatrix();
 
#ifdef Xyce_DEBUG_MPDE
  if (mpdeMgrRCPtr_->debugLevel > 1)
  {
    Xyce::dout() << "MPDE bX:" << std::endl;
    bX.printPetraObject(std::cout);
    Xyce::dout() << "MPDE bdQdx:" << std::endl;
    bdQdx.printPetraObject(std::cout);
    Xyce::dout() << "MPDE bdFdx:" << std::endl;
    bdFdx.printPetraObject(std::cout);
#ifdef Xyce_FLEXIBLE_DAE_LOADS
    Xyce::dout() << "MPDE bS:" << std::endl;
    bS.printPetraObject(std::cout);
    Xyce::dout() << "MPDE dSdt:" << std::endl;
    bdSdt.printPetraObject(std::cout);
    Xyce::dout() << "MPDE bStore:" << std::endl;
    bStore.printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

#endif // Xyce_DEBUG_MPDE
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_MPDELoader::updateState
// Purpose       :
// Special Notes : ERK.  This function needs to be a no-op.  The reason
//                 is that the state information needs to be the same
//                 at the time of updateState, loadDAEVectors and 
//                 loadDAEMatrices.  Thus, they have to all happen inside
//                 of the same "fast time" loop.  So, this functionality
//                 has been moved up into loadDAEVectors.
//
//                 Note: for similar reasons, loadDAEMatrices is called from
//                 within that function as well.
//
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::updateState      
 (N_LAS_Vector * nextSolVectorPtr, 
  N_LAS_Vector * currSolVectorPtr,
  N_LAS_Vector * lastSolVectorPtr,
  N_LAS_Vector * nextStaVectorPtr,
  N_LAS_Vector * currStaVectorPtr,
  N_LAS_Vector * lastStaVectorPtr,
  N_LAS_Vector * nextStoVectorPtr,
  N_LAS_Vector * currStoVectorPtr,
  N_LAS_Vector * lastStoVectorPtr
  )
{
  bool bsuccess = true;

  // For MPDE case, this needs to be a no-op.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::loadDAEVectors( N_LAS_Vector * X,
                                    N_LAS_Vector * currX,
                                    N_LAS_Vector * lastX,
                                    N_LAS_Vector * S,
                                    N_LAS_Vector * currS,
                                    N_LAS_Vector * lastS,
                                    N_LAS_Vector * dSdt,
                                    N_LAS_Vector * Store,
                                    N_LAS_Vector * currStore,
                                    N_LAS_Vector * lastStore,
                                    N_LAS_Vector * storeLeadCurrQComp,
                                    N_LAS_Vector * Q,
                                    N_LAS_Vector * F,
                                    N_LAS_Vector * B,
                                    N_LAS_Vector * dFdxdVp,
                                    N_LAS_Vector * dQdxdVp )
{
#ifdef Xyce_DEBUG_MPDE
  if (mpdeMgrRCPtr_->debugLevel > 0)
  {
    Xyce::dout() << std::endl
           << Xyce::section_divider << std::endl
           << "  N_MPDE_Loader::loadDAEVectors" << std::endl
           << "warpMPDE flag = " << (warpMPDE_ ? "true" : "false") << std::endl;
  }
#endif // Xyce_DEBUG_MPDE
  //Zero out vectors
  appNextVecPtr_->putScalar(0.0);
  appCurrVecPtr_->putScalar(0.0);
  appLastVecPtr_->putScalar(0.0);

  appNextStaVecPtr_->putScalar(0.0);
  appCurrStaVecPtr_->putScalar(0.0);
  appLastStaVecPtr_->putScalar(0.0);
  N_LAS_Vector appdSdt( *appNextStaVecPtr_ );
  appNextStoVecPtr_->putScalar(0.0);
  appCurrStoVecPtr_->putScalar(0.0);
  appLastStoVecPtr_->putScalar(0.0);
  appStoLeadCurrQCompVecPtr_->putScalar(0.0);

  N_LAS_Vector appQ( *appNextVecPtr_ );
  N_LAS_Vector appF( *appNextVecPtr_ );
  N_LAS_Vector appB( *appNextVecPtr_ );
  N_LAS_Vector appdFdxdVp( *appNextVecPtr_ );
  N_LAS_Vector appdQdxdVp( *appNextVecPtr_ );

  // This is a temporary load storage vector.
  N_LAS_Vector dQdt2( *appNextVecPtr_ ); 

  // 12/8/06 tscoffe:   Note:  "b" at beginning of variable name means N_LAS_BlockVector
  N_LAS_BlockVector & bX = *dynamic_cast<N_LAS_BlockVector*>(X);
  N_LAS_BlockVector & bcurrX = *dynamic_cast<N_LAS_BlockVector*>(currX);
  N_LAS_BlockVector & blastX = *dynamic_cast<N_LAS_BlockVector*>(lastX);
  N_LAS_BlockVector & bS = *dynamic_cast<N_LAS_BlockVector*>(S);
  N_LAS_BlockVector & bcurrS = *dynamic_cast<N_LAS_BlockVector*>(currS);
  N_LAS_BlockVector & blastS = *dynamic_cast<N_LAS_BlockVector*>(lastS);
  N_LAS_BlockVector & bdSdt = *dynamic_cast<N_LAS_BlockVector*>(dSdt);
  N_LAS_BlockVector & bStore = *dynamic_cast<N_LAS_BlockVector*>(Store);
  N_LAS_BlockVector & bcurrStore = *dynamic_cast<N_LAS_BlockVector*>(currStore);
  N_LAS_BlockVector & blastStore = *dynamic_cast<N_LAS_BlockVector*>(lastStore);
  N_LAS_BlockVector & bQ = *dynamic_cast<N_LAS_BlockVector*>(Q);
  N_LAS_BlockVector & bF = *dynamic_cast<N_LAS_BlockVector*>(F);

  N_LAS_BlockVector & bdFdxdVp = *dynamic_cast<N_LAS_BlockVector*>(dFdxdVp);
  N_LAS_BlockVector & bdQdxdVp = *dynamic_cast<N_LAS_BlockVector*>(dQdxdVp);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
  bmdQdxPtr_->put(0.0);
  bmdFdxPtr_->put(0.0);
#endif
    
  int BlockCount = bQ.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Processing vectors for block " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    mpdeDevInterfacePtr_->setFastTime( times_[i] );
    
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Calling updateSources on the appLoader" << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    //Update the sources
    appLoaderPtr_->updateSources();  // this is here to handle "fast" sources.

    *appNextVecPtr_ = bX.block(i);
    *appCurrVecPtr_ = bcurrX.block(i);
    *appLastVecPtr_ = blastX.block(i);

    *appNextStaVecPtr_ = bS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    *appLastStoVecPtr_ = blastStore.block(i);
    *appStoLeadCurrQCompVecPtr_ = blastStore.block(i);  // need to update to correct block

#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Updating State for block " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE

    // Note: This updateState call is here (instead of in the 
    // N_MPDE_Loader::updateState function) because it has to be called
    // for the same fast time point.
    appLoaderPtr_->updateState 
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_,
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_ ,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_ , &*appLastStoVecPtr_ );

    bS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;
    blastStore.block(i) = *appLastStoVecPtr_;

#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Calling loadDAEVectors on the appLoader" << std::endl;
    }
#endif // Xyce_DEBUG_MPDE

    // This has to be done because the app loader does NOT zero these vectors out.
    appQ.putScalar(0.0);
    appF.putScalar(0.0);
    appB.putScalar(0.0);
    appdFdxdVp.putScalar(0.0);
    appdQdxdVp.putScalar(0.0);

    appLoaderPtr_->loadDAEVectors
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_, 
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_, &*appLastStaVecPtr_, 
        &appdSdt, 
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_, &*appLastStoVecPtr_, &*appStoLeadCurrQCompVecPtr_, 
        &appQ, &appF, &appB,
        &appdFdxdVp, &appdQdxdVp );

    bQ.block(i) = appQ;
    bF.block(i) = appF;
    bdFdxdVp.block(i) = appdFdxdVp;
    bdQdxdVp.block(i) = appdQdxdVp;

#ifndef Xyce_FLEXIBLE_DAE_LOADS
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Processing matrices for block " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE

    // This has to be done because the app loader does NOT zero these out.
    appdQdxPtr_->put(0.0);
    appdFdxPtr_->put(0.0);

    appLoaderPtr_->loadDAEMatrices( &*appNextVecPtr_, &*appNextStaVecPtr_, &appdSdt, &*appNextStoVecPtr_, 
                                    &*appdQdxPtr_, &*appdFdxPtr_);

#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Copying diagonal block into bmdQdx" << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    bmdQdxPtr_->block(i,i).add( *appdQdxPtr_ );
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Copying diagonal block into bmdFdx" << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    bmdFdxPtr_->block(i,i).add( *appdFdxPtr_ );
#endif
  }
  
  int phiGID = -1;
  int phiLID = -1;
  if (warpMPDE_)
  {
    // Add \dot{phi(t_1)} = omega(t_1) term to Q
    phiGID = warpMPDEPhasePtr_->getPhiGID();
    phiLID = bQ.pmap()->globalToLocalIndex(phiGID);
    if (phiLID >= 0)
    {
      double phiValue = bX[phiLID];

#ifdef Xyce_DEBUG_MPDE
      if (mpdeMgrRCPtr_->debugLevel > 0)
      {
        Xyce::dout() << "Inserting phi with value = " << phiValue << " into q" << std::endl;
      }
#endif // Xyce_DEBUG_MPDE

      bQ.setElementByGlobalIndex( phiGID, phiValue );
    }
  }

  //-------------------------------------------------------------
  // Now do the fast time scale time derivative -----------------
  //-------------------------------------------------------------

  int Start = fastTimeDiscRCPtr_->Start();
  int Width = fastTimeDiscRCPtr_->Width();
  const std::vector<double> & Coeffs = fastTimeDiscRCPtr_->Coeffs();
    
  double omega = 1.0;
  int omegaGID = -1;
  if (warpMPDE_)
  {
    omegaGID = warpMPDEPhasePtr_->getOmegaGID();

    // tscoffe/tmei 08/04/05:  zero out temporary vector for use in Jacobian
    bOmegadQdt2Ptr_->putScalar(0.0);
    // tscoffe/tmei 08/11/05:  Convert global indices to local indices:

    int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID);

    double tmpOmega = 0.0; 
    if (omegaLID >= 0)
    {
      tmpOmega = bX[omegaLID];
    }
    bX.pmap()->pdsComm()->sumAll( &tmpOmega, &omega, 1 );
  }

  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_MPDE
    if (mpdeMgrRCPtr_->debugLevel > 0)
    {
      Xyce::dout() << "Processing dQdt2 block " << i << " of " << BlockCount-1 << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
    dQdt2.putScalar(0.0);
    
    int Loc;
    int indexT1 = i + Start + periodicTimesOffset_;
    int indexT2 = indexT1 + Width - 1;
    double invh2 = 1.0 / (periodicTimes_[indexT2] - periodicTimes_[indexT1]);
    
    for( int j = 0; j < Width; ++j )
    {
      Loc = i + (j+Start);
      
      if( Loc < 0 )
      {
        Loc += BlockCount;
      }
      else if( Loc > (BlockCount-1) )
      {
        Loc -= BlockCount;
      }
      dQdt2.update( Coeffs[j]*invh2, bQ.block(Loc), 1.0 );
    }
    
    if (warpMPDE_)
    {
      // 12/8/06 tscoffe:  save dQdt2 term for Jacobian evaluation
      bOmegadQdt2Ptr_->block(i) = dQdt2; // This is not scaled by omega
    }
    // Update F with omega*dq/dt2 term
    bF.block(i).update( omega, dQdt2, 1.0 );
  }

  //  tscoffe/tmei 08/04/05:  Add omega equation to end of f:
  if ( warpMPDE_ )
  {
    double phaseValue = warpMPDEPhasePtr_->getPhaseCondition(&bX,times_);
    if ( bF.pmap()->globalToLocalIndex(omegaGID) >= 0 )
    {
#ifdef Xyce_DEBUG_MPDE
      if (mpdeMgrRCPtr_->debugLevel > 0)
      {
        Xyce::dout() << "Inserting phase condition = " << phaseValue << " into f" << std::endl;
        Xyce::dout() << "Inserting omega = " << omega << " into f" << std::endl;
      }
#endif // Xyce_DEBUG_MPDE

      bF.setElementByGlobalIndex( omegaGID, phaseValue );
      // Add \dot{phi(t_1)} = omega(t_1) term to F
      bF.setElementByGlobalIndex( phiGID, -omega );
    }
  }

  // Now that the vector loading is finished, synchronize the global copy of the block vector
  bX.assembleGlobalVector();
  bS.assembleGlobalVector();
  bdSdt.assembleGlobalVector();
  bStore.assembleGlobalVector();
  bQ.assembleGlobalVector();
  bF.assembleGlobalVector();
  bdFdxdVp.assembleGlobalVector();
  bdQdxdVp.assembleGlobalVector();
  bcurrS.assembleGlobalVector();
  blastS.assembleGlobalVector();
  bcurrStore.assembleGlobalVector();
  blastStore.assembleGlobalVector();

#ifdef Xyce_DEBUG_MPDE
  if (mpdeMgrRCPtr_->debugLevel > 1)
  {
    Xyce::dout() << "MPDE X Vector" << std::endl;
    bX.printPetraObject(std::cout);
    Xyce::dout() << "MPDE S Vector" << std::endl;
    bS.printPetraObject(std::cout);
    Xyce::dout() << "MPDE dSdt Vector" << std::endl;
    bdSdt.printPetraObject(std::cout);
    Xyce::dout() << "MPDE Store Vector" << std::endl;
    bStore.printPetraObject(std::cout);
    Xyce::dout() << "MPDE Q Vector" << std::endl;
    bQ.printPetraObject(std::cout);
    Xyce::dout() << "MPDE F Vector" << std::endl;
    bF.printPetraObject(std::cout);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
    bmdQdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "MPDE bmdQdx_" << std::endl;
    bmdQdxPtr_->printPetraObject(std::cout);

    bmdFdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "MPDE bmdFdx_" << std::endl;
    bmdFdxPtr_->printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS

    Xyce::dout() << Xyce::section_divider << std::endl;
  }
#endif // Xyce_DEBUG_MPDE

  return true;
}

