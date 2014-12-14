//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTInterfaceDecl.hpp,v $
//
// Purpose        : This class acts as an abstract interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 11/11/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2014/04/21 21:59:35 $
//
// Current Owner  : $Author: tmei $
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_DECL_HPP
#define N_UTL_FFTINTERFACE_DECL_HPP


// ---------- Standard Includes ----------

#include <N_UTL_DFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterfaceDecl
// Purpose       : This class acts as an abstract interface to an FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  
// Special Notes :
// Creator       : Heidi Thornquist
// Creation Date : 11/11/10
// Last Modified : 11/11/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTInterfaceDecl : public N_UTL_DFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_FFTInterfaceDecl( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      { signalLength_ = length; numberSignals_ = numSignals; stride_ = reqStride; overwrite_ = overwrite; }
    
    // Basic destructor 
    virtual ~N_UTL_FFTInterfaceDecl() {};

    // Return signal length.
    virtual int getSignalLength()
    { return signalLength_; }

    virtual int getScalar()
    { return signalLength_; }

  protected:
    // this is the length of the real array.  The complex result will be signalLength+2 long
    // if signalLength is even and signalLength+1 if signalLength is odd
    int signalLength_;

    // this is the number of signals on which we will take an fft/ift
    int numberSignals_;

    // If the signals are grouped by blocks at the same time (say x0, x1 ... xn at t0) and
    // then (x0, x1 ... xn at t1). Then stride is the spacing from one x0 at t0 to the next
    // x0 at t1.  This lets one take ffts/ifts of data that is blocked by time
    int stride_;

    // Whether the input and output vectors should be expected to be the same, save space if possible.
    bool overwrite_;
};

#endif

