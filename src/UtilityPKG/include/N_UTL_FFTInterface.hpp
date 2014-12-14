//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTInterface.hpp,v $
//
// Purpose        : This class acts as an interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  
//
// Special Notes  : 
//
// Creator        : Richard Schiek 
//
// Creation Date  : 5/27/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/04/08 16:24:05 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_HPP
#define N_UTL_FFTINTERFACE_HPP 1


// ---------- Standard Includes ----------
#include <Xyce_config.h>

#ifdef Xyce_USE_INTEL_FFT
#include <N_UTL_IntelFFT_Interface.hpp>
#endif

#ifdef Xyce_USE_FFTW
#include <N_UTL_FFTW_Interface.hpp>
#endif

#include <N_UTL_FFTInterfaceDecl.hpp>

#include <N_ERH_ErrorMgr.h>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface
// Purpose       : This class acts as a templated interface to any FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  It is originally
//                 implemented for Intel's Math Library but may be extended
//                 to FFTW at some time in the future.
// Special Notes :
// Creator       : Richard Schiek (templated by Heidi Thornquist)
// Creation Date : 5/27/08
// Last Modified : 11/17/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTInterface
{
  public:
    N_UTL_FFTInterface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
    {
#ifdef Xyce_USE_INTEL_FFT
      fftInterface_ = Teuchos::rcp( new N_UTL_IntelFFT_Interface<VectorType>( length, numSignals, reqStride, overwrite ) );
#elif defined(Xyce_USE_FFTW)
      fftInterface_ = Teuchos::rcp( new N_UTL_FFTW_Interface<VectorType>( length, numSignals, reqStride, overwrite ) );
#else
      std::string msg = "Xyce has not been configured with an FFT library! Please reconfigure with FFT enabled to perform any frequency analysis!\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
#endif
    }
    
    virtual ~N_UTL_FFTInterface() {}
   
    // Register new vectors for the FFT/IFT interface to use.
    void registerVectors( VectorType& fftInData, VectorType* fftOutData,
                          VectorType& iftInData, VectorType* iftOutData )
    {
      fftInterface_->registerVectors( Teuchos::rcp( &fftInData, false ), Teuchos::rcp( fftOutData, false ),
                                      Teuchos::rcp( &iftInData, false ), Teuchos::rcp( iftOutData, false ) );
    }

    // Return the vectors that were registered for the FFT interface to use.
    void getFFTVectors( Teuchos::RCP<VectorType>& fftInData,  Teuchos::RCP<VectorType>& fftOutData )
    {
      fftInterface_->getDFTVectors( fftInData, fftOutData );
    }

    // Return the vectors that were registered for the IFT interface to use.
    // NOTE:  Consider using the already registered vectors to avoid recreation of the FFT object.
    void getIFTVectors( Teuchos::RCP<VectorType>& iftInData,  Teuchos::RCP<VectorType>& iftOutData )
    {
      fftInterface_->getIFTVectors( iftInData, iftOutData );
    }
 
    void calculateFFT( VectorType& inData, VectorType* outResult)
    {
      fftInterface_->calculateDFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
    }
    
    void calculateIFT( VectorType& inData, VectorType* outResult)
    {
      fftInterface_->calculateIFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
    }
   
    void calculateFFT()
    {
      fftInterface_->calculateDFT();
    }
    
    void calculateIFT()
    {
      fftInterface_->calculateIFT();
    }

    Teuchos::RCP< N_UTL_FFTInterfaceDecl<VectorType> > getFFTInterface()
    {
      return fftInterface_;
    }

  private:

  // A pointer to the FFT interface used by this class to compute the forward and backward transforms.
  Teuchos::RCP< N_UTL_FFTInterfaceDecl<VectorType> > fftInterface_;
};

#endif

