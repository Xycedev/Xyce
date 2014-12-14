//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_DFTInterfaceDecl.hpp,v $
//
// Purpose        : This class acts as an abstract interface to an DFT library
//                  for DFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given DFT library so 
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/07/22 19:43:04 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------
#ifndef N_UTL_DFTINTERFACE_DECL_HPP
#define N_UTL_DFTINTERFACE_DECL_HPP


// ---------- Standard Includes ----------

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : DFTInterfaceDecl
// Purpose       : This class acts as an abstract interface to an DFT library
//                 for DFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given DFT library so 
//                 that multiple libraries and/or specialized implementations
//                 can be used.  
// Special Notes :
// Creator       : Heidi Thornquist
// Creation Date : 4/7/14
// Last Modified : 4/7/14
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_DFTInterfaceDecl
{
  public:

    // Basic destructor 
    virtual ~N_UTL_DFTInterfaceDecl() {};

    // Register new vectors for the DFT/IFT interface to use.
    virtual void registerVectors( const Teuchos::RCP<VectorType>& dftInData, const Teuchos::RCP<VectorType>& dftOutData,
                                  const Teuchos::RCP<VectorType>& iftInData, const Teuchos::RCP<VectorType>& iftOutData )
    { dftInData_ = dftInData; dftOutData_ = dftOutData; iftInData_ = iftInData; iftOutData_ = iftOutData; }
  
    // Return the vectors that were registered for the DFT interface to use.
    // NOTE:  Consider using the already registered vectors to avoid recreation of the DFT object. 
    virtual void getDFTVectors( Teuchos::RCP<VectorType>& dftInData,  Teuchos::RCP<VectorType>& dftOutData )
    {  dftInData = dftInData_; dftOutData = dftOutData_; }
 
    // Return the vectors that were registered for the IFT interface to use.
    // NOTE:  Consider using the already registered vectors to avoid recreation of the DFT object. 
    virtual void getIFTVectors( Teuchos::RCP<VectorType>& iftInData,  Teuchos::RCP<VectorType>& iftOutData )
    {  iftInData = iftInData_; iftOutData = iftOutData_; }
 
    // Calculate DFT with new vectors, not the ones that have been registered or used in the constructor.
    virtual void calculateDFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // In some cases, it doesn't matter if the vectors have changed for the DFT library
      dftInData_ = inData;
      dftOutData_ = outData;
      this->calculateDFT();
    }

    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    virtual void calculateIFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // In some cases, it doesn't matter if the vectors have changed for the DFT library
      iftInData_ = inData;
      iftOutData_ = outData;
      this->calculateIFT();
    }

    // Calculate DFT with the vectors that have been registered.
    virtual void calculateDFT() = 0;
    // Calculate IFT with the vectors that have been registered.
    virtual void calculateIFT() = 0;

    virtual int getScalar() = 0;

  protected:

    // DFT/IFT vectors
    Teuchos::RCP<VectorType> dftInData_, iftInData_;
    Teuchos::RCP<VectorType> dftOutData_, iftOutData_;

};

#endif

