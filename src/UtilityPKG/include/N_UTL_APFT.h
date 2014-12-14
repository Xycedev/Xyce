//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_APFT.h,v $
//
// Purpose        : This class performs DFT and IFT for APFT algorithm
//
// Special Notes  : 

//
// Creator        : Ting Mei
//
// Creation Date  : 4/09/14
//
// Revision Information:
//-------------------------------------------------------------------------
#ifndef N_UTL_APFT_H
#define N_UTL_APFT_H


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : N_UTL_APFT
// Purpose       : This class performs DFT and IFT 
// Special Notes :
// Creator       : Ting Mei
// Creation Date : 4/15/14
// Last Modified : 4/15/14
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_APFT: public N_UTL_DFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.

    N_UTL_APFT(const Teuchos::SerialDenseMatrix<int,double>& idftMatrix, const Teuchos::SerialDenseMatrix<int,double>& dftMatrix )
    {

      idftMatrix_ = idftMatrix;
      dftMatrix_ = dftMatrix;
     
//      std::cout << "checking the IDFT matrix in FT interface"  << std::endl;
//      idftMatrix_.print(std::cout);

//      std::cout << "checking the dftmatrix in FT interface"  << std::endl;
//      dftMatrix_.print(std::cout); 

    }

    // Basic destructor 
    virtual ~N_UTL_APFT() {}

    // Calculate FFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateDFT();
    // Calculate IFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateIFT();

    void calculateDFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      N_UTL_DFTInterfaceDecl<VectorType>::calculateDFT( inData, outData );
    }

    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateIFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      N_UTL_DFTInterfaceDecl<VectorType>::calculateIFT( inData, outData );
    }

    int getScalar()
    {return 1;}

  private:

  // Fourier matrices
  Teuchos::SerialDenseMatrix<int,double> idftMatrix_, dftMatrix_;

};

#endif
