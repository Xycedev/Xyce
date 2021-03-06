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
// Filename       : $RCSfile: N_IO_Objective.C,v $
//
// Purpose        : Objective class, for Dakota
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.38 $
//
// Revision Date  : $Date: 2014/07/15 19:03:34 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>
#include <fstream>
#include <string>

#include <sstream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <Epetra_SerialDenseVector.h>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

//------- Xyce Includes ----------
#include <N_IO_Objective.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>

#include <N_TOP_Topology.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Version.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_FFTInterface.hpp>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : Objective::Objective
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 02/13/06
//-----------------------------------------------------------------------------
Objective::Objective()
  : value(NULL),
    weight(NULL),
    valueSetup(false),
    lastResult(0),
    minResult(0),
    maxResult(0),
    magResult(0),
    rmsResult(0),
    lastV1(0.0),
    lastWeight(0.0),
    lastInd2(0),
    lastIndInterpolate(0),
    lastN1(0),
    lastN2(0),
    n1(0),
    n2(0),
    lastResultValid(false)
{}

Objective::Objective(const Objective &objective)
  : value(objective.value),
    weight(objective.weight),
    valueSetup(false),
    lastResult(objective.lastResult),
    minResult(objective.minResult),
    maxResult(objective.maxResult),
    magResult(objective.magResult),
    rmsResult(objective.rmsResult),
    lastV1(objective.lastV1),
    lastWeight(objective.lastWeight),
    lastInd2(objective.lastInd2),
    lastIndInterpolate(objective.lastIndInterpolate),
    lastN1(objective.lastN1),
    lastN2(objective.lastN2),
    n1(objective.n1),
    n2(objective.n2),
    lastResultValid(objective.lastResultValid)
{}

//-----------------------------------------------------------------------------
// Function      : Objective::~Objective
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 02/13/06
//-----------------------------------------------------------------------------
Objective::~Objective()
{
  if (value != NULL)
    delete value;
  if (weight != NULL)
    delete weight;
}

//-----------------------------------------------------------------------------
// Function      : Objective:parse
// Purpose       : Parse input string, and set objective
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 04/06/06
//-----------------------------------------------------------------------------

bool Objective::parse(const std::string & in)
{
  std::string line (in);
  std::string fileStr, functionStr, valueStr, weightStr;
  std::string::size_type n, begin, end;

  if (line.find_first_of('=') == std::string::npos)
    valueStr = line;
  else
  {
    while ((n = line.find_first_of(' ')) != std::string::npos)
      line.erase(n,1);
    while ((n = line.find_first_of('\t')) != std::string::npos)
      line.erase(n,1);
    ExtendedString uLine(line);
    uLine.toUpper();

    // get the filename
    n = uLine.find(std::string("FILE="));
    if (n != std::string::npos)
    {
      begin = n + 6;
      if (line.size() < begin+1 || line[begin-1] != '"')
        return false;
      end = begin + uLine.substr(begin).find_first_of('"');
      if (end == std::string::npos)
        return false;
      fileStr = line.substr(begin,end-begin);
    }

    // get the function we will apply
    n = uLine.find(std::string("FUNCTION="));
    if (n != std::string::npos) {
      begin = n + 10;
      if (line.size() < begin+1 || line[begin-1] != '"')
        return false;
      end = begin + uLine.substr(begin).find_first_of('"');
      if (end == std::string::npos)
        return false;
      functionStr = line.substr(begin,end-begin);
    }

    // get the value expression
    n = uLine.find(std::string("VALUE="));
    if (n != std::string::npos)
    {
      begin = n + 7;
      if (line.size() < begin+1 || line[begin-1] != '{')
        return false;
      end = begin + uLine.substr(begin).find_first_of('}');
      if (end == std::string::npos)
        return false;
      valueStr = line.substr(begin,end-begin);
    }

    // get the weighting expression
    n = uLine.find(std::string("WEIGHT="));
    if (n != std::string::npos)
    {
      begin = n + 8;
      if (line.size() < begin+1 || line[begin-1] != '{')
        return false;
      end = begin + uLine.substr(begin).find_first_of('}');
      if (end == std::string::npos)
        return false;
      weightStr = line.substr(begin,end-begin);
    }
  }
  return initializeInternal(fileStr, functionStr, valueStr, weightStr);
}

//-----------------------------------------------------------------------------
// Function      : Objective:initialize
// Purpose       : Parse input string, and set objective
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 04/07/06
//-----------------------------------------------------------------------------
bool Objective::initialize(const Util::OptionBlock &option_block)
{
  std::string fileStr, functionStr, valueStr, weightStr;
  for (std::list<Util::Param>::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end(); it != end; ++it)
  {
    if ((*it).uTag() == "NAME")
    {
    }
    else if ((*it).uTag() == "FUNCTION" )
    {
      functionStr = (*it).stringValue();
    }
    else if ((*it).uTag() == "VALUE")
    {
      valueStr = (*it).stringValue();
    }
    else if ((*it).uTag() == "WEIGHT")
    {
      weightStr = (*it).stringValue();
    }
    else if ((*it).uTag() == "FILE")
    {
      fileStr = (*it).stringValue();
    }
  }
  
  return initializeInternal(fileStr, functionStr, valueStr, weightStr);
}

//-----------------------------------------------------------------------------
// Function      : Objective:initializeInternal
// Purpose       : Internal intialization from .objective or objective name
//                 from Dakota
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley
// Creation Date : 04/07/06
//-----------------------------------------------------------------------------
bool
Objective::initializeInternal(
  std::string & fileStr,
  std::string & functionStr,
  std::string & valueStr,
  std::string & weightStr)
{
  if (fileStr.size() > 0)
  {
    file = fileStr;
    readData();
  }
  if (functionStr.size() == 0)
  {
    // apply default behavior
    function = "VALUE";
  }
  else if (functionStr.size() > 0)
  {
    if( functionStr == "VALUE" ||
        functionStr == "L1NORM" ||
        functionStr == "L2NORM" ||
        functionStr == "L1NORM_POWERSP" ||
        functionStr == "L2NORM_POWERSP" ||
        functionStr == "INFNORM_POWERSP" ||
        functionStr == "INFNORM" ||
        functionStr == "MAX" ||
        functionStr == "MIN" ||
        functionStr == "MAG" ||
        functionStr == "RMS" ||
        functionStr == "NMAX" ||
        functionStr == "NMIN" ||
        functionStr == "NMAG" ||
        functionStr == "NRMS" )
    {
      // found a valid choice
      function = functionStr;
    }
    else
    {
      // issue a fatal error
      std::string msg("Objective: Unknows function type requested \"" + functionStr + "\"");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL,msg);
    }
  }
  if (valueStr.size() > 0)
  {
    value = new Util::ExpressionData(valueStr);
  }
  if (weightStr.size() > 0)
  {
    weight = new Util::ExpressionData(weightStr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Objective::readData
// Purpose       : Reads in data from external file.  Needs improvement in
//                 getting a small amount of feedback to user that the operation
//                 was or was not successful. -- RLS
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void Objective::readData()
{
  int i;
  double v;

  if (file == "")
    return;

  std::ifstream inputFileStream;
  inputFileStream.open( file.c_str() );

  if( !inputFileStream )
  {
    std::string msg("Could not open file of experimental data: " + file);
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL,msg);
  }
  std::string inputLine;

  getline( inputFileStream, inputLine );
  std::stringstream line1;
  line1 << inputLine;
  line1 >> var1;
  getline( inputFileStream, inputLine );
  std::stringstream line2;
  line2 << inputLine;
  line2 >> var2;
  while (line2 >> v)
  {
    //cout << "Val = " << v << endl;
    var2Vals.push_back(v);
  }
  n2 = var2Vals.size();
  if (n2 <= 1)
    var2 = "";
  dataVal.resize(n2);
  while (getline( inputFileStream, inputLine ))
  {
    std::stringstream converterStream;
    converterStream << inputLine;
    converterStream >> v;
    var1Vals.push_back(v);
    for (i=0 ; i<n2 ; ++i)
    {
      converterStream >> v;
      dataVal[i].push_back(v);
    }
  }
  n1 = var1Vals.size();
  simVal.resize(n2);
  simValValid.resize(n2);
  for (i=0 ; i<n2 ; ++i)
  {
    simVal[i].resize(n1,0);
    simValValid[i].resize(n1,false);
  }
  weightVal.resize(n2);
  for (i=0 ; i<n2 ; ++i)
  {
    weightVal[i].resize(n1,0);
  }

  //printData();

  return;
}

//-----------------------------------------------------------------------------
// Function      : Objective::printData
// Purpose       : Debug function to print out read and stored data
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Sandia National Lab
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
void Objective::printData()
{
  int len = var1Vals.size();
  for(int i=0; i<len; i++ )
  {
    Xyce::dout() << "var1Vals[ " << i << " ] = " << var1Vals[i] << std::endl;
  }

  len = var2Vals.size();
  for(int i=0; i<len; i++ )
  {
    Xyce::dout() << "var2Vals[ " << i << " ] = " << var2Vals[i] << std::endl;
  }

  len = dataVal.size();
  for(int i=0; i<len; i++)
  {
    int len2 = dataVal[i].size();
    Xyce::dout() << "dataVal[ " << i << " ] = { ";
    for(int j=0; j<len2; j++)
    {
      Xyce::dout() << "   " << dataVal[i][j];
    }
    Xyce::dout() << " }" << std::endl;
  }

  len = simVal.size();
  for(int i=0; i<len; i++)
  {
    int len2 = simVal[i].size();
    Xyce::dout() << "simVal[ " << i << " ] = { ";
    for(int j=0; j<len2; j++)
    {
      Xyce::dout() << "   " << simVal[i][j];
    }
    Xyce::dout() << " }" << std::endl;
  }

  len = simValValid.size();
  for(int i=0; i<len; i++)
  {
    int len2 = simValValid[i].size();
    Xyce::dout() << "simValValid[ " << i << " ] = { ";
    for(int j=0; j<len2; j++)
    {
      Xyce::dout() << "   " << simValValid[i][j];
    }
    Xyce::dout() << " }" << std::endl;
  }

  len = weightVal.size();
  for(int i=0; i<len; i++)
  {
    int len2 = weightVal[i].size();
    Xyce::dout() << "weightVal[ " << i << " ] = { ";
    for(int j=0; j<len2; j++)
    {
      Xyce::dout() << "   " << weightVal[i][j];
    }
    Xyce::dout() << " }" << std::endl;
  }

}
//-----------------------------------------------------------------------------
// Function      : Objective::reset
// Purpose       : Reset objective to initial state
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------
void Objective::reset()
{
  int i, j;
  lastResultValid = false;

  for (i=0 ; i<n2 ; ++i)
  {
    for (j=0 ; j<n1 ; ++j)
    {
      simValValid[i][j] = false;
    }
  }

  lastResultValid = false;
  maxResult = 0.0;
  minResult = 0.0;
  magResult = 0.0;
  rmsResult = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : Objective::save
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
double Objective::save(Parallel::Machine comm, double current_circuit_time, const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  if (value)
    lastResult = value->evaluate(comm, current_circuit_time, solnVecPtr, stateVecPtr, storeVecPtr);
  else
    lastResult = 0.0;

  if( lastResultValid == false )
  {
    // first time through so save lastResult as our max and min
    maxResult = lastResult;
    minResult = lastResult;

    // set flag to indicate that we've sampled lastResult at least once
    lastResultValid = true;
  }

  // track max, min and mag
  if( lastResult > maxResult )
  {
    maxResult = lastResult;
    magResult = maxResult - minResult;
  }
  if( lastResult < minResult )
  {
    minResult = lastResult;
    magResult = maxResult - minResult;
  }

  return lastResult;
}

void
Objective::setup(Parallel::Machine comm, const OutputMgr &output_manager)
{
  if (!valueSetup) {
    valueSetup = true;

    if (value)
      value->setup(comm, output_manager);
    if (weight)
      weight->setup(comm, output_manager);
  }
}


//-----------------------------------------------------------------------------
// Function      : Objective::save
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
double
Objective::save(
  Parallel::Machine     comm,
  double                current_circuit_time,
  double                v1,
  double                v2,
  const N_LAS_Vector *  solnVecPtr,
  const N_LAS_Vector *  stateVecPtr,
  const N_LAS_Vector *  storeVecPtr)
{
  int i, j, ind1, ind2, indInterpolate;
  double val, wgt, sum;

  if (value)
    val = value->evaluate(comm, current_circuit_time, solnVecPtr, stateVecPtr, storeVecPtr);
  else
    val = 0;
  if (weight)
    wgt = weight->evaluate(comm, current_circuit_time, solnVecPtr, stateVecPtr, storeVecPtr);
  else
    wgt = 1;

   //Xyce::dout() << "Objective::save() ------------------------" << std::endl
   //    << "\t" << var1 << "\tv1 = " << v1 << std::endl
   //    << "\t" << var2 << "\tv2 = " << v2 << std::endl
   //    << "\tval = " << val << std::endl
   //    << "\twgt = " << wgt << std::endl;

  ind1 = -1;
  ind2 = -1;
  indInterpolate = -1;
  if (var1 != "")
  {
    for (i=lastN1 ; i<n1 ; ++i)
    {
      if (i<n1-1 && (v1-var1Vals[i])*(v1-var1Vals[i+1]) <= 0)
        indInterpolate = i;
      if (fabs(v1-var1Vals[i]) < 0.000001*(fabs(var1Vals[n1-1]-var1Vals[0])/n1))
      {
        ind1 = i;
        lastN1 = i;
        break;
      }
    }
    if( ind1 < 0 )
    {
      // reset lastN1 and try again
      // we use lastN1 to speed saves on big data sets, but
      // we don't have a good way to reset it when starting a new set
      lastN1=0;
      for (i=lastN1 ; i<n1 ; ++i)
      {
        if (i<n1-1 && (v1-var1Vals[i])*(v1-var1Vals[i+1]) <= 0)
          indInterpolate = i;
        if (fabs(v1-var1Vals[i]) < 0.000001*(fabs(var1Vals[n1-1]-var1Vals[0])/n1))
        {
          ind1 = i;
          lastN1 = i;
          break;
        }
      }
    }
  }
  if (var2 != "")
  {
    for (i=0 ; i<n2 ; ++i)
    {
      if (fabs(v2-var2Vals[i]) < 0.000001*(fabs(var2Vals[n2-1]-var2Vals[0])/n2))
      {
        ind2 = i;
        break;
      }
    }
  }
  else
  {
    ind2 = 0;
  }
  if (ind2 >= 0)
  {
    if (ind1 != -1 && ind2 != -1)
    {
      simVal[ind2][ind1] = val;
      simValValid[ind2][ind1] = true;
      weightVal[ind2][ind1] = wgt;
      indInterpolate = ind1;
    }
    else if (lastResultValid && lastInd2 == ind2)
    {
      if (lastIndInterpolate != indInterpolate)
      {
        int lo, hi;
        if (lastIndInterpolate < indInterpolate)
        {
          lo = lastIndInterpolate;
          hi = indInterpolate;
        }
        else
        {
          lo = indInterpolate;
          hi = lastIndInterpolate;
        }
        for (i=lo+1 ; i<= hi ; ++i)
        {
          double frac = (var1Vals[i] - lastV1)/(v1 - lastV1);
          simVal[ind2][i] = val*frac+ lastResult*(1-frac);
          simValValid[ind2][i] = true;
          weightVal[ind2][i] = wgt*frac+ lastWeight*(1-frac);
        }
      }
    }

    lastResultValid = true;
    lastResult = val;
    lastV1 = v1;
    lastWeight = wgt;
    lastInd2 = ind2;
    lastIndInterpolate = indInterpolate;
  }
  else
  {
    Xyce::dout() << "Objective::save() unable to save data:" << std::endl;
    // Being here indicates that we have data we are unable to save
  }

  return val;
}

//-----------------------------------------------------------------------------
// Function      : Objective::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
double Objective::evaluate() const
{
  double result = 0.0;

  //Xyce::dout() << "Objective::evaluate---------------------------------------" << std::endl;
  //printData();
  // Xyce::dout() << "In Objective::evaluate function \"" << function << "\"" << std::endl;
  if( function == "VALUE" )
  {
    // this is the default behavior.  Just
    // return the last result
    result = lastResult;
  }
  else if( function == "L1NORM" ||
           function == "L2NORM" ||
           function == "INFNORM" )
  {
    // make an epetra vector to hold difference values
    Epetra_SerialDenseVector differenceVector(n1*n2);
    int numValid = 0;
    // load the difference vector
    for (int i=0 ; i<n2 ; ++i)
    {
      for (int j=0 ; j<n1 ; ++j)
      {
        if (simValValid[i][j])
        {
          numValid++;
          differenceVector[i+j] = weightVal[i][j]*(simVal[i][j]-dataVal[i][j]);
        }
      }
    }

    // should issue a warning if numValid is much less than n1 or n1*n2
    // Xyce::dout() << "Used " << numValid << " points from a set of " << n1*n2 << " L2Norm = "<< result << std::endl;
    // printData();
    if( function == "L1NORM" )
    {
      result = differenceVector.Norm1();
    }
    else if( function == "L2NORM" )
    {
      result = differenceVector.Norm2();
    }
    else
    {
      result = differenceVector.NormInf();
    }

    /*
    int i, j, ind1, ind2;
    double val, wgt, sum, num, norm;

    if (var1.empty() && var2.empty())
    {
      return lastResult;
    }

    sum = 0;
    num = 0;
    norm = 0;

    for (i=0 ; i<n2 ; ++i)
    {
      for (j=0 ; j<n1 ; ++j)
      {
        if (simValValid[i][j])
        {
          sum += weightVal[i][j]*(simVal[i][j]-dataVal[i][j])*(simVal[i][j]-dataVal[i][j]);
          norm += (weightVal[i][j]*dataVal[i][j])*(weightVal[i][j]*dataVal[i][j]);
          num ++;
        }
      }
    }
    if (num > 0 && norm > 0)
    {
      sum = sum/norm;
    }
    else
    {
      sum = 0;
    }
    result = sum;
    */
  }
  else if( function == "L1NORM_POWERSP" ||
           function == "L2NORM_POWERSP" ||
           function == "INFNORM_POWERSP" )
  {
    // in this case we will calculate the power spectra of the two signals
    // and then take the appropriate norm of the differnce
    std::vector<double> simValPowerSpecResult;

    std::vector<double> dataValPowerSpecResult;

    // set up the fft engine to calculate the fft in place to conserve memory
    int sigLen = var1Vals.size();
    N_UTL_FFTInterface<std::vector<double> > fftEngine( sigLen, 1, 0, true );

    // copy simulation and external data to working vectors
    // for the calculation.  Note: we do filtering and weighting here as needed
    int numSig = simVal.size();
    std::vector<double> simValFFTResult;
    std::vector<double> dataValFFTResult;
    simValFFTResult.resize(2 * sigLen);
    dataValFFTResult.resize(2 * sigLen );

    for (int i=0 ; i<sigLen; i++)
    {
      for (int j=0 ; j<numSig; j++)
      {
        if (simValValid[j][i])
        {
          simValFFTResult[i] = weightVal[j][i]*simVal[j][i];
          dataValFFTResult[i] = weightVal[j][i]* dataVal[j][i];
        }
      }
    }

    // calculate the FFT
    fftEngine.calculateFFT( dataValFFTResult, &dataValFFTResult);
    fftEngine.calculateFFT( simValFFTResult, &simValFFTResult);

    // calculate Power spectra and difference vector
    Epetra_SerialDenseVector differenceVector(sigLen);
    for (int i=0, j=0 ; i<sigLen; i++, j+=2)
    {
      result = simValFFTResult[j]*simValFFTResult[j] + simValFFTResult[j+1]*simValFFTResult[j+1];
      simValFFTResult[i] = result;
      result = dataValFFTResult[j]*dataValFFTResult[j] + dataValFFTResult[j+1]*dataValFFTResult[j+1];
      dataValFFTResult[i] = result;
      differenceVector[i] = simValFFTResult[i] - dataValFFTResult[i];
    }

    // now calculate the appropriate norm
    if( function == "L1NORM_POWERSP" )
    {
      result = differenceVector.Norm1();
    }
    else if( function == "L2NORM_POWERSP" )
    {
      result = differenceVector.Norm2();
    }
    else
    {
      result = differenceVector.NormInf();
    }

    // when we go out of scope, the fftEngine will be de-allocated
  }
  else if( function == "MAX" )
  {
    result = maxResult;
  }
  else if( function == "MIN" )
  {
    result = minResult;
  }
  else if( function == "MAG" )
  {
    result = magResult;
  }
  else if( function == "RMS" )
  {
    result = rmsResult;
  }
  else if( function == "NMAX" )
  {
    result = -maxResult;
  }
  else if( function == "NMIN" )
  {
    result = -minResult;
  }
  else if( function == "NMAG" )
  {
    result = -magResult;
  }
  else if( function == "NRMS" )
  {
    result = -rmsResult;
  }

  return result;
}

} // namespace IO
} // namespace Xyce
