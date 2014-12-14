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
// Filename       : $RCSfile: N_IO_DistributionTool.C,v $
//
// Purpose        : Defines the DistributionTool class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.97.2.2 $
//
// Revision Date  : $Date: 2014/08/29 21:26:34 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_IO_fwd.h>
#include <N_IO_DistributionTool.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_PDS_Comm.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistributionTool::DistributionTool
// Purpose       : ctor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistributionTool::DistributionTool( CircuitBlock * cktBlk,
                                    CmdParse & cp,
                                    N_PDS_Comm * pdsCommPtr = NULL )
: cktBlk_( cktBlk ),
  numProcs_(1),
  currProc_(0),
  commandLine_(cp),
  pdsCommPtr_( pdsCommPtr ),
  procDeviceCount_(0),
  devices_(0),
  deviceLinesSent_(0),
  charBufferSize_(250000),
  charBufferPos_(0),
  charBuffer_(0),
  circuitContextReady_(false),
  circuitOptionsReady_(false)
{ 
  if ( pdsCommPtr_ )
  {
    numProcs_ = pdsCommPtr_->numProc();

    // procID == 0 may not be true if we're running in a hierarchical parallel context
    // so check if numProc() == 1 too.
    if (pdsCommPtr_->procID() == 0)
    {
      ( numProcs_ == 1) ? currProc_ = 0 : currProc_ = 1;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::~DistributionTool
// Purpose       : dtor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistributionTool::~DistributionTool()
{
  if( charBuffer_ != 0 )
  {
    delete [] charBuffer_;
  }
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool DistributionTool::registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::deviceCount
// Purpose       : set total device count and determine device/proc ratio
// Special Notes : a noop in serial
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::deviceCount( int devices )
{
  // determine how many devices to send to each far proc
  devices_ = devices;
  procDeviceCount_ = devices / numProcs_;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::setCircuitContext
// Purpose       : send circuit context to all procs
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::setCircuitContext( CircuitContext * const
 circuitContexts )
{
  // store until explicit send
  circuitContexts_ = circuitContexts;

#ifdef Xyce_PARALLEL_MPI

  // resize buffer if necessary
  int tmpSize = circuitContexts_->packedByteCount() + sizeof( int );

  if ( tmpSize > charBufferSize_ )
  {
    charBufferSize_ = tmpSize;
  }

  // set flag to tx global data
  circuitContextReady_ = true;

#endif

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::packCircuitContext
// Purpose       : send circuit context to all procs
// Special Notes : 
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int DistributionTool::packCircuitContext()
{
  int bsize = 0;

#ifdef Xyce_PARALLEL_MPI

  if ( circuitContextReady_ )
  {
    int pos = 0;

    // pack circuit contexts
    circuitContexts_->pack( charBuffer_, charBufferSize_, pos, pdsCommPtr_ );

    bsize=pos;
  } 

#endif

  return bsize;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::unpackCircuitContext
// Purpose       : receive circuit context from proc 0
// Special Notes : 
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::unpackCircuitContext(int bsize)
{
#ifdef Xyce_PARALLEL_MPI

  if ( !circuitContextReady_ )
  {
    int pos = 0;

    // unpack circuit contexts; set top level parent to NULL
    circuitContexts_ = cktBlk_->getCircuitContextPtr();

    circuitContexts_->setParentContextPtr( NULL );
    circuitContexts_->unpack( charBuffer_, bsize, pos, pdsCommPtr_ );

    // parse circuit context on far end
    cktBlk_->receiveCircuitContext( *circuitContexts_ );
  }

#endif

  return true;

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::setCircuitOptions
// Purpose       : stage options for tx
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::setCircuitOptions( const std::list< Util::OptionBlock > &
  options )
{
  // store until explicit send
  options_ = options;

#ifdef Xyce_PARALLEL_MPI

  // resize buffer if necessary
  int tmpSize = sizeof( int );
  std::list< Util::OptionBlock >::const_iterator it_obL = options_.begin();
  std::list< Util::OptionBlock >::const_iterator it_oeL = options_.end();

  for( ; it_obL != it_oeL; ++it_obL )
  {
    tmpSize += it_obL->packedByteCount();
  }

  if ( tmpSize > charBufferSize_ )
  {
    charBufferSize_ = tmpSize;
  }

  // set flag to tx global data
  circuitOptionsReady_ = true;

#endif

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::packCircuitOptions
// Purpose       : send option blocks to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int DistributionTool::packCircuitOptions()
{
  int bsize = 0;

#ifdef Xyce_PARALLEL_MPI

  if ( circuitOptionsReady_ )
  {
    int pos = 0;

    std::list< Util::OptionBlock >::iterator it_obL = options_.begin();
    std::list< Util::OptionBlock >::iterator it_oeL = options_.end();

    // pack options
    int count = options_.size();
    pdsCommPtr_->pack( &count, 1, charBuffer_, charBufferSize_, pos );
    for( ; it_obL != it_oeL; ++it_obL )
    {
      it_obL->pack( charBuffer_, charBufferSize_, pos, pdsCommPtr_ );
    }

    bsize = pos;
  }

#endif

  return bsize;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::unpackCircuitOptions
// Purpose       : unpack option blocks from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::unpackCircuitOptions(int bsize)
{
#ifdef Xyce_PARALLEL_MPI

  if ( !circuitOptionsReady_ )
  {
    int pos = 0;
    int size;

    // unpack options
    pdsCommPtr_->unpack( charBuffer_, bsize, pos, &size, 1 );
    for( int i = 0; i < size; ++i )
    {
      Util::OptionBlock anOptionBlock;
      anOptionBlock.unpack( charBuffer_, bsize, pos, pdsCommPtr_ );
      options_.push_back( anOptionBlock );
    }
  }

#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::registerCircuitOptions
// Purpose       : register options
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::registerCircuitOptions()
{

  std::string netListFile("");
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  // validate the ref count pointer we're going to use
  if( pkgOptMgrPtr_ == 0 )
  {
    Report::DevelFatal().in("DistributionTool::registerOptions()") << "Called with null PkgOptionsMgr";
  }

  std::list<N_UTL_OptionBlock>::iterator iterOB = options_.begin();
  std::list<N_UTL_OptionBlock>::iterator endOB = options_.end();

  //register options with pkg option mgr (new options control technique)
  for( ; iterOB != endOB; ++iterOB )
  {
    if( (iterOB->getName() == "GLOBAL") && (iterOB->getStatus() != Xyce::Util::PROCESSED_STATE) )
    {
      pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
      iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
    }
  }

  for( iterOB = options_.begin(); iterOB != endOB; ++iterOB )
  {
    if( iterOB->getName() != "GLOBAL"
        && iterOB->getName() != "PARALLEL"
        && iterOB->getName() != "PRINT" 
        && iterOB->getName() != "SENS" )
    {
      if(iterOB->getStatus() != Xyce::Util::PROCESSED_STATE)
      {
        pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
        iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
      }
    }
  }

  // Yet another hack so that SENS comes before print.
  for( iterOB = options_.begin(); iterOB != endOB; ++iterOB )
  {
    if( iterOB->getName() == "SENS"
        && (iterOB->getStatus() != Xyce::Util::PROCESSED_STATE))
    {
      pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
      iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
    }
  }

  for( iterOB = options_.begin(); iterOB != endOB; ++iterOB )
  {
    if( iterOB->getName() == "PRINT"
        && (iterOB->getStatus() != Xyce::Util::PROCESSED_STATE))
    {
      pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
      iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitDeviceLine
// Purpose       : Send a circuit device line to current proc
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitDeviceLine(std::vector< SpiceSeparatedFieldTool::StringToken > & deviceLine )
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    int size = 0;

    // count line type
    size += sizeof(char);

    // count line size
    size += sizeof(int);

    // count device line
    int deviceLineSize = deviceLine.size();
    for (int i = 0; i < deviceLineSize; ++i)
    {
      size += deviceLine[i].packedByteCount();
    }

    // flush buffer as needed
    send(size);

    // pack line type; "d"evice line
    char lineType = 'd';
    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack device line size
    pdsCommPtr_->pack( &deviceLineSize, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack device line
    for (int i = 0; i < deviceLineSize; ++i)
    {
      deviceLine[ i ].pack( charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
    }

    // increment device line counter
    deviceLinesSent_++;

    // manage proc change if n/p reached
    if (deviceLinesSent_ >= procDeviceCount_)
    {
      int minus1 = -1;

      // flush buffer
      send();

      // allow this proc to exit receive loop and begin processing
      pdsCommPtr_->send( &minus1, 1, currProc_ );

      // reset device line counter
      deviceLinesSent_ = 0;

      // move to next processor
      currProc_++;

      // switch to proc 0 if nearly complete
      if (currProc_ == numProcs_)
      {
        currProc_ = 0;
      }

      // otherwise, prepare currProc for next set of device lines
      else
      {
        int length = fileName_.size();
        lineType = 'f';

        // flush buffer as needed
        send(sizeof(char) + sizeof(int) + length);

        // pack the filename
        pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
        pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
        pdsCommPtr_->pack( fileName_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

        // pack the subcircuit names
        int subcircuitNamesSize = subcircuitNames_.size();
        for (int i = 0; i < subcircuitNamesSize; ++i)
        {
          if( currProc_ != 0 )
          {
            packSubcircuitData(subcircuitNames_[i], subcircuitNodes_[i],
                               subcircuitPrefixes_[i], subcircuitParams_[i]);
          }
        }
      }
    }

    return true;
  }

  else

#endif

    // let proc 0 parse line
    return false;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::endDeviceLines()
// Purpose       : Make sure other processors have been told that device processing
//                 is ended.  This can be a problem with small device count circuits
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/06/06
//-----------------------------------------------------------------------------
void DistributionTool::endDeviceLines()
{

#ifdef Xyce_PARALLEL_MPI

  if (currProc_ > 0)
  {
    int minus1 = -1;

    // flush buffer
    send();

    // end parsing on current node
    pdsCommPtr_->send( &minus1, 1, currProc_ );
    ++currProc_;

    // end parsing on remaining nodes
    for ( ; currProc_ < numProcs_ ; ++currProc_)
    {
      pdsCommPtr_->send( &minus1, 1, currProc_ );
    }

    // complete any remaining parser action on node 0
    currProc_ = 0;
  }

#endif

}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitStart
// Purpose       : change current subcircuit context after a new subcircuit is
//               : started
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitStart( std::string const & subcircuitName,
                                          std::list<std::string> const & nodes,
					  std::string const & prefix,
                                          std::vector<N_DEV_Param> const & params )
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    // save new context data on stacks
    subcircuitNames_.push_back( subcircuitName );
    subcircuitPrefixes_.push_back( prefix );
    subcircuitNodes_.push_back( nodes );
    subcircuitParams_.push_back( params );

    // pack data into buffer
    packSubcircuitData( subcircuitName, nodes, prefix, params );
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::setFileName
// Purpose       : Change name of netlist file
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DistributionTool::setFileName(std::string const & fileNameIn)
{

#ifdef Xyce_PARALLEL_MPI

  fileName_ = fileNameIn;

  char lineType = 'f';
  int length = fileName_.size();

  // flush buffer as needed
  send(sizeof(char) + sizeof(int) + length);

  pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( fileName_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

#endif

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::packSubcircuitData
// Purpose       : pack subcircuit data into buffer; helper for circuitStart()
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::packSubcircuitData( std::string const & subcircuitName,
                                                std::list<std::string> const & nodes,
                                                std::string const & prefix,
                                                std::vector<N_DEV_Param> const & params )
{
#ifdef Xyce_PARALLEL_MPI

  int size = 0;

  // count line type, name size, and name chars
  size += sizeof(char) + sizeof(int) + subcircuitName.length();

  // count # of nodes,
  size += sizeof(int);

  // count node sizes and node chars
  std::list<std::string>::const_iterator it_sbL = nodes.begin();
  std::list<std::string>::const_iterator it_seL = nodes.end();
  for( ; it_sbL != it_seL; ++it_sbL  )
  {
    size += sizeof(int) + it_sbL->length();
  }

  // count prefix size and prefix chars
  size += sizeof(int) + prefix.length();

  // count # of params
  size += sizeof(int);

  // count params
  std::vector<N_DEV_Param>::const_iterator it_pbL = params.begin();
  std::vector<N_DEV_Param>::const_iterator it_peL = params.end();
  for( ; it_pbL != it_peL; ++it_pbL  )
  {
    size += it_pbL->packedByteCount();
  }

  // flush buffer as needed
  send(size);


  // pack line type; subcircuit "s"tart
  char lineType = 's';
  pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack name size and name chars
  int length = subcircuitName.length();
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( subcircuitName.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack # of nodes
  size = nodes.size();
  pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack nodes name sizes and chars
  it_sbL = nodes.begin();
  for( ; it_sbL != it_seL; ++it_sbL  )
  {
    length = it_sbL->length();
    pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( it_sbL->c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );
  }

  // pack prefix size and prefix chars
  length = prefix.length();
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( prefix.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack # of params
  size = params.size();
  pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack params
  it_pbL = params.begin();
  for( ; it_pbL != it_peL; ++it_pbL  )
  {
    it_pbL->pack( charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
  }

#endif

  return true;

}



//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitEnd
// Purpose       : change current subcircuit context to previous context
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitEnd()
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    // remove latest context data on stacks
    subcircuitNames_.pop_back();
    subcircuitPrefixes_.pop_back();
    subcircuitNodes_.pop_back();
    subcircuitParams_.pop_back();

    // flush buffer as needed
    send(sizeof(char));

    // pack line type; subcircuit "e"nd
    char lineType = 'e';
    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  }

#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::receiveCircuitData
// Purpose       : receive all data from proc 0 and store it in device buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::receiveCircuitData()
{
#ifdef Xyce_PARALLEL_MPI

  int bsize;
  char *currBuffer;

  // receive the device lines first
  while (true)
  {
    pdsCommPtr_->recv( &bsize, 1, 0 );

    // check for halt request
    if (bsize < 0)
    {
      break;
    }

    currBuffer = new char[bsize];
    bufs_.push_back(currBuffer);
    bufSize_.push_back(bsize);

    pdsCommPtr_->recv(currBuffer, bsize, 0);
  }

  return processDeviceBuffer();
#else
  return true;
#endif
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::send
// Purpose       : send buffer from proc 0 and adjust buffer limits
// Special Notes : size == -1, and no args send(), will force a transmission
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::send(int size)
{
#ifdef Xyce_PARALLEL_MPI

  // transmit buffer if next set of object will fill it completely, or forced
  if ((charBufferPos_ + size >= charBufferSize_) || (size == -1))
  {

#ifdef Xyce_DEBUG_DISTRIBUTED_PARSER

    cerr << "node " << pdsCommPtr_->procID() << " packed "
         << charBufferPos_ << " bytes into " << charBufferSize_ << " byte buffer"
         << " [data size:  " << size << "]" << endl
         << "node " << currProc_ << " was sent "
         << deviceLinesSent_ << " of " << procDeviceCount_ << " devices "
         << (float)deviceLinesSent_ / procDeviceCount_ * 100.0
         << ( size == -1 ? "% complete  [flushed]" : "% complete"  ) << endl;

#endif

    // tx buffer size
    pdsCommPtr_->send( &charBufferPos_, 1, currProc_ );

    // tx buffer contents
    pdsCommPtr_->send( charBuffer_, charBufferPos_, currProc_ );

    // reset counters
    charBufferPos_ = 0;

    // adjust buffer size if objects larger than existing buffer
    if (size > charBufferSize_)
    {

      // free memory
      // delete [] charBuffer_;

      // set new buffer size
      charBufferSize_ = size;

      // reserve memory; including offset
      charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

#ifdef Xyce_DEBUG_DISTRIBUTED_PARSER

      cerr << "node " << pdsCommPtr_->procID()
           << " resized buffer to " << charBufferSize_
           << " with offset " << sizeof(char) + sizeof(int) << endl;

#endif

    }
  }

#endif

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::processDeviceBuffer
// Purpose       : process all data in the buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::processDeviceBuffer()
{
#ifdef Xyce_PARALLEL_MPI

  int size, bsize, length, pos, i;
  char lineType;

  char *currBuffer;
  unsigned int ib;

  // process received device lines, contexts info, and exit calls
  for (ib=0 ; ib<bufs_.size() ; ++ib)
  {
    currBuffer = bufs_[ib];
    bsize = bufSize_[ib];

    // get ready to process buffer
    pos = 0;

    // break apart buffered lines and parse
    while( pos < bsize )
    {
      // get the linetype marker that indicates incoming line
      pdsCommPtr_->unpack( currBuffer, bsize, pos, &lineType, 1 );

      // process line according to type
      switch( lineType )
      {
        case 'd': // process a device line
        {
          // unpack the device line
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
          std::vector< SpiceSeparatedFieldTool::StringToken > deviceLine( size );

          for( i = 0; i < size; ++i )
          {
            deviceLine[i].unpack( currBuffer, bsize, pos, pdsCommPtr_ );
          }

          // hand to circuit block for processing
          cktBlk_->handleDeviceLine( deviceLine );

          break;
        }

        case 's': // subcircuit start found
        {
          // get the subcircuit name
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
          std::string subcircuitName(std::string( ( currBuffer + pos ), length ));
          pos += length;

          // get the nodes for this subcircuit call
          std::list<std::string> nodes;
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
            nodes.push_back( std::string( ( currBuffer + pos ), length ) );
            pos += length;
          }

          // get the prefix
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
          std::string prefix(std::string( ( currBuffer + pos ), length ));
          pos += length;

          // get params
          std::vector<N_DEV_Param> params;
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            Device::Param param;
            param.unpack( currBuffer, bsize, pos, pdsCommPtr_ );
            params.push_back( param );
          }

          // send to cktblk
          circuitContexts_->setContext( subcircuitName, prefix, nodes );
          circuitContexts_->resolve( params );

          break;
        }

        case 'e': // subcircuit end found
        {
          // adjust the circuit context pointer
          circuitContexts_->restorePreviousContext();

          break;
        }

        case 'f': // change netlist file name
        {
          pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
          fileName_ = std::string( ( currBuffer + pos ), length );
          pos += length;
          cktBlk_->setFileName(fileName_);

          break;
        }

        default:  // something went wrong
        {
          Report::DevelFatal().in("DistributionTool::processDeviceBuffer")
            << "Node " << pdsCommPtr_->procID() << " received invalid message type \"" << lineType;
        }
      }
    }
    delete [] bufs_[ib];
  }

  bufs_.clear();
  bufSize_.clear();
#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::broadcastGlobalData
// Purpose       : send bufsize, options, metatdata, and context to procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::broadcastGlobalData()
{
#ifdef Xyce_PARALLEL_MPI

  // tx buffer size
  pdsCommPtr_->bcast( &charBufferSize_, 1, 0 );

  // check for halt request
  if (charBufferSize_ < 0)
  {
    return false;
  }

  // reserve memory; including offset
  charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

  // make sure data is ready to tx
  if (numProcs_ > 1)
  {
    // tx global data
    int bsize = packCircuitOptions();

    // broadcast options
    pdsCommPtr_->bcast( &bsize, 1, 0 );
    pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

    // receive global data
    unpackCircuitOptions( bsize );
 
    // tx global data
    bsize = packCircuitContext();
    
    // broadcast context
    pdsCommPtr_->bcast( &bsize, 1, 0 );
    pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

    unpackCircuitContext( bsize );
  }
 
  // transform booleans into integers and transmit them
  bool netlistcopy = commandLine_.getHangingResistor().getNetlistCopy();
  bool oneterm = commandLine_.getHangingResistor().getOneTerm();
  bool nodcpath = commandLine_.getHangingResistor().getNoDCPath();

  int nlcopyint = 0;
  int otint = 0;
  int nodcint = 0;

  if (netlistcopy)
    nlcopyint = 1;
  if (oneterm)
    otint = 1;
  if (nodcpath)
    nodcint = 1;

  // tx global data
  pdsCommPtr_->bcast(&nlcopyint,1,0);
  pdsCommPtr_->bcast(&otint,1,0);
  pdsCommPtr_->bcast(&nodcint,1,0);

  std::string onetermres(commandLine_.getHangingResistor().getOneTermRes());
  std::string nodcpathres(commandLine_.getHangingResistor().getNoDCPathRes());
  int otreslength = onetermres.size();
  int nodcreslength = nodcpathres.size();

  // tx global data
  pdsCommPtr_->bcast(&otreslength,1,0);
  onetermres.resize(otreslength);
  pdsCommPtr_->bcast(&onetermres[0],otreslength,0);
  pdsCommPtr_->bcast(&nodcreslength,1,0);
  nodcpathres.resize(nodcreslength);
  pdsCommPtr_->bcast(&nodcpathres[0],nodcreslength,0);

  //Set the appropriate booleans in commandLine_
  if (nlcopyint == 1)
    commandLine_.getHangingResistor().setNetlistCopy(true);
  if (otint == 1)
    commandLine_.getHangingResistor().setOneTerm(true);
  if (nodcint == 1)
    commandLine_.getHangingResistor().setNoDCPath(true);

  commandLine_.getHangingResistor().setOneTermRes(onetermres);
  commandLine_.getHangingResistor().setNoDCPathRes(nodcpathres);

#endif

  // register the global parameters
  cktBlk_->registerGlobalParams();

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());   
#endif

  // register circuit options 
  return registerCircuitOptions();
}


} // namespace IO
} // namespace Xyce
