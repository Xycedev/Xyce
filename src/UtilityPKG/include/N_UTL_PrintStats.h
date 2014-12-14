/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Xyce_N_UTL_PrintStats_h
#define Xyce_N_UTL_PrintStats_h

#include <iosfwd>

#include <N_UTL_Stats.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Stats {

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream in XML format
///
/// @param os           output stream
/// @param metrics_mask make of metrics to write
/// @param checkpoint   true if these values are checkpointed
///
/// @return output stream
///
///
std::ostream &printXML(std::ostream& os, MetricsMask metrics_mask, bool checkpoint);

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream as a table, single processor.
///
/// @param os           output stream
/// @param root_timer   root timer
/// @param metrics_mask make of metrics to write
/// @param checkpoint   true if these values are checkpointed
///
/// @return output stream
///
std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint);

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream as a table, accumulated and reported on processor rank 0.
///
/// @param os                   output stream
/// @param root_timer           root timer
/// @param metrics_mask         make of metrics to write
/// @param checkpoint           true if these values are checkpointed
/// @param parallel_machine     communicator to accumulate stats
///
/// @return output stream
///
std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint, Parallel::Machine parallel_machine);

} // namespace Stats
} // namespace Xyce

#endif // Xyce_N_UTL_PrintStats_h
