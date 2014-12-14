/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>

#include <N_UTL_StatMetricTraits.h>
#include <N_UTL_CPUTime.h>
#include <N_UTL_WallTime.h>
#include <N_UTL_MallocUsed.h>
#include <N_UTL_FormatTime.h>
#include <N_UTL_FormatMemorySize.h>

namespace Xyce {
namespace Stats {

namespace {

int s_timeFormat = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS;

} // namespace <empty>


int
getStatTimeFormat() 
{
  return s_timeFormat;
}

void
setStatTimeFormat(
  int           time_format)
{
  s_timeFormat = time_format;
}


MetricTraits<LapCount>::Type
MetricTraits<LapCount>::value_now()
{
  return 1;
}

MetricTraits<CPUTime>::Type
MetricTraits<CPUTime>::value_now()
{
  return Xyce::cpu_time();
}

MetricTraits<WallTime>::Type
MetricTraits<WallTime>::value_now()
{
  return Xyce::wall_time();
}

MetricTraits<MPICount>::Type
MetricTraits<MPICount>::value_now()
{
  return 0;
}

MetricTraits<MPIByteCount>::Type
MetricTraits<MPIByteCount>::value_now()
{
  return 0;
}

MetricTraits<HeapAlloc>::Type
MetricTraits<HeapAlloc>::value_now()
{
  return ::malloc_used();
}

std::string
MetricTraits<LapCount>::table_header() {
  return "Count";
}

std::string
MetricTraits<CPUTime>::table_header() {
  return "CPU Time";
}

std::string
MetricTraits<WallTime>::table_header() {
  return "Wall Time";
}

std::string
MetricTraits<MPICount>::table_header() {
  return "MPI Count";
}

std::string
MetricTraits<MPIByteCount>::table_header() {
  return "MPI Byte Count";
}

std::string
MetricTraits<HeapAlloc>::table_header() {
  return "Heap Allocated";
}


std::string
MetricTraits<CPUTime>::format(
  MetricTraits<CPUTime>::Type           time)
{
  return formatTime(time, getStatTimeFormat());
}


std::string
MetricTraits<WallTime>::format(
  MetricTraits<WallTime>::Type          time)
{
  return formatTime(time, getStatTimeFormat());
}


std::string
MetricTraits<MPICount>::format(
  MetricTraits<MPICount>::Type          count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}


std::string
MetricTraits<MPIByteCount>::format(
  MetricTraits<MPIByteCount>::Type      count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}

std::string
MetricTraits<HeapAlloc>::format(
  MetricTraits<HeapAlloc>::Type         count)
{
  return formatMemorySize(count);
}

} // namespace Stats
} // namespace Xyce
