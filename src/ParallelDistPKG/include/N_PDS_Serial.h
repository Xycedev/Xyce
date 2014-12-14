#ifndef XYCE_UTIL_Serial_h
#define XYCE_UTIL_Serial_h

#include <stdexcept>
#include <string>
#include <vector>

#include <N_PDS_fwd.h>
#include <N_PDS_ParallelMachine.h>

#ifndef Xyce_PARALLEL_MPI

static const int MPI_MAX = 0;
static const int MPI_MIN = 0;
static const int MPI_SUM = 0;
static const int MPI_PROD = 0;
static const int MPI_LAND = 0;
static const int MPI_BAND = 0;
static const int MPI_LOR = 0;
static const int MPI_BOR = 0;
static const int MPI_LXOR = 0;
static const int MPI_BXOR = 0;

namespace Xyce {
namespace Parallel {

inline int op_identifier_compare_op() {
  return 0;
}

template<class T>
inline void
AllReduce(Machine comm, int op, T *src_dest, size_t size)
{}

template<class T>
inline void
AllReduce(Machine comm, int op, const T *source, T *dest, size_t size)
{
  std::copy(source, source + size, dest);
}

template<class T>
inline void
AllReduce(Machine comm, int op, std::vector<T> &src_dest)
{}

template<class T>
inline void
AllGather(Machine mpi_comm, const std::vector<T> &source, std::vector<T> &dest)
{
  if (source.size() != dest.size())
    throw std::runtime_error("Xyce::Serial::AllGather(MPI_Comm mpi_comm, std::vector<T> &source, std::vector<T> &dest) vector lengths not equal");

  dest = source;
}

template<class T>
inline void
AllGather(Machine mpi_comm, const T &source, std::vector<T> &dest)
{
  if (dest.size() != 1)
    throw std::runtime_error("Xyce::Serial::AllGather(MPI_Comm mpi_comm, const T &source, std::vector<T> &dest) vector lengths not equal");

  dest[0] = source;
}

template<class T>
inline void
Broadcast(Machine comm, T *src_dest, size_t len, int root)
{}

inline void
AllWriteString(
  Machine               comm,
  std::ostream &        os,
  const std::string &   message)
{
  os << message;
}

inline void Barrier(Machine comm)
{}

inline void
GatherV(
  Machine                       comm,
  unsigned                      root,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  dest.resize(1);
  dest[0] = src;
}

template<class T>
inline void
GatherV(
  Machine                       mpi_comm,
  unsigned                      root,
  const std::vector<T> &        src,
  std::vector<T> &              dest)
{
  dest = src;
}

inline void
AllGatherV(
  Machine                       comm,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  dest.resize(1);
  dest[0] = src;
}

template<class T>
inline void
AllGatherV(
  Machine                       mpi_comm,
  const std::vector<T> &        src,
  std::vector<T> &              dest)
{
  dest = src;
}

} // namespace Parallel
} // namespace Xyce

#endif // Xyce_PARALLEL_MPI

#endif // Xyce_UTIL_Serial_h
