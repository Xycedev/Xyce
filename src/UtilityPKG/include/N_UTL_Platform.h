#ifndef Xyce_N_UTL_Platform_h
#define Xyce_N_UTL_Platform_h

#include <iosfwd>

namespace Xyce {

///
/// @addtogroup EnvDetail
/// @{
///

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>hostname</b> returns the hostname of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the host name obtained from
 *				the operating system.
 */
std::string hostname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>domainname</b> returns the domainname of the domain running the
 * application.
 *
 * @return			a <b>String</b> value of the domain name obtained from
 *				the operating system.
 */
std::string domainname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>username</b> returns the username of the user running the
 * application.
 *
 * @return			a <b>String</b> value of the username obtained from
 *				the operating system.
 */
std::string username();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>hardware</b> returns the hardware type of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>machine</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string hardware();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>osname</b> returns the operating system nameof the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>sysname</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string osname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>osversion</b> returns the hardware type of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>release</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string osversion();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>pid</b> returns the process id of the process running the
 * application.
 *
 * @return			a <b>int</b> value of the process id obtained from
 *				the operating system.
 */
int pid();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>pgrp</b> returns the process group id of the process running
 * the application.
 *
 * @return			a <b>int</b> value of the process group id obtained from
 *				the operating system.
 */
int pgrp();


/**
 * @brief Member function <b>get_heap_info</b> returns the amount of heap
 * memory used in bytes and the largest free block of memory in bytes.
 *
 * @param heap_size		a <b>size_t</b> returns the amount of heap
 * memory used in bytes.
 *
 * @param largest_free		a <b>size_t</b> returns the largest free block
 * of memory.
 *
 */
void get_heap_info(size_t &heap_size, size_t &largest_free);


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_heap_usage</b> returns the number of bytes used by the heap.
 *
 * @return			a <b>size_t</b> value of the number of bytes used by
 *				the heap.
 */
inline size_t get_heap_usage() {
  size_t heap_size;
  size_t largest_free;
  get_heap_info(heap_size, largest_free);

  return heap_size;
}

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_available_memory</b> returns an estimation of the amount of memory available to the process.
 *
 * @return			a <b>size_t</b> value of the number of bytes available to the process.
 */
size_t get_available_memory();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_memory_info</b> returns the total memory usage of the
 * process and the number of page faults accumulated by the process.
 *
 * @param memory_usage		a <b>size_t</b> reference to receive the number of
 *				bytes currently used by the process.
 *
 * @param faults		a <b>size_t</b> reference to treceive the number of
 *				page faults incurred by the process.
 *
 */
void get_memory_info(size_t &memory_usage, size_t &faults);

///
/// @}
///

} // namespace Xyce

#endif // Xyce_N_UTL_Platform_h
