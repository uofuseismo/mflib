#ifndef MFLIB_MPI_HPP
#define MFLIB_MPI_HPP
#include <mpi.h>
namespace MFLib
{
class WaveformTemplate;
namespace MPI
{
/*!
 * @brief Broadcasts a waveform template.
 * @param[in,out] tplate  On input, this is the waveform template defined on
 *                        the root process.  
 *                        On exit, this is the waveform template defined on
 *                        all processes on the communicator.
 * @param[in] root        The root process on the communicator.
 * @param[in] comm        The communicator.
 * @throws std::runtime_error if MPI_FAILURE was detected.
 */
void Broadcast(WaveformTemplate &tplate, int root, MPI_Comm comm);
}
}
#endif
