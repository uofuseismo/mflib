#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mpi.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/mpi.hpp"

using namespace MFLib;

/// Broadcast the waveform template
void MFLib::MPI::Broadcast(WaveformTemplate &tplate, const int root,
                           const MPI_Comm comm)
{
    // Figure out who I am
    int myid;
    MPI_Comm_rank(root, &myid);
    // Variables to broadcast
    std::vector<double> signal;
    double *signalPtr = nullptr;
    double samplingRate = 0;
    double weight =-1;
    double onsetTime =-1;
    int nSamples = 0;
    // Have root extract data
    if (root == myid)
    {
        if (tplate.haveSamplingRate())
        {
            samplingRate = tplate.getSamplingRate();
        }
        if (tplate.haveSignal())
        {
            nSamples = tplate.getSignalLength();
            signal.resize(nSamples);
            signalPtr = signal.data();
            tplate.getSignal(nSamples, &signalPtr);
        }
        if (tplate.havePhaseOnsetTime())
        {
            onsetTime = tplate.getPhaseOnsetTime();
        }
        weight = tplate.getShiftAndStackWeight();
    }
    MPI_Bcast(&samplingRate, 1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&weight,       1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&onsetTime,    1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&nSamples,     1, MPI_INTEGER, root, comm);
    if (nSamples > 0)
    {
        if (myid != root){signal.resize(nSamples);}
        signalPtr = signal.data();
        MPI_Bcast(signalPtr, nSamples, MPI_DOUBLE_PRECISION, root, comm);
    }
    // Now set the data on slaves
    if (myid != root)
    {
        if (samplingRate > 0){tplate.setSamplingRate(samplingRate);}
        tplate.setShiftAndStackWeight(weight);
        if (nSamples > 0){tplate.setSignal(nSamples, signal.data());}
        if (onsetTime > 0){tplate.setPhaseOnsetTime(onsetTime);}
    }
}
