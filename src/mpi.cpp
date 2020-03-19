#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <limits>
#include <mpi.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/networkStationPhase.hpp"
#include "mflib/mpi.hpp"

using namespace MFLib;

namespace
{

void broadcastString(std::string &value, const int root, const MPI_Comm comm)
{
    int myid;
    int len = 0;
    MPI_Comm_rank(root, &myid);
    if (myid == root)
    {
        len = static_cast<int> (value.size());
    } 
    MPI_Bcast(&len, 1, MPI_INTEGER, root, comm);
    if (myid != root)
    {
        value.resize(len);
    }
    char *ptr = value.data();
    MPI_Bcast(ptr, len, MPI_CHAR, root, comm);
}

void broadcastIdentifier(std::pair<MFLib::NetworkStationPhase, uint64_t> *id,
                         const int root, const MPI_Comm comm)
{
    int myid;
    uint64_t idval;
    std::string network, station, phase;
    MPI_Comm_rank(root, &myid);
    if (myid == root)
    {
        network = id->first.getNetwork();
        station = id->first.getStation();
        phase   = id->first.getPhase();
        idval   = id->second;
    } 
    broadcastString(network, root, comm);
    broadcastString(station, root, comm);
    broadcastString(phase,   root, comm);
    MPI_Bcast(&idval, 1, MPI_UINT64_T, root, comm); 
    if (myid != root)
    {
        id->first.setNetwork(network);
        id->first.setStation(station);
        id->first.setPhase(phase);
        id->second = idval;
    }
}

}

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
    double travelTime =-1;
    double magnitude = std::numeric_limits<double>::lowest();
    std::pair<MFLib::NetworkStationPhase, uint64_t> identifier;
    bool haveIdentifier;
    int nSamples = 0;
    int ipolarity = 0;
    MFLib::Polarity polarity = MFLib::Polarity::UNKNOWN;
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
        if (tplate.haveTravelTime())
        {
            travelTime = tplate.getTravelTime();
        }
        if (tplate.haveMagnitude())
        {
            magnitude = tplate.getMagnitude();
        }
        haveIdentifier = tplate.haveIdentifier();
        if (haveIdentifier)
        {
            identifier = tplate.getIdentifier();
        }
        weight = tplate.getShiftAndStackWeight();
        ipolarity = static_cast<int> (tplate.getPolarity());
    }
    MPI_Bcast(&samplingRate, 1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&weight,       1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&onsetTime,    1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&travelTime,   1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&magnitude,    1, MPI_DOUBLE_PRECISION, root, comm);
    MPI_Bcast(&nSamples,     1, MPI_INTEGER, root, comm);
    MPI_Bcast(&haveIdentifier, 1, MPI_C_BOOL, root, comm);
    if (nSamples > 0)
    {
        if (myid != root){signal.resize(nSamples);}
        signalPtr = signal.data();
        MPI_Bcast(signalPtr, nSamples, MPI_DOUBLE_PRECISION, root, comm);
    }
    if (haveIdentifier)
    {
        broadcastIdentifier(&identifier, root, comm);
    }
    MPI_Bcast(&ipolarity,   1, MPI_INTEGER, root, comm);
    polarity = static_cast<MFLib::Polarity> (ipolarity);
    // Now set the data on slaves
    if (myid != root)
    {
        if (samplingRate > 0){tplate.setSamplingRate(samplingRate);}
        tplate.setShiftAndStackWeight(weight);
        if (nSamples > 0){tplate.setSignal(nSamples, signal.data());}
        if (onsetTime > 0){tplate.setPhaseOnsetTime(onsetTime);}
        if (travelTime >= 0){tplate.setTravelTime(travelTime);}
        if (magnitude > std::numeric_limits<double>::lowest())
        {
            tplate.setMagnitude(magnitude);
        }
        tplate.setPolarity(polarity);
        if (haveIdentifier)
        {
            tplate.setIdentifier(identifier);
        }
    }
}
