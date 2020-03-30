#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <vector>
#include <string>
#include <limits>
#include <mpi.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/networkStationPhase.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/mpi.hpp"

using namespace MFLib;

namespace
{

struct MPIDetection
{
    double mAmplitude1;// = 0;
    double mAmplitude2;// = 0;
    double mCorrelationCoefficient;// = 0;
    double mDetectionTime;// = 0;

    double mInterpolatedDetectionTime;// = 0;
    double mPhaseOnsetTime;// = 0;
    double mInterpolatedPhaseOnsetTime;// = 0;
    double mTravelTime;// = 0;

    char mNetwork[64];
    char mStation[64];
    char mPhase[64];
    uint64_t mTemplateID;// = 0;

    int mSignalLength;// = 0;

    bool mHaveAmplitude1;
    bool mHaveAmplitude2;

    bool mHaveCorrelationCoefficient;// = false;
    bool mHaveDetectedSignal;// = false;
    bool mHaveDetectionTime;// = false;
    bool mHaveInterpolatedDetectionTime;// = false;

    bool mHavePhaseOnsetTime;// = false;
    bool mHaveInterpolatedPhaseOnsetTime;// = false;
    bool mHaveTravelTime;// = false;
    bool mHaveTemplateID;// = false;
};

MPI_Datatype createMPIDetectionType()
{
    const int nItems = 23;
    MPI_Datatype types[23] = 
    {
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
        MPI_CHAR, MPI_CHAR, MPI_CHAR, 
        MPI_UINT64_T,
        MPI_INT,
        MPI_C_BOOL, MPI_C_BOOL,
        MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL,
        MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL
    };
    int blockLengths[23] = 
    {
        1, 1, 1, 1,
        1, 1, 1, 1, 
        64, 64, 64,
        1,
        1,
        1, 1, 
        1, 1, 1, 1,
        1, 1, 1, 1
    };        
    MPI_Aint offset[23];
    offset[0]  = offsetof(MPIDetection, mAmplitude1);
    offset[1]  = offsetof(MPIDetection, mAmplitude2);
    offset[2]  = offsetof(MPIDetection, mCorrelationCoefficient);
    offset[3]  = offsetof(MPIDetection, mDetectionTime);

    offset[4]  = offsetof(MPIDetection, mInterpolatedDetectionTime);
    offset[5]  = offsetof(MPIDetection, mPhaseOnsetTime);
    offset[6]  = offsetof(MPIDetection, mInterpolatedPhaseOnsetTime);
    offset[7]  = offsetof(MPIDetection, mTravelTime);

    offset[8]  = offsetof(MPIDetection, mNetwork);
    offset[9]  = offsetof(MPIDetection, mStation);
    offset[10] = offsetof(MPIDetection, mPhase);
    offset[11] = offsetof(MPIDetection, mTemplateID);
    offset[12] = offsetof(MPIDetection, mSignalLength);

    offset[13]  = offsetof(MPIDetection, mAmplitude1);
    offset[14]  = offsetof(MPIDetection, mAmplitude2);

    offset[15]  = offsetof(MPIDetection, mHaveCorrelationCoefficient);
    offset[16]  = offsetof(MPIDetection, mHaveDetectedSignal);
    offset[17]  = offsetof(MPIDetection, mHaveDetectionTime);
    offset[18]  = offsetof(MPIDetection, mHaveInterpolatedDetectionTime);

    offset[19]  = offsetof(MPIDetection, mHavePhaseOnsetTime);
    offset[20]  = offsetof(MPIDetection, mHaveInterpolatedPhaseOnsetTime);
    offset[21]  = offsetof(MPIDetection, mHaveTravelTime);
    offset[22]  = offsetof(MPIDetection, mHaveTemplateID);

    // Chars and bools can create potential holes in alignment
    MPI_Datatype mpiDetection, tempType; 
    MPI_Aint lb, extent;
    MPI_Type_create_struct(nItems, blockLengths, offset, types, &tempType);
    MPI_Type_get_extent(tempType, &lb, &extent);
    MPI_Type_create_resized(tempType, lb, extent, &mpiDetection);
    MPI_Type_commit(&mpiDetection);
    MPI_Type_free(&tempType);
    return mpiDetection;
}

template<class T>
MPIDetection convertDetection(const MFLib::SingleChannel::Detection<T> &det)
{
    MPIDetection detOut;
    std::memset(&detOut, 0, sizeof(MPIDetection));
    std::memset(detOut.mNetwork, 0, 64*sizeof(char));
    std::memset(detOut.mStation, 0, 64*sizeof(char));
    std::memset(detOut.mPhase,   0, 64*sizeof(char));
    
    detOut.mHaveAmplitude1 = det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
    detOut.mHaveAmplitude2 = det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
    detOut.mHaveCorrelationCoefficient = det.haveCorrelationCoefficient();
    detOut.mHaveDetectedSignal = det.haveDetectedSignal();
    detOut.mHaveDetectionTime = det.haveDetectionTime();
    detOut.mHaveInterpolatedDetectionTime = det.haveInterpolatedDetectionTime();
    detOut.mHavePhaseOnsetTime = det.havePhaseOnsetTime();
    detOut.mHaveInterpolatedPhaseOnsetTime
        = det.haveInterpolatedPhaseOnsetTime();
    detOut.mHaveTravelTime = det.haveTravelTime();
    detOut.mHaveTemplateID = det.haveTemplateIdentifier();
    if (detOut.mHaveDetectedSignal)
    {
        detOut.mSignalLength = det.getDetectedSignalLength();
    }
    if (detOut.mHaveAmplitude1)
    {
        detOut.mAmplitude1 = det.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
    }
    if (detOut.mHaveAmplitude2)
    {
        detOut.mAmplitude2 = det.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
    }
    if (detOut.mHaveCorrelationCoefficient)
    {
        detOut.mCorrelationCoefficient = det.getCorrelationCoefficient();
    }
    if (detOut.mHaveDetectionTime)
    {
        detOut.mDetectionTime = det.getDetectionTime();
    }
    if (detOut.mHaveInterpolatedDetectionTime)
    {
        detOut.mInterpolatedDetectionTime = det.getInterpolatedDetectionTime();
    } 
    if (detOut.mHavePhaseOnsetTime)
    {
        detOut.mPhaseOnsetTime = det.getPhaseOnsetTime();
    }
    if (detOut.mHaveInterpolatedPhaseOnsetTime)
    {
        detOut.mInterpolatedPhaseOnsetTime
            = det.getInterpolatedPhaseOnsetTime();
    }
    if (detOut.mHaveTravelTime)
    {
        detOut.mTravelTime = det.getTravelTime();
    }
    if (detOut.mHaveTemplateID)
    {
        auto id = det.getTemplateIdentifier();
        auto network = id.first.getNetwork();
        auto station = id.first.getStation();
        auto phase   = id.first.getPhase();
        detOut.mTemplateID = id.second;
        auto lenc = std::min(static_cast<size_t> (64), network.size()); 
        strncpy(detOut.mNetwork, network.c_str(), lenc);
        lenc = std::min(static_cast<size_t> (64), station.size());
        strncpy(detOut.mStation, station.c_str(), lenc);
        lenc = std::min(static_cast<size_t> (64), phase.size());
        strncpy(detOut.mPhase,   phase.c_str(), lenc);
    }
    return detOut;
}

template<typename T>
void convertDetection(const MPIDetection &det,
                      MFLib::SingleChannel::Detection<T> *detOut)
{
    if (det.mHaveAmplitude1)
    {
        detOut->setAmplitudeScalingFactor(det.mAmplitude1,
            MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
    }
    if (det.mHaveAmplitude2)
    {
        detOut->setAmplitudeScalingFactor(det.mAmplitude2,
            MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014); 
    }
    if (det.mHaveCorrelationCoefficient)
    {
        detOut->setCorrelationCoefficient(det.mCorrelationCoefficient);
    }
    if (det.mHaveDetectionTime)
    {
        detOut->setDetectionTime(det.mDetectionTime);
    }
    if (det.mHaveInterpolatedDetectionTime)
    {
        detOut->setInterpolatedDetectionTime(det.mInterpolatedDetectionTime);
    }
    if (det.mHavePhaseOnsetTime)
    {
        detOut->setPhaseOnsetTime(det.mPhaseOnsetTime);
    }
    if (det.mHaveInterpolatedPhaseOnsetTime)
    {
        detOut->setInterpolatedPhaseOnsetTime(det.mInterpolatedPhaseOnsetTime);
    } 
    if (det.mHaveTravelTime)
    {
        detOut->setTravelTime(det.mTravelTime);
    }
    if (det.mHaveTemplateID)
    {
        MFLib::NetworkStationPhase nsp;
        nsp.setNetwork(std::string(det.mNetwork));
        nsp.setStation(std::string(det.mStation));
        nsp.setPhase(std::string(det.mPhase));
        detOut->setTemplateIdentifier(std::pair(nsp, det.mTemplateID));
    }
}


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

template<>
void MFLib::MPI::Send(const MFLib::SingleChannel::Detection<double> &detection,
                      const int dest, const int tag, const MPI_Comm comm)
{
    auto mpiDetection = convertDetection(detection);
    auto MPI_DETECTION_TYPE = createMPIDetectionType(); 
    MPI_Send(&mpiDetection, 1, MPI_DETECTION_TYPE,  dest, tag, comm);
    if (mpiDetection.mHaveDetectedSignal)
    {
        auto signalPtr = detection.getDetectedSignalPointer();
        MPI_Send(signalPtr, mpiDetection.mSignalLength, MPI_DOUBLE,
                 dest, tag, comm);
    }
}

template<>
void MFLib::MPI::Send(const MFLib::SingleChannel::Detection<float> &detection,
                      const int dest, const int tag, const MPI_Comm comm)
{
    auto mpiDetection = convertDetection(detection);
    auto MPI_DETECTION_TYPE = createMPIDetectionType();
    MPI_Send(&mpiDetection, 1, MPI_DETECTION_TYPE,  dest, tag, comm);
    if (mpiDetection.mHaveDetectedSignal)
    {
        auto signalPtr = detection.getDetectedSignalPointer();
        MPI_Send(signalPtr, mpiDetection.mSignalLength, MPI_FLOAT,
                 dest, tag, comm);
    }
}

template<>
void MFLib::MPI::Recv(MFLib::SingleChannel::Detection<double> *detection,
                      const int source, const int tag,
                      const MPI_Comm comm, MPI_Status *status)
{
    MPIDetection mpiDetection;// = convertDetection(detection);
    auto MPI_DETECTION_TYPE = createMPIDetectionType();
    MPI_Recv(&mpiDetection, 1, MPI_DETECTION_TYPE, source, tag, comm, status);
    convertDetection(mpiDetection, detection);
    if (mpiDetection.mHaveDetectedSignal)
    {
        std::vector<double> signal(mpiDetection.mSignalLength);
        MPI_Recv(signal.data(), mpiDetection.mSignalLength, MPI_DOUBLE,
                 source, tag, comm, status);
        detection->setDetectedSignal(signal.size(), signal.data());
    }
}

template<>
void MFLib::MPI::Recv(MFLib::SingleChannel::Detection<float> *detection,
                      const int source, const int tag,
                      const MPI_Comm comm, MPI_Status *status)
{
    MPIDetection mpiDetection;
    auto MPI_DETECTION_TYPE = createMPIDetectionType();
    MPI_Recv(&mpiDetection, 1, MPI_DETECTION_TYPE, source, tag, comm, status);
    convertDetection(mpiDetection, detection);
    if (mpiDetection.mHaveDetectedSignal)
    {
        std::vector<float> signal(mpiDetection.mSignalLength);
        MPI_Recv(signal.data(), mpiDetection.mSignalLength, MPI_FLOAT,
                 source, tag, comm, status);
        detection->setDetectedSignal(signal.size(), signal.data());
    }
}
