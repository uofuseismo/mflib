#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL 1
#endif
#include "mflib/singleChannel/associator.hpp"
#include "mflib/singleChannel/associatorParameters.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/event.hpp"
#include "mflib/networkStationPhase.hpp"
#include "private/dbscan.hpp"
#include "private/weightedStatistics.hpp"

using namespace MFLib::SingleChannel;

namespace
{

struct TemplateIndices
{
    uint64_t evid = 0;
    std::vector<size_t> indices;
};

}

template<class T>
class Associator<T>::AssociatorImpl
{
public:
    //MFLib::DBSCAN mDBSCAN;
    MFLib::SingleChannel::AssociatorParameters mParameters;
    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    std::vector<TemplateIndices> mTemplateIndices;
    int mAssociations = 0;
    bool mInitialized = false;
};

/// Constructor
template<class T>
Associator<T>::Associator() :
    pImpl(std::make_unique<AssociatorImpl> ())
{
}

/// Copy c'tor
template<class T>
Associator<T>::Associator(const Associator &associator)
{
    *this = associator;
}

/// Move c'tor
template<class T>
Associator<T>::Associator(Associator &&associator) noexcept
{
    *this = std::move(associator);
}

/// Copy assignment operator
template<class T>
Associator<T>& Associator<T>::operator=(const Associator &associator)
{
    if (&associator == this){return *this;}
    pImpl = std::make_unique<AssociatorImpl> (*associator.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
Associator<T>& Associator<T>::operator=(Associator &&associator) noexcept
{
    pImpl->mAssociations = 0;
    if (&associator == this){return *this;}
    pImpl= std::move(associator.pImpl);
    return *this;
}

/// Destructor
template<class T>
Associator<T>::~Associator() = default;

/// Releases memory on the class and resets it
template<class T>
void Associator<T>::clear() noexcept
{
    pImpl->mParameters.clear();
    pImpl->mDetections.clear();
    pImpl->mTemplateIndices.clear();
    pImpl->mAssociations = 0;
    pImpl->mInitialized = false;
}

template<class T>
void Associator<T>::addDetection(const MFLib::SingleChannel::Detection<T> &det)
{
    pImpl->mAssociations = 0;
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    if (!det.haveTravelTime())
    {
        throw std::invalid_argument("Travel time must be set\n");
    }
    if (!det.haveTemplateIdentifier())
    {
        throw std::invalid_argument("Template identifier must be set\n");
    }
    if (pImpl->mParameters.useCorrelationCoefficientWeighting())
    {
        if (!det.haveCorrelationCoefficient())
        {
            throw std::invalid_argument("XC coefficient must be set\n");
        }
    }
    pImpl->mDetections.push_back(det);
}

template<class T>
void Associator<T>::clearDetections() noexcept
{
    pImpl->mDetections.clear();
}

template<class T>
int Associator<T>::getNumberOfDetections() const noexcept
{
    return static_cast<int> (pImpl->mDetections.size());
}

/// Initialize the class
template<class T>
void Associator<T>::initialize(const AssociatorParameters &parameters)
{
    clear();
    pImpl->mParameters = parameters;
    pImpl->mInitialized = true;
}

/// Determines if the class is intitialized
template<class T>
bool Associator<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// 
template<class T>
void Associator<T>::associate()
{
    pImpl->mAssociations = 0;
    pImpl->mTemplateIndices.clear();
    if (!isInitialized()){throw std::runtime_error("Class is not intialized\n");}
    // Sort the detections based on event IDs
    printf("Sorting...\n");
    std::sort(
#ifdef USE_PSTL
              std::execution::par,
#endif
              pImpl->mDetections.begin(), pImpl->mDetections.end(),
              [] (const MFLib::SingleChannel::Detection<T> &a,
                  const MFLib::SingleChannel::Detection<T> &b)
              {
                  return a.getTemplateIdentifier().second <
                         b.getTemplateIdentifier().second;
              });
    // Get the unique event IDs
    printf("Computing unique IDs...\n");
    std::vector<uint64_t> evids;
    evids.reserve(pImpl->mDetections.size());
    for (const auto &det : pImpl->mDetections)
    {
        evids.push_back(det.getTemplateIdentifier().second);
    }
    auto uniqueEvids = evids;
#ifdef USE_PSTL
    auto it = std::unique(std::execution::par,
                          uniqueEvids.begin(), uniqueEvids.end());
#else
    auto it = std::unique(uniqueEvids.begin(), uniqueEvids.end());
#endif
    uniqueEvids.resize(std::distance(uniqueEvids.begin(), it));
    // Fire up DBSCAN
    const auto useXCWeights
        = pImpl->mParameters.useCorrelationCoefficientWeighting();
    MFLib::DBSCAN dbscan;
    dbscan.initialize(pImpl->mParameters.getOriginTimeTolerance(),
                      pImpl->mParameters.getMinimumNumberOfPicksInEvent());
    printf("Clustering...\n");
    // Cluster based on each unique event ID
    for (const auto &evid : uniqueEvids)
    {
        auto i1 = static_cast<size_t> (
                    std::lower_bound(evids.begin(), evids.end(), evid)
                  - evids.begin());
        auto i2 = static_cast<size_t> (
                    std::upper_bound(evids.begin(), evids.end(), evid)
                  - evids.begin());
        if (i2 - i1 < 1){continue;} // No data
        // Do a fine sort of the detections for this event based on 
        // the arrival times.  When we extract the extract the events
        // it will then be naturally sorted on arrival time.
        std::sort(pImpl->mDetections.begin() + i1,
                  pImpl->mDetections.begin() + i2,
                  [] (const MFLib::SingleChannel::Detection<T> &a,
                      const MFLib::SingleChannel::Detection<T> &b)
                  {
                      auto ao = a.getPhaseOnsetTime();
                      auto bo = b.getPhaseOnsetTime();
                      return ao < bo;
                  });
                // And associate this event
        // Extract the origin times - we'll cluster based on this 
        std::vector<double> originTimes(i2 - i1);
        for (auto i=i1; i<i2; ++i)
        {
            originTimes[i-i1] = pImpl->mDetections[i].getPhaseOnsetTime()
                              - pImpl->mDetections[i].getTravelTime();
        }
        // Set the data
        const int nObs = static_cast<int> (originTimes.size()); 
        const int nFeatures = 1;
        if (useXCWeights)
        {
            std::vector<double> weights(i2 - i1);
            for (auto i=i1; i<i2; ++i)
            {
                weights[i-i1]
                    = pImpl->mDetections[i].getCorrelationCoefficient();
            }
            dbscan.setWeightedData(nObs, nFeatures,
                                   originTimes.data(), weights.data());
        }
        else
        {
            dbscan.setData(nObs, nFeatures, originTimes.data());
        }
        // Cluster
        dbscan.cluster();
        auto nDetections = dbscan.getNumberOfClusters();
        if (nDetections < 1){continue;} // No detections 
        printf("Event %ld created %d events\n", evid, nDetections);
        auto labels = dbscan.getLabels();
        // For each cluster extract the indices
        for (int ic=0; ic<nDetections; ++ic)
        {
            TemplateIndices indices;
            indices.evid = evid;
            for (size_t i=i1; i<i2; ++i)
            {
                if (labels[i-i1] == ic){indices.indices.push_back(i);}
            }
            pImpl->mTemplateIndices.push_back(indices);
        }
    } // Loop on unique event IDs
    pImpl->mAssociations = static_cast<int> (pImpl->mTemplateIndices.size());
}

/// Gets the number of associated events
template<class T>
int Associator<T>::getNumberOfEvents() const noexcept
{
    return static_cast<int> (pImpl->mAssociations);
}

/// Gets the waveforms for a given detection
template<class T>
std::vector<MFLib::SingleChannel::Detection<T>> 
Associator<T>::getDetectionsInEvent(const int iev) const 
{
    int nDetections = getNumberOfDetections();
    if (nDetections < 1)
    {
        throw std::runtime_error("No detections have been associated\n");
    }
    if (iev < 0 || iev >= nDetections)
    {
        throw std::invalid_argument("iev = " + std::to_string(iev)
                                  + " must be in range [0,"
                                  + std::to_string(nDetections-1) + "]\n");
    }
    std::vector<MFLib::SingleChannel::Detection<T>> detections;
    detections.resize(pImpl->mTemplateIndices[iev].indices.size());
    for (int i=0; i<static_cast<int> (detections.size()); ++i)
    {
        auto idx = pImpl->mTemplateIndices[iev].indices[i];
        detections[i] = pImpl->mDetections[idx];
    }
    return detections;
}

/// Implementation
template class MFLib::SingleChannel::Associator<double>;
template class MFLib::SingleChannel::Associator<float>;
