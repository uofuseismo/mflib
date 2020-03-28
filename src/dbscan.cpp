#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#ifdef __ICC
  #pragma warning(push)
  #pragma warning(disable: 68 654 1125 1225)
#endif
#ifdef __GNUG__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wunused-variable"
  #pragma GCC diagnostic ignored "-Wsign-compare"
  #pragma GCC diagnostic ignored "-Wextra"
#endif
#include <daal.h>
#ifdef __ICC
  #pragma warning(pop)
#endif
#ifdef __GNUG__
  #pragma GCC diagnostic pop
#endif
#include "private/dbscan.hpp"

using namespace MFLib;
using namespace daal;
namespace dDBSCAN = daal::algorithms::dbscan;
namespace dDM = daal::data_management::interface1;

class DBSCAN::DBSCANImpl
{
public:
    DBSCANImpl()
    {
    }
    // This is a table of features.
    daal::services::SharedPtr<dDM::HomogenNumericTable<double> > mX;
    // This is a table of weights
    daal::services::SharedPtr<dDM::HomogenNumericTable<double> > mW;
    // Labels
    std::vector<int> mLabels; 
    // Size of range query in seconds
    double mEpsilon = 3;
    // Minimum number of points for each cluster
    int mMinObservations = 5;
    // The number of clusters found
    int mNumberOfClusters = 0;
    /// Determines if the weights were set 
    bool mHaveWeights = false;
    // Flag indicating that the data was set
    bool mHaveData = false;
    // Flag if the clusters are computed
    bool mHaveClusters = false;
    // Flag indiating that this is intiialized.
    bool mInitialized = false;
};

/// C'tor
DBSCAN::DBSCAN() : 
    pImpl(std::make_unique<DBSCANImpl> ())
{
}

/// Destructor
DBSCAN::~DBSCAN() = default;

/// Resets the class
void DBSCAN::clear() noexcept
{
    if (pImpl->mX){pImpl->mX->resize(1);}
    if (pImpl->mW){pImpl->mW->resize(1);}
    pImpl->mEpsilon = 3;
    pImpl->mMinObservations = 5;
    pImpl->mNumberOfClusters = 0;
    pImpl->mHaveWeights = false;
    pImpl->mHaveData = false;
    pImpl->mHaveClusters = false;
    pImpl->mInitialized = false;
}

/// Initialize
void DBSCAN::initialize(const double epsilon, const int minObservations)
{
    clear();
    if (epsilon < 0)
    {
        throw std::invalid_argument("epsilon = " + std::to_string(epsilon)
                                  + " must be non-negative\n");
    }
    if (minObservations < 1)
    {
        throw std::invalid_argument("minObservations = "
                                  + std::to_string(minObservations)
                                  + " must be positive\n");
    }
    pImpl->mEpsilon = epsilon;
    pImpl->mMinObservations = minObservations; 
    pImpl->mInitialized = true;
}

/// Sets the data
void DBSCAN::setData(const int nObservations, const int nFeatures,
                     const double *X)
{
    pImpl->mHaveData = false;
    pImpl->mHaveClusters = false;
    pImpl->mHaveWeights = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    if (nObservations < 1){throw std::invalid_argument("No observations\n");}
    if (nFeatures < 1){throw std::invalid_argument("No features\n");}
    if (X == nullptr){throw std::invalid_argument("X is NULL\n");}
    // Allocate space
    auto allocate = dDM::NumericTable::doAllocate;
    auto nf = static_cast<size_t> (nFeatures);
    auto no = static_cast<size_t> (nObservations);
    pImpl->mX = dDM::HomogenNumericTable<double>::create(nf, no, allocate);
    // Copy the data onto X
    dDM::BlockDescriptor<double> block;
    pImpl->mX->getBlockOfRows(0, no, daal::data_management::writeOnly, block);
    auto ncopy = nf*no;
    double *xPtr = block.getBlockPtr();
    std::copy(X, X+ncopy, xPtr);
    pImpl->mX->releaseBlockOfRows(block); 
    pImpl->mHaveData = true;
}

/// Sets the data with observations weights
void DBSCAN::setWeightedData(const int nObservations, const int nFeatures,
                             const double *X, const double *weights)
{
    pImpl->mHaveData = false;
    pImpl->mHaveClusters = false;
    pImpl->mHaveWeights = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    if (nObservations < 1){throw std::invalid_argument("No observations\n");}
    if (nFeatures < 1){throw std::invalid_argument("No features\n");}
    if (X == nullptr){throw std::invalid_argument("X is NULL\n");}
    bool copyWeights = false;
    if (weights != nullptr)
    {
        auto minElement = std::min_element(weights, weights+nObservations);
        if (*minElement <= 0)
        {
            throw std::invalid_argument("All weights must be positive\n");
        }
        copyWeights = true;
    }
    // Allocate space
    auto allocate = dDM::NumericTable::doAllocate;
    auto nf = static_cast<size_t> (nFeatures);
    auto no = static_cast<size_t> (nObservations);
    pImpl->mX = dDM::HomogenNumericTable<double>::create(nf, no, allocate);
    pImpl->mW = dDM::HomogenNumericTable<double>::create(1, no, allocate, 1.0);
    // Copy the data onto X
    dDM::BlockDescriptor<double> block;
    pImpl->mX->getBlockOfRows(0, no, daal::data_management::writeOnly, block);
    auto ncopy = nf*no;
    double *xPtr = block.getBlockPtr();
    std::copy(X, X+ncopy, xPtr);
    pImpl->mX->releaseBlockOfRows(block);
    // Copy the weights
    if (copyWeights)
    {
        pImpl->mW->getBlockOfRows(0, no, daal::data_management::writeOnly, block); 
        double *wPtr = block.getBlockPtr();
        std::copy(weights, weights+nObservations, wPtr);
        pImpl->mW->releaseBlockOfRows(block);
        pImpl->mHaveWeights = true;
    }
    pImpl->mHaveData = true;
}

/// Cluster the data
void DBSCAN::cluster()
{
    // Check the data is set
    pImpl->mHaveClusters = false;
    pImpl->mLabels.resize(0);
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    if (!haveData()){throw std::runtime_error("Data not yet set\n");}
    dDBSCAN::Batch<double> dbscan(pImpl->mEpsilon, pImpl->mMinObservations);
    dbscan.input.set(dDBSCAN::data,    pImpl->mX);
    if (pImpl->mHaveWeights){dbscan.input.set(dDBSCAN::weights, pImpl->mW);}
    dbscan.compute(); 
    // Save the number of clusters
    auto clusters = dbscan.getResult()->get(dDBSCAN::nClusters);
    auto nRows = clusters->getNumberOfRows();
    dDM::BlockDescriptor<int> block;
    clusters->getBlockOfRows(0, nRows, daal::data_management::readOnly, block);
    pImpl->mNumberOfClusters = block.getBlockPtr()[0];
    clusters->releaseBlockOfRows(block);
    // Save the labels
    auto labels = dbscan.getResult()->get(dDBSCAN::assignments);
    nRows = labels->getNumberOfRows();
    labels->getBlockOfRows(0, nRows, daal::data_management::readOnly, block);
    int *labelsPtr = block.getBlockPtr();
    pImpl->mLabels.resize(nRows);
    std::copy(labelsPtr, labelsPtr+nRows, pImpl->mLabels.begin());
    labels->releaseBlockOfRows(block);
    pImpl->mHaveClusters = true;
}

/// Initialized?
bool DBSCAN::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Have the data?
bool DBSCAN::haveData() const noexcept
{
    return pImpl->mHaveData;
}

/// Get the number of clusters
int DBSCAN::getNumberOfClusters() const noexcept
{
    return pImpl->mNumberOfClusters;
}

/// Get the cluster assignments
std::vector<int> DBSCAN::getLabels() const
{
    if (!haveLabels()){throw std::runtime_error("clustering not yet done\n");}
    return pImpl->mLabels;
} 

bool DBSCAN::haveLabels() const noexcept
{
    return pImpl->mHaveClusters;
}

/*
int main()
{
    int nObs = 6;
    int nFeatures = 2;
    std::vector<double> X(nObs*nFeatures);
    X[0] = 1;
    X[1] = 2;

    X[2] = 2;
    X[3] = 2;

    X[4] = 2;
    X[5] = 2;

    X[6] = 8;
    X[7] = 7;

    X[8] = 8;
    X[9] = 8;

    X[10] = 25;
    X[11] = 80;

    DBSCAN dbscan;
    double epsilon = 3;
    int minSamples = 2;
    dbscan.initialize(epsilon, minSamples);
printf("setting data...\n");
    dbscan.setData(nObs, nFeatures, X.data());
printf("clustering...\n");
    dbscan.cluster();
    auto nClusters = dbscan.getNumberOfClusters();
    printf("nClusters: %d\n", nClusters);
    auto labels = dbscan.getLabels();
for (const auto &label : labels)
{
printf("%d\n", label); 
}

    dbscan.setData(nObs-1, nFeatures, X.data());
    dbscan.cluster();
    printf("nClusters: %d\n", nClusters);

    labels = dbscan.getLabels();
    for (const auto &label : labels)
    {
        printf("%d\n", label);
    }

}
*/
