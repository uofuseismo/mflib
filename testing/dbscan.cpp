#include <cstdio>
#include <cstdlib>
#include <vector>
#include "private/dbscan.hpp"
#include <gtest/gtest.h>
namespace
{
TEST(dbscan, dbscan)
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

    std::vector<int> labelsRef({0, 0, 0, 1, 1, -1});

    MFLib::DBSCAN dbscan;
    double epsilon = 3;
    int minSamples = 2;
    EXPECT_NO_THROW(dbscan.initialize(epsilon, minSamples));
    EXPECT_NO_THROW(dbscan.setData(nObs, nFeatures, X.data()));
    EXPECT_NO_THROW(dbscan.cluster());
    auto nClusters = dbscan.getNumberOfClusters();
    EXPECT_EQ(nClusters, 2);
    //printf("nClusters: %d\n", nClusters);
    auto labels = dbscan.getLabels();
    EXPECT_EQ(labels.size(), labelsRef.size());
    for (int i=0; i<static_cast<int> (labels.size()); ++i)
    {
        EXPECT_EQ(labels[i], labelsRef[i]);
    }
    // Remove the outlier
    EXPECT_NO_THROW(dbscan.setData(nObs-1, nFeatures, X.data()));
    EXPECT_NO_THROW(dbscan.cluster());
    nClusters = dbscan.getNumberOfClusters();
    //printf("nClusters: %d\n", nClusters);
    EXPECT_EQ(nClusters, 2);
    labels = dbscan.getLabels();
    EXPECT_EQ(labels.size(), labelsRef.size()-1);
    for (int i=0; i<static_cast<int> (labels.size()); ++i)
    {
        EXPECT_EQ(labels[i], labelsRef[i]);
    }
}
}
