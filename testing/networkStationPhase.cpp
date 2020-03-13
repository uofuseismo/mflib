#include <cstdio>
#include <cstdlib>
#include <string>
#include <array>
#include "mflib/networkStationPhase.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace MFLib;
TEST(networkStationPhase, basic)
{
    NetworkStationPhase nsp;
    std::string network = " uU ";
    std::string networkRef = "UU";
    std::string station = " DuG  "; 
    std::string stationRef = "DUG";
    std::string phase = " s";
    std::string phaseRef = "S";

    nsp.setNetwork(network);
    nsp.setStation(station);
    nsp.setPhase(phase);
    EXPECT_EQ(nsp.getNetwork(), networkRef);
    EXPECT_EQ(nsp.getStation(), stationRef);
    EXPECT_EQ(nsp.getPhase(), phaseRef);
    // Test copy c'tor
    NetworkStationPhase nspCopy(nsp);
    EXPECT_TRUE(nspCopy == nsp);
    // Perturb
    nspCopy.setPhase("P");
    EXPECT_FALSE(nspCopy == nsp);
    nsp.clear();
    nspCopy.clear();
    EXPECT_EQ(nspCopy, nsp);
}
}
