#include "SDOT/Distances/LineQuadrature.h"
#include "SDOT/Assert.h"

using namespace sdot::distances;


std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::Get(unsigned int degree)
{
  SDOT_ASSERT(degree>0);
  SDOT_ASSERT(degree<8);

  switch(degree) {
    case 1 : return GaussLegendre0();
    case 2 : return GaussLegendre1();
    case 3 : return GaussLegendre2();
    case 4 : return GaussLegendre3();
    case 5 : return GaussLegendre4();
    case 6 : return GaussLegendre5();
    case 7 : return GaussLegendre6();
  }
  return std::make_pair(Eigen::Matrix2Xd(), Eigen::VectorXd());
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre0()
{
  Eigen::VectorXd pts_00(1);
  Eigen::VectorXd wts_00(1);
  pts_00 << 0.500000000000000;
  wts_00 << 1.000000000000000;
  return std::make_pair(pts_00,wts_00);
}


std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre1()
{
  Eigen::VectorXd pts_01(2);
  Eigen::VectorXd wts_01(2);
  pts_01 << 0.211324865405187, 0.788675134594813;
  wts_01 << 0.500000000000000, 0.500000000000000;
  return std::make_pair(pts_01,wts_01);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre2()
{
  Eigen::VectorXd pts_02(3);
  Eigen::VectorXd wts_02(3);
  pts_02 << 0.112701665379258, 0.500000000000000, 0.887298334620741;
  wts_02 << 0.277777777777778, 0.444444444444444, 0.277777777777778;
  return std::make_pair(pts_02,wts_02);
}


std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre3()
{
  Eigen::VectorXd pts_03(4);
  Eigen::VectorXd wts_03(4);
  pts_03 << 0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026;
  wts_03 << 0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727;
  return std::make_pair(pts_03,wts_03);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre4()
{
  Eigen::VectorXd pts_04(5);
  Eigen::VectorXd wts_04(5);
  pts_04 << 0.046910077030668, 0.230765344947159, 0.500000000000000, 0.769234655052842, 0.953089922969332;
  wts_04 << 0.118463442528094, 0.239314335249683, 0.284444444444444, 0.239314335249684, 0.118463442528094;
  return std::make_pair(pts_04,wts_04);
}


std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre5()
{
  Eigen::VectorXd pts_05(6);
  Eigen::VectorXd wts_05(6);
  pts_05 << 0.033765242898424, 0.169395306766868, 0.380690406958402, 0.619309593041598, 0.830604693233133, 0.966234757101576;
  wts_05 << 0.085662246189585, 0.180380786524069, 0.233956967286345, 0.233956967286345, 0.180380786524069, 0.085662246189585;
  return std::make_pair(pts_05,wts_05);
}


std::pair<Eigen::VectorXd, Eigen::VectorXd> LineQuadrature::GaussLegendre6()
{
  Eigen::VectorXd pts_06(7);
  Eigen::VectorXd wts_06(7);
  pts_06 << 0.025446043828620, 0.129234407200303, 0.297077424311302, 0.500000000000000, 0.702922575688699, 0.870765592799698, 0.974553956171379;
  wts_06 << 0.064742483084435, 0.139852695744638, 0.190915025252559, 0.208979591836735, 0.190915025252559, 0.139852695744639, 0.064742483084434;
  return std::make_pair(pts_06,wts_06);
}
