#include "SDOT/Distances/Distances.h"

#include <gtest/gtest.h>

using namespace sdot;
using namespace sdot::distances;

TEST(Distances, Wasserstein2)
{
  Eigen::Vector2d bottomLeft{0.0,0.0};
  Eigen::Vector2d bottomRight{1.0,0.0};
  Eigen::Vector2d topRight{1.0,1.0};
  Eigen::Vector2d topLeft{0.0,1.0};

  Eigen::Vector2d xi{0.25,0.5};
  double wi = 0.5;

  double val = Wasserstein2::RectangularIntegral(wi, xi, bottomLeft, topRight);
  double truth = wi - std::pow(xi[0],2.0) +xi[0] - std::pow(xi[1],2.0) + xi[1] - 2.0/3.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topRight);
  truth = 0.5*wi - 0.5*xi[0]*xi[0] + 2.0*xi[0]/3.0 -0.5*xi[1]*xi[1] +xi[1]/3.0 - 1.0/3.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topLeft);
  truth = 0.5*wi - 0.5*xi[0]*xi[0] + xi[0]/3.0 -0.5*xi[1]*xi[1] +xi[1]/3.0 - 1.0/6.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegralDeriv(wi,xi, bottomLeft, bottomRight, topLeft);
  EXPECT_DOUBLE_EQ(0.5, val);

  val = Wasserstein2::LineIntegralDeriv(wi, xi, bottomLeft, topRight);
  EXPECT_DOUBLE_EQ(std::sqrt(2), val);
}


TEST(Distances, QuadraticRegularization)
{
  Eigen::Vector2d bottomLeft{0.0,0.0};
  Eigen::Vector2d bottomRight{1.0,0.0};
  Eigen::Vector2d topRight{1.0,1.0};
  Eigen::Vector2d topLeft{0.0,1.0};

  Eigen::Vector2d xi{0.25,0.5};
  double wi = 0.5;

  // Numerical approximation of integral with center point rule
  int N = 5000;
  double dx = (topRight(0)-bottomLeft(0))/(N+1.0);
  double dy = (topRight(1)-bottomLeft(1))/(N+1.0);


  // Conjugate function over rectangle
  Eigen::Vector2d pt(2);
  Eigen::VectorXd vals(N*N);
  for(int xind=0; xind<N; ++xind){
    for(int yind=0; yind<N; ++yind){
      pt(0) = bottomLeft(0) + dx*(xind+0.5);
      pt(1) = bottomLeft(1) + dy*(yind+0.5);
      vals(yind+xind*N) = dx*dy*QuadraticRegularization::Evaluate(wi-(pt-xi).squaredNorm());
    }
  }
  double truth = vals.sum();

  // Call the quadratic regularization class
  double integral = QuadraticRegularization::RectangularIntegral(wi,xi,bottomLeft,topRight);

  EXPECT_NEAR(truth, integral, 2e-4);


  // First derivative volume integral over rectangle
  for(int xind=0; xind<N; ++xind){
    for(int yind=0; yind<N; ++yind){
      pt(0) = bottomLeft(0) + dx*(xind+0.5);
      pt(1) = bottomLeft(1) + dy*(yind+0.5);
      vals(yind+xind*N) = dx*dy*QuadraticRegularization::Derivative(wi-(pt-xi).squaredNorm());
    }
  }
  truth = vals.sum();

  integral = QuadraticRegularization::RectangularIntegralDeriv(wi,xi,bottomLeft,topRight);

  EXPECT_NEAR(truth, integral, 5e-4);


  // Second derivative volume integral over rectangle
  for(int xind=0; xind<N; ++xind){
    for(int yind=0; yind<N; ++yind){
      pt(0) = bottomLeft(0) + dx*(xind+0.5);
      pt(1) = bottomLeft(1) + dy*(yind+0.5);
      vals(yind+xind*N) = dx*dy*QuadraticRegularization::Derivative2(wi-(pt-xi).squaredNorm());
    }
  }
  truth = vals.sum();

  integral = QuadraticRegularization::RectangularIntegralDeriv2(wi,xi,bottomLeft,topRight);

  EXPECT_NEAR(truth, integral, 1e-3);


  // Derivative integral over triangle
  for(int xind=0; xind<N; ++xind){
    for(int yind=0; yind<N; ++yind){
      pt(0) = bottomLeft(0) + dx*(xind+0.5);
      pt(1) = bottomLeft(1) + dy*(yind+0.5);
      if(pt(1)<=1.0-pt(0)){
        vals(yind+xind*N) = dx*dy*QuadraticRegularization::Derivative(wi-(pt-xi).squaredNorm());
      }else{
        vals(yind+xind*N) = 0.0;
      }
    }
  }
  truth = vals.sum();

  integral = QuadraticRegularization::TriangularIntegralDeriv(wi,xi,bottomLeft, bottomRight, topLeft);

  EXPECT_NEAR(truth, integral, 5e-4);
}
