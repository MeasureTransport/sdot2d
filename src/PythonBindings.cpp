#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "SDOT/BoundingBox.h"
#include "SDOT/LaguerreDiagram.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/Distribution2d.h"
#include "SDOT/DiscretizedDistribution.h"
#include "SDOT/SemiDiscreteOT.h"
#include "SDOT/Distances/Distances.h"

namespace py = pybind11;
using namespace sdot;
using namespace sdot::distances;

PYBIND11_MODULE(_pysdot, m) {

  py::class_<BoundingBox, std::shared_ptr<BoundingBox>>(m, "BoundingBox")
    .def(py::init<double,double,double,double>())
    .def_readonly("xMin", &BoundingBox::xMin)
    .def_readonly("xMax", &BoundingBox::xMax)
    .def_readonly("yMin", &BoundingBox::yMin)
    .def_readonly("yMax", &BoundingBox::yMax);

  py::class_<Polygon, std::shared_ptr<Polygon>>(m,"Polygon")
    .def(py::init<Eigen::Matrix2Xd const&>())
    .def("GetVertices", &Polygon::GetVertices)
    .def("ConvexPartition",&Polygon::ConvexPartition)
    .def("IsConvex", &Polygon::IsConvex)
    .def("ConvexHull", &Polygon::ConvexHull)
    .def("Intersection", &Polygon::Intersection)
    .def("Intersects", &Polygon::Intersects)
    .def("Area",&Polygon::Area)
    .def("SecondMoment", &Polygon::SecondMoment)
    .def("Translate", &Polygon::Translate)
    .def("Rotate", &Polygon::Rotate);

  py::class_<RegularGrid, std::shared_ptr<RegularGrid>>(m,"RegularGrid")
    .def(py::init<double,double,double,double, unsigned int, unsigned int>())
    .def(py::init<BoundingBox const&, unsigned int, unsigned int>())
    .def("NumCells",&RegularGrid::NumCells, py::arg("dim") = -1)
    .def("NumNodes",&RegularGrid::NumNodes, py::arg("dim") = -1)
    .def("Center", &RegularGrid::Center)
    .def("LeftSide",&RegularGrid::LeftSide)
    .def("RightSide",&RegularGrid::RightSide)
    .def("TopSide",&RegularGrid::TopSide)
    .def("BottomSide",&RegularGrid::BottomSide)
    .def("LeftNode", &RegularGrid::LeftNode)
    .def("RightNode", &RegularGrid::RightNode)
    .def("TopNode", &RegularGrid::TopNode)
    .def("BottomNode", &RegularGrid::BottomNode)
    .def_readonly("xMin", &RegularGrid::xMin)
    .def_readonly("xMax", &RegularGrid::xMax)
    .def_readonly("yMin", &RegularGrid::yMin)
    .def_readonly("yMax", &RegularGrid::yMax)
    .def_readonly("Nx", &RegularGrid::Nx)
    .def_readonly("Ny", &RegularGrid::Ny)
    .def_readonly("dx", &RegularGrid::dx)
    .def_readonly("dy", &RegularGrid::dy);

  py::class_<Distribution2d, std::shared_ptr<Distribution2d>>(m,"Distribution2d")
    .def("Density", &Distribution2d::Density)
    .def("Grid", &Distribution2d::Grid);

  py::class_<DiscretizedDistribution, Distribution2d, std::shared_ptr<DiscretizedDistribution>>(m,"DiscretizedDistribution")
    .def(py::init<std::shared_ptr<RegularGrid> const&, Eigen::MatrixXd const&>());

  py::class_<LaguerreDiagram, std::shared_ptr<LaguerreDiagram>>(m, "LaguerreDiagram")
    .def(py::init<double, double, double, double, Eigen::Matrix2Xd const&, Eigen::VectorXd const&>())
    .def(py::init<BoundingBox const&, Eigen::Matrix2Xd const&, Eigen::VectorXd const&>())
    .def("NumCells", &LaguerreDiagram::NumCells)
    .def("SeedPts", &LaguerreDiagram::SeedPts)
    .def("Prices", &LaguerreDiagram::Prices)
    .def("GetCellVertices", &LaguerreDiagram::GetCellVertices)
    .def("CellCentroid", (Eigen::Vector2d (LaguerreDiagram::*)(unsigned int) const) &LaguerreDiagram::CellCentroid)
    .def("CellCentroid", (Eigen::Vector2d (LaguerreDiagram::*)(unsigned int, std::shared_ptr<Distribution2d> const&) const) &LaguerreDiagram::CellCentroid)
    .def("Centroids", &LaguerreDiagram::Centroids, py::arg("dist")=std::shared_ptr<Distribution2d>(nullptr))
    .def("CellArea", (double (LaguerreDiagram::*)(unsigned int) const) &LaguerreDiagram::CellArea)
    .def("CellArea", (double (LaguerreDiagram::*)(unsigned int, std::shared_ptr<Distribution2d> const&) const) &LaguerreDiagram::CellArea)
    .def("Areas", &LaguerreDiagram::Areas, py::arg("dist")=std::shared_ptr<Distribution2d>(nullptr))
    .def("BoundBox", &LaguerreDiagram::BoundBox)
    .def_static("LatinHypercubeSample", &LaguerreDiagram::LatinHypercubeSample)
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(BoundingBox const&, unsigned int, unsigned int, double, std::shared_ptr<Distribution2d> const&)) &LaguerreDiagram::BuildCentroidal, py::arg("bbox"), py::arg("numPts"), py::arg("maxIts")=200, py::arg("tol")=1e-3, py::arg("dist")=std::shared_ptr<Distribution2d>(nullptr))
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(BoundingBox const&, Eigen::Matrix2Xd const&, Eigen::VectorXd const&, unsigned int, double, std::shared_ptr<Distribution2d> const&)) &LaguerreDiagram::BuildCentroidal, py::arg("bbox"), py::arg("points"), py::arg("prices"), py::arg("maxIts")=200, py::arg("tol")=1e-3, py::arg("dist")=std::shared_ptr<Distribution2d>(nullptr));

  py::class_<SemidiscreteOT<Wasserstein2>, std::shared_ptr<SemidiscreteOT<Wasserstein2>>>(m, "SemidiscreteOT")
    .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&>())
    .def("Solve", &SemidiscreteOT<Wasserstein2>::Solve, py::arg("prices0"), py::arg("opts")=OptionList())
    .def("Diagram", &SemidiscreteOT<Wasserstein2>::Diagram)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<Wasserstein2>::PointGradient)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)() const) &SemidiscreteOT<Wasserstein2>::PointGradient)
    .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<Wasserstein2>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<Wasserstein2>::PointHessian)
    .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<Wasserstein2>::*)() const) &SemidiscreteOT<Wasserstein2>::PointHessian)
    .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<Wasserstein2>::LloydPointHessian)
    .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)() const) &SemidiscreteOT<Wasserstein2>::LloydPointHessian)
    .def("SetPoints", &SemidiscreteOT<Wasserstein2>::SetPoints)
    .def("Objective", &SemidiscreteOT<Wasserstein2>::Objective)
    .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<Wasserstein2>::MarginalCentroids)
    .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<Wasserstein2>::*)() const) &SemidiscreteOT<Wasserstein2>::MarginalCentroids)
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::Matrix2Xd const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<Wasserstein2>::BuildCentroidal, py::arg("dist"), py::arg("initialPoints"), py::arg("probs"),py::arg("opts")=OptionList())
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<Wasserstein2>::BuildCentroidal, py::arg("dist"), py::arg("probs"),py::arg("opts")=OptionList())
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, unsigned int, OptionList)) &SemidiscreteOT<Wasserstein2>::BuildCentroidal, py::arg("dist"),py::arg("numPts"), py::arg("opts")=OptionList());

  py::class_<SemidiscreteOT<QuadraticRegularization>, std::shared_ptr<SemidiscreteOT<QuadraticRegularization>>>(m, "SemidiscreteQR")
    .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&>())
    .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&, double>())
    .def("Solve", &SemidiscreteOT<QuadraticRegularization>::Solve, py::arg("prices0"), py::arg("opts")=OptionList())
    .def("Diagram", &SemidiscreteOT<QuadraticRegularization>::Diagram)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<QuadraticRegularization>::PointGradient)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)() const) &SemidiscreteOT<QuadraticRegularization>::PointGradient)
    .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<QuadraticRegularization>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<QuadraticRegularization>::PointHessian)
    .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<QuadraticRegularization>::*)() const) &SemidiscreteOT<QuadraticRegularization>::PointHessian)
    .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<QuadraticRegularization>::LloydPointHessian)
    .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)() const) &SemidiscreteOT<QuadraticRegularization>::LloydPointHessian)
    .def("SetPoints", &SemidiscreteOT<QuadraticRegularization>::SetPoints)
    .def("Objective", &SemidiscreteOT<QuadraticRegularization>::Objective)
    .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<QuadraticRegularization>::MarginalCentroids)
    .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<QuadraticRegularization>::*)() const) &SemidiscreteOT<QuadraticRegularization>::MarginalCentroids)
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::Matrix2Xd const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<QuadraticRegularization>::BuildCentroidal, py::arg("dist"), py::arg("initialPoints"), py::arg("probs"),py::arg("opts")=OptionList())
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<QuadraticRegularization>::BuildCentroidal, py::arg("dist"), py::arg("probs"),py::arg("opts")=OptionList())
    .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, unsigned int, OptionList)) &SemidiscreteOT<QuadraticRegularization>::BuildCentroidal, py::arg("dist"),py::arg("numPts"), py::arg("opts")=OptionList());

    py::class_<SemidiscreteOT<GHK>, std::shared_ptr<SemidiscreteOT<GHK>>>(m, "SemidiscreteGHK")
      .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&>())
      .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&, double>())
      .def("Solve", &SemidiscreteOT<GHK>::Solve, py::arg("prices0"), py::arg("opts")=OptionList())
      .def("Diagram", &SemidiscreteOT<GHK>::Diagram)
      .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<GHK>::PointGradient)
      .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)() const) &SemidiscreteOT<GHK>::PointGradient)
      .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<GHK>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<GHK>::PointHessian)
      .def("PointHessian", (Eigen::SparseMatrix<double> (SemidiscreteOT<GHK>::*)() const) &SemidiscreteOT<GHK>::PointHessian)
      .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<GHK>::LloydPointHessian)
      .def("LloydPointHessian", (Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)() const) &SemidiscreteOT<GHK>::LloydPointHessian)
      .def("SetPoints", &SemidiscreteOT<GHK>::SetPoints)
      .def("Objective", &SemidiscreteOT<GHK>::Objective)
      .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)(Eigen::VectorXd const&, LaguerreDiagram const&) const) &SemidiscreteOT<GHK>::MarginalCentroids)
      .def("MarginalCentroids",(Eigen::Matrix2Xd (SemidiscreteOT<GHK>::*)() const) &SemidiscreteOT<GHK>::MarginalCentroids)
      .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::Matrix2Xd const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<GHK>::BuildCentroidal, py::arg("dist"), py::arg("initialPoints"), py::arg("probs"),py::arg("opts")=OptionList())
      .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, Eigen::VectorXd const&, OptionList)) &SemidiscreteOT<GHK>::BuildCentroidal, py::arg("dist"), py::arg("probs"),py::arg("opts")=OptionList())
      .def_static("BuildCentroidal", (std::shared_ptr<LaguerreDiagram> (*)(std::shared_ptr<Distribution2d> const&, unsigned int, OptionList)) &SemidiscreteOT<GHK>::BuildCentroidal, py::arg("dist"),py::arg("numPts"), py::arg("opts")=OptionList());
}
