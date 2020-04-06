#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "SDOT/BoundingBox.h"
#include "SDOT/LaguerreDiagram.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/Distribution2d.h"
#include "SDOT/DiscretizedDistribution.h"
#include "SDOT/SemidiscreteOT.h"

namespace py = pybind11;
using namespace sdot;


PYBIND11_MODULE(_pysdot, m) {

  py::class_<BoundingBox, std::shared_ptr<BoundingBox>>(m, "BoundingBox")
    .def(py::init<double,double,double,double>())
    .def_readonly("xMin", &BoundingBox::xMin)
    .def_readonly("xMax", &BoundingBox::xMax)
    .def_readonly("yMin", &BoundingBox::yMin)
    .def_readonly("yMax", &BoundingBox::yMax);

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

  py::class_<SemidiscreteOT, std::shared_ptr<SemidiscreteOT>>(m, "SemidiscreteOT")
    .def(py::init<std::shared_ptr<Distribution2d> const&, Eigen::MatrixXd const&, Eigen::VectorXd const&>())
    .def("Solve", &SemidiscreteOT::Solve)
    .def("Diagram", &SemidiscreteOT::Diagram)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT::*)(LaguerreDiagram const&) const) &SemidiscreteOT::PointGradient)
    .def("PointGradient", (Eigen::Matrix2Xd (SemidiscreteOT::*)() const) &SemidiscreteOT::PointGradient);
}
