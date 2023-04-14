/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateRBEActor.h"

void init_GateRBEActor(py::module &m) {
  py::class_<GateRBEActor, std::unique_ptr<GateRBEActor, py::nodelete>,
             GateVActor>(m, "GateRBEActor")
      .def(py::init<py::dict &>())
      .def_readwrite("cpp_numerator_image", &GateRBEActor::cpp_numerator_image)
      .def_readwrite("cpp_denominator_image",
                     &GateRBEActor::cpp_denominator_image)
      .def_readwrite("fPhysicalVolumeName", &GateRBEActor::fPhysicalVolumeName);
}
