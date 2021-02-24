/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "G4ProductionCutsTable.hh"

void init_G4ProductionCutsTable(py::module &m) {

    py::class_<G4ProductionCutsTable>(m, "G4ProductionCutsTable")
        //.def(py::init())
        .def_static("GetProductionCutsTable", &G4ProductionCutsTable::GetProductionCutsTable,
                    py::return_value_policy::reference)
        .def("SetEnergyRange", &G4ProductionCutsTable::SetEnergyRange);
}

