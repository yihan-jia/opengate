/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GatePrimaryFilter_h
#define GatePrimaryFilter_h

#include "GateVFilter.h"
#include <pybind11/stl.h>

namespace py = pybind11;

class GatePrimaryFilter : public GateVFilter {

public:
  GatePrimaryFilter() : GateVFilter() {}

  void InitializeUserInfo(py::dict &user_info) override;

  // To avoid gcc -Woverloaded-virtual
  // https://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
  using GateVFilter::Accept;

  bool Accept(G4Step *step) const override;
  bool Accept(const G4Track *track) const override;
  

  std::string fPolicy;
};

int IsPrimary(const G4Step *step);
int IsPrimary(const G4Track *track);

#endif // GatePrimaryFilter_h
