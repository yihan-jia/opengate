/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateVActor.h"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "GateActorManager.h"
#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateMultiFunctionalDetector.h"

GateVActor::GateVActor(py::dict &user_info, bool MT_ready)
    : G4VPrimitiveScorer(DictGetStr(user_info, "name")) {
  // register this actor to the global list of actors
  GateActorManager::AddActor(this);
  // MT ?
  fMultiThreadReady = MT_ready;
  // Do not work (yet) with multi-thread
  if (!fMultiThreadReady && G4Threading::IsMultithreadedApplication()) {
    std::ostringstream oss;
    oss << "Sorry, the actor '" << GetName()
        << "' cannot (yet) be used in multi-threads mode. ";
    Fatal(oss.str());
  }
}

GateVActor::~GateVActor() {}

void GateVActor::InitializeUserInput(py::dict &user_info) {
  fMotherVolumeName = DictGetStr(user_info, "attached_to_volume");
  auto op = DictGetStr(user_info, "filters_boolean_operator");
  if (op == "and")
    fOperatorIsAnd = true;
  else
    fOperatorIsAnd = false;
}

void GateVActor::AddActions(std::set<std::string> &actions) {
  fActions.insert(actions.begin(), actions.end());
}

void GateVActor::PreUserTrackingAction(const G4Track *track) {
  for (auto f : fFilters) {
    if (!f->Accept(track))
      return;
  }
}

void GateVActor::PostUserTrackingAction(const G4Track *track) {
  for (auto f : fFilters) {
    if (!f->Accept(track))
      return;
  }
}

G4bool GateVActor::ProcessHits(G4Step *step, G4TouchableHistory *) {
  /*
   In the G4 docs:

   "The second argument is a G4TouchableHistory object for the Readout geometry
   described in the next section. The second argument is NULL if Readout
   geometry is not assigned to this sensitive detector. In this method, one or
   more G4VHit objects should be constructed if the current step is meaningful
   for your detector."

   "The second argument of FillHits() method, i.e. G4TouchableHistory, is
   obsolete and not used. If user needs to define an artificial second geometry,
   use Parallel Geometries."

    => so we decide to simplify and remove "touchable" in the following.
   */

  // if the operator is AND, we perform the SteppingAction only if ALL filters
  // are true (If only one is false, we stop and return)
  if (fOperatorIsAnd) {
    for (auto f : fFilters) {
      if (!f->Accept(step))
        return true;
    }
    SteppingAction(step);
    return true;
  }
  // if the operator is OR, we accept as soon as one filter is OK
  for (auto f : fFilters) {
    if (f->Accept(step)) {
      SteppingAction(step);
      return true;
    }
  }
  return true;
}

void GateVActor::RegisterSD(G4LogicalVolume *lv) {
  // Look is a SD already exist for this LV
  auto currentSD = lv->GetSensitiveDetector();
  GateMultiFunctionalDetector *mfd;
  if (!currentSD) {
    // This is the first time a SD is set to this LV
    auto f = new GateMultiFunctionalDetector("mfd_" + lv->GetName());
    G4SDManager::GetSDMpointer()->AddNewDetector(f);
    lv->SetSensitiveDetector(f);
    mfd = f;
  } else {
    // A SD already exist, we reused it
    mfd = dynamic_cast<GateMultiFunctionalDetector *>(currentSD);
    for (auto i = 0; i < mfd->GetNumberOfPrimitives(); i++) {
      if (mfd->GetPrimitive(i)->GetName() == GetName()) {
        // In that case the actor is already registered, we skip to avoid
        // G4 exception. It happens when the LogVol has several PhysVol
        // (repeater)
        return;
      }
    }
  }
  // Register the actor to the GateMultiFunctionalDetector
  mfd->RegisterPrimitive(this);
}

// void RegisterCallBack(std::string callback_name, std::function func) {
//     std::cout << "Register callback " << callback_name << " (not yet
//     implemented)" << std::endl;
//   if (fcallBacks.count(callback_name) > 0) {
//     std::ostringstream oss;
//     oss << "You are trying to register a callback function with the name "
//         << callback_name
//         << ", but a callback with this name is already registered.";
//     FatalKeyError(oss.str());
//   } else {
//     fallBacks.insert({callback_name, func});
//   }
// }

// std::string GetOutputPathString(std::string output_type, int run_index) {
//   CallbackMap::const_iterator pos =
//   fcallBacks.find("get_output_path_string");
//
//   if (pos == fcallBacks.end()) {
//     std::ostringstream oss;
//     oss << "No callback function 'get_output_path_string' found for output "
//            "type "
//         << output_type;
//     FatalKeyError(oss.str());
//   }
//   auto func = pos->second;
//   string::path path = func(output_type, run_index);
//   return path
// }
