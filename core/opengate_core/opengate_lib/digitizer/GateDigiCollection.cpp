/* --------------------------------------------------
   Copyright (C): OpenGate Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateDigiCollection.h"
#include "G4Step.hh"
#include "GateDigiAttributeManager.h"
#include "GateDigiCollectionIterator.h"
#include "GateDigiCollectionsRootManager.h"

GateDigiCollection::GateDigiCollection(const std::string &collName)
    : G4VHitsCollection("", collName), fDigiCollectionName(collName) {
  fTupleId = -1;
  fDigiCollectionTitle = "Digi collection";
  fFilename = "";
  fCurrentDigiAttributeId = 0;
  fWriteToRootFlag = true;
  threadLocalData.Get().fBeginOfEventIndex = 0;
}

GateDigiCollection::~GateDigiCollection() {}

size_t GateDigiCollection::GetBeginOfEventIndex() const {
  return threadLocalData.Get().fBeginOfEventIndex;
}

void GateDigiCollection::SetBeginOfEventIndex(size_t index) {
  threadLocalData.Get().fBeginOfEventIndex = index;
}

void GateDigiCollection::SetBeginOfEventIndex() {
  SetBeginOfEventIndex(GetSize());
}

void GateDigiCollection::SetWriteToRootFlag(bool f) { fWriteToRootFlag = f; }

void GateDigiCollection::SetFilename(std::string filename) {
  fFilename = filename;
  if (fFilename.empty())
    SetWriteToRootFlag(false);
  else
    SetWriteToRootFlag(true);
}

void GateDigiCollection::InitializeDigiAttributes(
    const std::vector<std::string> &names) {
  StartInitialization();
  for (const auto &name : names)
    InitializeDigiAttribute(name);
  FinishInitialization();
}

void GateDigiCollection::InitializeDigiAttributes(
    const std::set<std::string> &names) {
  StartInitialization();
  for (const auto &name : names)
    InitializeDigiAttribute(name);
  FinishInitialization();
}

void GateDigiCollection::StartInitialization() {
  if (!fWriteToRootFlag)
    return;
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  auto id = am->DeclareNewTuple(fDigiCollectionName);
  fTupleId = id;
}

void GateDigiCollection::InitializeRootTupleForMaster() {
  if (!fWriteToRootFlag)
    return;
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  am->CreateRootTuple(this);
}

void GateDigiCollection::InitializeRootTupleForWorker() {
  if (!fWriteToRootFlag)
    return;
  // no need if not multi-thread
  if (not G4Threading::IsMultithreadedApplication())
    return;
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  am->CreateRootTuple(this);
  SetBeginOfEventIndex();
}

void GateDigiCollection::FillToRootIfNeeded(bool clear) {
  /*
      Policy :
      - can write to root or not according to the flag
      - can clear every N calls
   */
  if (!fWriteToRootFlag) {
    // need to set the index before (in case we don't clear)
    if (clear)
      Clear();
    else
      SetBeginOfEventIndex();
    return;
  }
  FillToRoot();
}

void GateDigiCollection::FillToRoot() {
  /*
   * maybe not very efficient to loop that way (row then column)
   * but I don't manage to do elsewhere
   */
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  for (size_t i = 0; i < GetSize(); i++) {
    for (auto *att : fDigiAttributes) {
      att->FillToRoot(i);
    }
    am->AddNtupleRow(fTupleId);
  }
  // required ! Cannot fill without clear
  Clear();
}

void GateDigiCollection::Clear() {
  for (auto *att : fDigiAttributes) {
    att->Clear();
  }
  SetBeginOfEventIndex(0);
}

void GateDigiCollection::Write() const {
  if (!fWriteToRootFlag)
    return;
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  am->Write(fTupleId);
}

void GateDigiCollection::Close() const {
  if (!fWriteToRootFlag)
    return;
  auto *am = GateDigiCollectionsRootManager::GetInstance();
  am->CloseFile(fTupleId);
}

void GateDigiCollection::InitializeDigiAttribute(const std::string &name) {
  if (fDigiAttributeMap.find(name) != fDigiAttributeMap.end()) {
    std::ostringstream oss;
    oss << "Error the branch named '" << name
        << "' is already initialized. Abort";
    Fatal(oss.str());
  }
  auto *att = GateDigiAttributeManager::GetInstance()->NewDigiAttribute(name);
  InitializeDigiAttribute(att);
}

void GateDigiCollection::InitializeDigiAttribute(GateVDigiAttribute *att) {
  auto name = att->GetDigiAttributeName();
  if (fDigiAttributeMap.find(name) != fDigiAttributeMap.end()) {
    std::ostringstream oss;
    oss << "Error the branch named '" << name
        << "' is already initialized. Abort";
    Fatal(oss.str());
  }
  fDigiAttributes.push_back(att);
  fDigiAttributeMap[name] = att;
  att->SetDigiAttributeId(fCurrentDigiAttributeId);
  att->SetTupleId(fTupleId);
  fCurrentDigiAttributeId++;
  // special case for type=3
  if (att->GetDigiAttributeType() == '3')
    fCurrentDigiAttributeId += 2;
}

void GateDigiCollection::FinishInitialization() {
  // Finally, not useful (yet?)
}

void GateDigiCollection::FillHits(G4Step *step) {
  for (auto *att : fDigiAttributes) {
    att->ProcessHits(step);
  }
}

void GateDigiCollection::FillDigiWithEmptyValue() {
  for (auto *att : fDigiAttributes) {
    att->FillDigiWithEmptyValue();
  }
}

size_t GateDigiCollection::GetSize() const {
  if (fDigiAttributes.empty())
    return 0;
  return fDigiAttributes[0]->GetSize();
}

GateVDigiAttribute *
GateDigiCollection::GetDigiAttribute(const std::string &name) {
  // Sometimes it is faster to apologize instead of asking permission ...
  try {
    return fDigiAttributeMap.at(name);
  } catch (std::out_of_range &) {
    std::ostringstream oss;
    oss << "Error the branch named '" << name << "' does not exist. Abort";
    Fatal(oss.str());
  }
  return nullptr; // fake to avoid warning
}

bool GateDigiCollection::IsDigiAttributeExists(const std::string &name) const {
  return (fDigiAttributeMap.count(name) != 0);
}

std::set<std::string> GateDigiCollection::GetDigiAttributeNames() const {
  std::set<std::string> list;
  for (auto *att : fDigiAttributes)
    list.insert(att->GetDigiAttributeName());
  return list;
}

GateDigiCollection::Iterator GateDigiCollection::NewIterator() {
  return {this, 0};
}

std::string GateDigiCollection::DumpLastDigi() const {
  std::ostringstream oss;
  auto n = GetSize() - 1;
  for (auto *att : fDigiAttributes) {
    oss << att->GetDigiAttributeName() << " = " << att->Dump(n) << "  ";
  }
  return oss.str();
}
