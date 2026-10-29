/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef PSM_STORE_H_
#define PSM_STORE_H_

#include <memory>
#include <vector>

#include "PSMDescription.h"

/*
 * PsmStore owns PSMDescription objects in fixed-size blocks. Their
 * addresses stay valid while more PSMs are added, and there is no heap
 * allocation per PSM.
 */
class PsmStore {
 public:
  PsmStore() = default;
  PsmStore(const PsmStore&) = delete;
  PsmStore& operator=(const PsmStore&) = delete;

  // moves psm into the store and returns its address in the store
  PSMDescription* add(PSMDescription&& psm) {
    PSMDescription* slot;
    if (!freeSlots_.empty()) {
      slot = freeSlots_.back();
      freeSlots_.pop_back();
    } else {
      if (numUsed_ == blocks_.size() * kBlockSize) {
        blocks_.emplace_back(new PSMDescription[kBlockSize]);
      }
      slot = &blocks_.back()[numUsed_ % kBlockSize];
      ++numUsed_;
    }
    *slot = std::move(psm);
    return slot;
  }

  // empties the slot of a PSM so that a later add can reuse it
  void release(PSMDescription* psm) {
    *psm = PSMDescription();
    freeSlots_.push_back(psm);
  }

  void clear() {
    blocks_.clear();
    freeSlots_.clear();
    numUsed_ = 0;
  }

 private:
  static const size_t kBlockSize = 1u << 16;
  std::vector<std::unique_ptr<PSMDescription[]> > blocks_;
  std::vector<PSMDescription*> freeSlots_;
  size_t numUsed_ = 0;
};

#endif /* PSM_STORE_H_ */
