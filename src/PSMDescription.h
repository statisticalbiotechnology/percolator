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
#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "Enzyme.h"

/*
 * PSMDescription
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 *
 */
class PSMDescription {
   public:
    PSMDescription() = default;
    PSMDescription(const std::string& peptide) : peptide(peptide) {}


    // TODO: move these static functions somewhere else
    static std::string removePTMs(const std::string& peptideSeq);
    static std::string removeFlanks(const std::string& peptideSeq) {
        return peptideSeq.substr(2, peptideSeq.size() - 4);
    }
    static bool ptrLess(PSMDescription* one, PSMDescription* other) {
        return *one < *other;
    }
    static bool ptrEqual(PSMDescription* one, PSMDescription* other) {
        return *one == *other;
    }

    std::string getPeptideSequence() { return peptide.substr(2, peptide.size() - 4); }
    // same as getPeptideSequence, without copying
    std::string_view getPeptideSequenceView() const {
        return std::string_view(peptide).substr(2, peptide.size() - 4);
    }
    std::string& getFullPeptideSequence() { return peptide; }
    std::string getFlankN() { return peptide.substr(0, 1); }
    std::string getFlankC() { return peptide.substr(peptide.size() - 1, peptide.size()); }

    friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
    void printProteins(std::ostream& out);

    bool operator<(const PSMDescription& other) const {
        return (peptide < other.peptide) ||
               (peptide == other.peptide && getRetentionTime() < other.getRetentionTime());
    }

    bool operator==(const PSMDescription& other) const {
        return (peptide == other.peptide);
    }

    inline void setId(std::string id) { id_ = std::move(id); }
    inline std::string& getId() { return id_; }

    std::string& getFullPeptide() { return peptide; }
    void setPeptide(std::string pep_seq) { peptide = std::move(pep_seq); }
    static inline void setProteinNameSeparator(const std::string sep) {
        proteinNameSeparator_ = sep;
    }
    static inline const std::string& getProteinNameSeparator() { return proteinNameSeparator_; }

    void setSpectrumFileName(const std::string& fileName) {
        // consecutive PSMs usually come from the same file
        size_t& lastIndex = lastSpectrumFileIndex_;
        if (lastIndex >= spectraFileNames_.size() || spectraFileNames_[lastIndex] != fileName) {
            auto specFilePos = std::find(spectraFileNames_.begin(), spectraFileNames_.end(), fileName);
            lastIndex = specFilePos - spectraFileNames_.begin();
            if (specFilePos == spectraFileNames_.end()) {
                spectraFileNames_.push_back(fileName);
            }
        }
        specFileNr = static_cast<unsigned int>(lastIndex);
    }
    // forgets the spectrum file names of an earlier run in the same process
    static void resetSpectrumFileNames() {
        spectraFileNames_.clear();
        lastSpectrumFileIndex_ = 0;
    }
    inline const std::string& getSpectrumFileName() {
        static const std::string noFileName;
        return hasSpectrumFileName() ? spectraFileNames_.at(specFileNr) : noFileName;
    }
    inline bool static hasSpectrumFileName() { return !spectraFileNames_.empty(); }

    inline void setRetentionTime(const double retentionTime) {
        retentionTime_ = retentionTime;
    }
    inline double getRetentionTime() const { return retentionTime_; }

    double* features = NULL;  // owned by a FeatureMemoryPool instance, no need to delete
    double expMass = 0.0, calcMass = 0.0;
    double retentionTime_ = std::numeric_limits<double>::quiet_NaN();
    unsigned int scan = 0u;
    unsigned int specFileNr = 0u;
    std::vector<std::string> proteinIds;

   protected:
    std::string id_;
    std::string peptide;
    static std::string proteinNameSeparator_;
    static vector<std::string> spectraFileNames_;
    static size_t lastSpectrumFileIndex_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
    out << "Peptide: " << psm.peptide << endl;
    out << "Spectrum scan number: " << psm.scan << endl;
    out << endl;
    return out;
}

#endif /*PSMDESCRIPTION_H_*/
