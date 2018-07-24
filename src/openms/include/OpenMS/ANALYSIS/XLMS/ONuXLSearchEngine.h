// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/XLMS/ONuXLDeisotoper.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLModificationsGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLReport.h>
#include <OpenMS/ANALYSIS/XLMS/MorpheusScore.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLMarkerIonExtractor.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLFragmentAnnotationHelper.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLParameterParsing.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/HyperScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <set>

namespace OpenMS
{
namespace OpenNuXL
{



class OPENMS_DLLAPI ONuXLSearchEngine
  : public DefaultParamHandler, public ProgressLogger
{
  public:
    ONuXLSearchEngine();

    ~ONuXLSearchEngine() override;

    /// Exit codes
    enum class ExitCodes
    {
    EXECUTION_OK,
    ILLEGAL_PARAMETERS,
    INPUT_FILE_EMPTY,
    UNEXPECTED_RESULT,
    UNKNOWN_ERROR
    };

    /*
     * @brief run OpenNuXL search.
     * Note: input must contain centroided and sorted MS2 spectra
     *        fasta_db is only changed if decoy generation is requested
     */
    ExitCodes run(
      const MSExperiment& spectra, 
      std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const String& database_filename,
      std::vector<ProteinIdentification>& protein_ids,
      std::vector<PeptideIdentification>& peptide_ids,
      TextFile& csv_file);

  private:
    void updateMembers_() override;

    /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
    struct AnnotatedHit_
    {
      StringView sequence;
      SignedSize peptide_mod_index = 0; // enumeration index of the non-RNA peptide modification
      Size rna_mod_index = 0; // index of the RNA modification
      int isotope_error = 0; // wheter the hit has been matched with isotopic misassignment

      static constexpr const char UNKNOWN_NUCLEOTIDE = '?';
      char cross_linked_nucleotide = UNKNOWN_NUCLEOTIDE;
      // main score
      float score = 0;

      // total loss morpheus related subscores
      float MIC = 0;
      float err = 0;
      float Morph = 0;

      // partial loss morpheus related subscores
      float pl_MIC = 0;
      float pl_err = 0;
      float pl_Morph = 0;

      // complete TIC fraction of explained peaks
      float total_MIC = 0;

      // subscores
      float immonium_score = 0;
      float precursor_score = 0;
      float a_ion_score = 0;
      float marker_ions_score = 0;
      float partial_loss_score = 0;

      float best_localization_score = 0;
      String localization_scores;
      String best_localization;  
      std::vector<PeptideHit::PeakAnnotation> fragment_annotations;

      static bool hasBetterScore(const AnnotatedHit_& a, const AnnotatedHit_& b)
      {
        return a.score > b.score;
      }
    };

  // determine main score and sub scores of peaks without shifts
  void scoreTotalLossFragments_(const PeakSpectrum &exp_spectrum,
                                const PeakSpectrum &total_loss_spectrum,
                                double fragment_mass_tolerance,
                                bool fragment_mass_tolerance_unit_ppm,
                                const PeakSpectrum &a_ion_sub_score_spectrum,
                                const PeakSpectrum &precursor_sub_score_spectrum,
                                const PeakSpectrum &immonium_sub_score_spectrum,
                                float &total_loss_score,
                                float &tlss_MIC,
                                float &tlss_err,
                                float &tlss_Morph,
                                float &immonium_sub_score,
                                float &precursor_sub_score,
                                float &a_ion_sub_score) const;

  void scorePartialLossFragments_(const PeakSpectrum &exp_spectrum,
                                  double fragment_mass_tolerance,
                                  bool fragment_mass_tolerance_unit_ppm,
                                  const PeakSpectrum &partial_loss_spectrum_z1,
                                  const PeakSpectrum &partial_loss_spectrum_z2,
                                  const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                                  float &partial_loss_sub_score,
                                  float &marker_ions_sub_score,
                                  float &plss_MIC, float &plss_err, float &plss_Morph) const;

    // fast or all-shifts scoring mode
    bool fast_scoring_ = true;

    // if localization should be performed
    bool localization_ = false;

    // nucleotides can form cross-link
    std::set<char> can_xl_;

  /* @brief Localization step of the cross-link identification engine.
   * Given a top scoring candidate (based on total loss spectrum) it:
   *  - generates all fragment adducts based on the attached precursor adduct
   *  - annotated peaks
   *  - calculates an additive score that considers the presence or absence of evidence for a cross-linking site
   *  - the maximum score is reported
   */
  void postScoreHits_(const PeakMap& exp, 
                      std::vector<std::vector<AnnotatedHit_> >& annotated_hits, 
                      Size top_hits, 
                      const ONuXLModificationMassesResult& mm, 
                      const std::vector<ResidueModification>& fixed_modifications, 
                      const std::vector<ResidueModification>& variable_modifications, 
                      Size max_variable_mods_per_peptide, 
                      const TheoreticalSpectrumGenerator& partial_loss_spectrum_generator, 
                      double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, 
                      const ONuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts);

  /// Filter by top scoring hits, reconstruct original peptide from memory efficient structure, and add additional meta information.
  void postProcessHits_(const PeakMap& exp, 
    std::vector<std::vector<AnnotatedHit_> >& annotated_hits, 
    std::vector<ProteinIdentification>& protein_ids, 
    std::vector<PeptideIdentification>& peptide_ids, 
    Size top_hits, 
    const ONuXLModificationMassesResult& mm, 
    const std::vector<ResidueModification>& fixed_modifications, 
    const std::vector<ResidueModification>& variable_modifications, 
    Size max_variable_mods_per_peptide,
    const String& database_filename);

  void mapPrecursorMassesToScans_(const Int min_precursor_charge,
                                 const Int max_precursor_charge,
                                 const IntList &precursor_isotopes,
                                 const double small_peptide_mass_filter_threshold,
                                 const Size peptide_min_size,
                                 const PeakMap & spectra,
                                 std::multimap<double, std::pair<Size, int>> & multimap_mass_2_scan_index) const;

  void initializeSpectrumGenerators(TheoreticalSpectrumGenerator &total_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &partial_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &a_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &immonium_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &precursor_ion_sub_score_spectrum_generator) const;


}; 

}
}