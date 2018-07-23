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

#include <OpenMS/ANALYSIS/XLMS/ONuXLFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLFragmentAnnotationHelper.h>

using namespace std;

namespace OpenMS
{
namespace OpenNuXL
{

void ONuXLFragmentIonGenerator::addMS2MarkerIons(
  const vector<ONuXLFragmentAdduct> &marker_ions, 
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
  for (auto const & m : marker_ions)
  {
    const double mz = m.mass + Constants::PROTON_MASS_U;

    spectrum.emplace_back(mz, 1.0);
    spectrum_charge.emplace_back(1);
    spectrum_annotation.emplace_back(ANNOTATIONS_MARKER_ION_PREFIX + m.name);  // add name (e.g., MI:U-H2O)
  }
}

void ONuXLFragmentIonGenerator::addSpecialLysImmonumIons(
  const String& unmodified_sequence,
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
   if (unmodified_sequence.has('K'))
   {
      const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight(); 
      spectrum.emplace_back(immonium_ion2_mz, 1.0);
      spectrum_charge.emplace_back(1);
      spectrum_annotation.emplace_back(String("iK(C5H10N1)"));

      // usually only observed without shift (A. Stuetzer)
      const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight(); 
      spectrum.emplace_back(immonium_ion3_mz, 1.0);
      spectrum_charge.emplace_back(1);
      spectrum_annotation.emplace_back(String("iK(C6H13N2O)"));
    }
}

void ONuXLFragmentIonGenerator::addShiftedImmoniumIons(const String &unmodified_sequence,
                                                                    const String &fragment_shift_name,
                                                                    const double fragment_shift_mass,
                                                                    PeakSpectrum &partial_loss_spectrum,
                                                                    PeakSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
                                                                    PeakSpectrum::StringDataArray &partial_loss_spectrum_annotation) 
{

  if (unmodified_sequence.hasSubstring("Y"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10NO").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('Y', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("W"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C10H11N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('W', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("F"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('F', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("H"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H8N3").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('H', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("C"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C2H6NS").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('C', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("P"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C4H8N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('P', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("L") || unmodified_sequence.hasSubstring("I"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H12N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('L', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("K"))
  {
    // classical immonium ion
    const double immonium_ion_mz = EmpiricalFormula("C5H13N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('K', fragment_shift_name));

    // TODO: check if only DNA specific and if also other shifts are observed
    // according to A. Stuetzer mainly observed with C‘-NH3 (94.0167 Da)
    const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion2_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C5H10N1)" + fragment_shift_name));

    // usually only observed without shift (A. Stuetzer)
    const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion3_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C6H13N2O)" + fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("M"))
  {
    const double immonium_ion_mz = 104.05285 + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(ONuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('M', fragment_shift_name));
  }
}

  /* 
  * Add peaks with shifts induced by the RNA/DNA:
  *   - Precursor with complete NA-oligo for charge 1..z
  *   - Partial shifts (without complete precursor adduct)
  *     - Add shifted immonium ions for charge 1 only
  *     - Add shifted b,y,a ions + precursors for charge 1..z (adding the unshifted version and performing the shift)
  */
void ONuXLFragmentIonGenerator::generatePartialLossSpectrum(const String &unmodified_sequence,
                                                                         const AASequence &fixed_and_variable_modified_peptide,
                                                                         const double &fixed_and_variable_modified_peptide_weight,
                                                                         const String &precursor_rna_adduct,
                                                                         const double &precursor_rna_weight,
                                                                         const int &precursor_charge,
                                                                         const vector<ONuXLFragmentAdduct> &partial_loss_modification,
                                                                         const TheoreticalSpectrumGenerator &partial_loss_spectrum_generator,
                                                                         PeakSpectrum &partial_loss_spectrum)
{
  partial_loss_spectrum.getIntegerDataArrays().resize(1);
  PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge = partial_loss_spectrum.getIntegerDataArrays()[0];

  partial_loss_spectrum.getStringDataArrays().resize(1); // annotation
  PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation = partial_loss_spectrum.getStringDataArrays()[0];

  // ADD: (mainly for ETD) MS2 precursor peaks of the MS1 adduct (total RNA) carrying peptide for all z <= precursor charge
  for (int charge = 1; charge <= static_cast<int>(precursor_charge); ++charge)
  {
    addPrecursorWithCompleteRNA_(fixed_and_variable_modified_peptide_weight,
                                 precursor_rna_adduct,
                                 precursor_rna_weight,
                                 charge,
                                 partial_loss_spectrum,
                                 partial_loss_spectrum_charge,
                                 partial_loss_spectrum_annotation);
  }

  for (Size i = 0; i != partial_loss_modification.size(); ++i)
  {
    // get name and mass of fragment adduct
    const String& fragment_shift_name = partial_loss_modification[i].name; // e.g. U-H2O
    const double fragment_shift_mass = partial_loss_modification[i].mass;

    // ADD: shifted immonium ion peaks of charge 1 (if the amino acid is present in the sequence)
    ONuXLFragmentIonGenerator::addShiftedImmoniumIons(
      unmodified_sequence,
      fragment_shift_name,
      fragment_shift_mass,
      partial_loss_spectrum,
      partial_loss_spectrum_charge,
      partial_loss_spectrum_annotation);

    // annotate generated a,b,y ions with fragment shift name
    PeakSpectrum shifted_series_peaks;
    shifted_series_peaks.getStringDataArrays().resize(1); // annotation
    shifted_series_peaks.getIntegerDataArrays().resize(1); // charge

    PeakSpectrum::StringDataArray& shifted_series_annotations = shifted_series_peaks.getStringDataArrays()[0];
    PeakSpectrum::IntegerDataArray& shifted_series_charges = shifted_series_peaks.getIntegerDataArrays()[0];

    // For every charge state
    for (int z = 1; z <= precursor_charge; ++z)
    {
      // 1. create unshifted peaks (a,b,y, MS2 precursor ions up to pc charge)
      PeakSpectrum tmp_shifted_series_peaks;
      partial_loss_spectrum_generator.getSpectrum(tmp_shifted_series_peaks, fixed_and_variable_modified_peptide, z, z);

      PeakSpectrum::StringDataArray& tmp_shifted_series_annotations = tmp_shifted_series_peaks.getStringDataArrays()[0];
      PeakSpectrum::IntegerDataArray& tmp_shifted_series_charges = tmp_shifted_series_peaks.getIntegerDataArrays()[0];

      // 2. shift peaks 
      for (Size i = 0; i != tmp_shifted_series_peaks.size(); ++i) 
      { 
        Peak1D& p = tmp_shifted_series_peaks[i];
        p.setMZ(p.getMZ() + fragment_shift_mass / static_cast<double>(z));         
      } 

      // 3. add shifted peaks to shifted_series_peaks
      shifted_series_peaks.insert(shifted_series_peaks.end(), tmp_shifted_series_peaks.begin(), tmp_shifted_series_peaks.end());
      shifted_series_annotations.insert(
        shifted_series_annotations.end(),
        tmp_shifted_series_annotations.begin(),
        tmp_shifted_series_annotations.end()
      );
      shifted_series_charges.insert(
        shifted_series_charges.end(),
        tmp_shifted_series_charges.begin(),
        tmp_shifted_series_charges.end()
      );
    }

    // 4. add fragment shift name to annotation of shifted peaks
    for (Size j = 0; j != shifted_series_annotations.size(); ++j)
    {
      shifted_series_annotations[j] = shifted_series_annotations[j] + " " + fragment_shift_name;
    }

    // append shifted and annotated ion series to partial loss spectrum
    partial_loss_spectrum.insert(partial_loss_spectrum.end(), shifted_series_peaks.begin(), shifted_series_peaks.end());
    partial_loss_spectrum_annotation.insert(
      partial_loss_spectrum_annotation.end(),
      shifted_series_annotations.begin(),
      shifted_series_annotations.end()
    );
    partial_loss_spectrum.getIntegerDataArrays()[0].insert(
      partial_loss_spectrum_charge.end(),
      shifted_series_charges.begin(),
      shifted_series_charges.end()
    );
  }

  partial_loss_spectrum.sortByPosition();
}

void ONuXLFragmentIonGenerator::addPrecursorWithCompleteRNA_(
  const double fixed_and_variable_modified_peptide_weight, 
  const String &precursor_rna_adduct,
  const double precursor_rna_weight, 
  const int charge, 
  PeakSpectrum &partial_loss_spectrum,
  MSSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
  MSSpectrum::StringDataArray &partial_loss_spectrum_annotation)
{
  const double xl_mz = (fixed_and_variable_modified_peptide_weight + precursor_rna_weight +
                  static_cast<double>(charge) * Constants::PROTON_MASS_U)
                 / static_cast<double>(charge);
  partial_loss_spectrum.push_back(Peak1D(xl_mz, 1.0));
  partial_loss_spectrum_charge.push_back(charge);
  partial_loss_spectrum_annotation.push_back(String("[M+") + precursor_rna_adduct + "]");
}

}
}