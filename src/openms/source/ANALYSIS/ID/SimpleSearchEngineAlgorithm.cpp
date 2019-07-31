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

#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/METADATA/SpectrumSettings.h>

#include <map>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif


using namespace std;

namespace OpenMS
{
  SimpleSearchEngineAlgorithm::SimpleSearchEngineAlgorithm() :
    DefaultParamHandler("SimpleSearchEngineAlgorithm"),
    ProgressLogger()
  {
    defaults_.setValue("precursor:mass_tolerance", 10.0, "Width of precursor mass tolerance window");

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered.");

    defaults_.setSectionDescription("precursor", "Precursor (Parent Ion) Options");

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    defaults_.setValue("precursor:isotopes", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)");

    defaults_.setValue("fragment:mass_tolerance", 10.0, "Fragment mass tolerance");

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("modifications:fixed", ListUtils::create<String>("Carbamidomethyl (C)", ','), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", all_mods);
    defaults_.setValue("modifications:variable", ListUtils::create<String>("Oxidation (M)", ','), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", all_mods);
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    defaults_.setValue("enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("enzyme", all_enzymes);

    defaults_.setValue("peptide:min_size", 7, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");
    defaults_.setValue("peptide:motif", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.");
    defaults_.setSectionDescription("peptide", "Peptide Options");

    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setSectionDescription("report", "Reporting Options");

    defaultsToParam_();

    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural20");
    size_t index(0);
    for (auto & r : aas)
    {
      mapResidue2Index_[r->getOneLetterCode()[0]] = index;
      ++index;
    }
  }

  void SimpleSearchEngineAlgorithm::updateMembers_()
  {
    precursor_mass_tolerance_ = param_.getValue("precursor:mass_tolerance");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit");

    precursor_min_charge_ = param_.getValue("precursor:min_charge");
    precursor_max_charge_ = param_.getValue("precursor:max_charge");

    precursor_isotopes_ = param_.getValue("precursor:isotopes");

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");

    fragment_mass_tolerance_unit_ = param_.getValue("fragment:mass_tolerance_unit");

    modifications_fixed_ = param_.getValue("modifications:fixed");

    modifications_variable_ = param_.getValue("modifications:variable");

    modifications_max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    enzyme_ = param_.getValue("enzyme");

    peptide_min_size_ = param_.getValue("peptide:min_size");
    peptide_max_size_ = param_.getValue("peptide:max_size");
    peptide_missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_motif_ = param_.getValue("peptide:motif");

    report_top_hits_ = param_.getValue("report:top_hits");
  }

  // static
  void SimpleSearchEngineAlgorithm::preprocessResidueEvidence_(PeakMap& exp,
    vector<ResidueEvidenceMatrix>& rems, 
    vector<CumScoreHistogram>& cums)
  {
    rems.clear();
    cums.clear();
    rems = vector<ResidueEvidenceMatrix>(exp.size(), ResidueEvidenceMatrix());
    cums = vector<CumScoreHistogram>(exp.size(), CumScoreHistogram());

    // sort by rt
    exp.sortSpectra(false);

    // remove 0 or negative intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    // 2. replace by square root
    SqrtMower sqrt_mower_filter;
    sqrt_mower_filter.filterPeakMap(exp);

    // 3. filter smaller 5%
    Normalizer normalizer;
    auto p_n = normalizer.getDefaults();
    p_n.setValue("method", "to_one");
    normalizer.setParameters(p_n);
    normalizer.filterPeakMap(exp);
    auto p_t = threshold_mower_filter.getDefaults();
    p_t.setValue("threshold", 0.05);
    threshold_mower_filter.setParameters(p_t);
    threshold_mower_filter.filterPeakMap(exp);
    
    // determine maximum mass bin and remove peaks outside of possible range
    double max_precursorMass(0);    
    for (auto & spectrum : exp)
    {
      if (spectrum.getPrecursors().empty()) continue;
      const double precurMz = spectrum.getPrecursors()[0].getMZ();
      const int precurCharge = spectrum.getPrecursors()[0].getCharge();
      const double precursorMass = precurMz * precurCharge - precurCharge * Constants::PROTON_MASS_U;

      // there should be no peaks right to neutral precursor mass (+ 50.0 as a safeguard)
      const double experimentalMassCutoff = precursorMass + 50.0;
      spectrum.erase(spectrum.MZBegin(experimentalMassCutoff), spectrum.end());
      
      if (precursorMass > max_precursorMass) max_precursorMass = precursorMass;
    }
    size_t maxPrecurMassBin = mass2bin_(max_precursorMass + 50.0); // no peak beyond this m/z after filtering

    // TODO: optionally remove precursor peaks? disabled by default in tide but in publication on residue evidence score

    // 4. divide spectrum into 10 segments, normalize intensities to 50
    for (auto & spec : exp)
    {
      if (spec.size() < 2) { spec.clear(false); continue; };
      double min_mz(spec.front().getMZ());
      double max_mz(spec.back().getMZ());
      double segment_mz = (max_mz - min_mz) / 10.0;
      if (segment_mz == 0) continue;
      for (double i = 0.0; i < 10.0; i += 1.0)
      { 
        auto left = spec.MZBegin(min_mz + i * segment_mz);
        auto right = spec.MZEnd(min_mz + (i + 1) * segment_mz);

        auto max_intensity_it = std::max_element(left, right, 
          [](const Peak1D& a, const Peak1D& b) 
          { 
            return a.getIntensity() < b.getIntensity(); 
          });
        double max_intensity = max_intensity_it->getIntensity();
        std::transform(left, right, left, [max_intensity](Peak1D& p) { p.setIntensity(p.getIntensity() / max_intensity * 50.0); return p; });        
      }    

      // calculate ranks and assign intensity/ranks as new intensities

      // initialize original index locations
      vector<size_t> idx(spec.size());
      std::iota(idx.begin(), idx.end(), 1);

      // sort indexes based on comparing intensity values (1 = highest intensity)
      sort(idx.begin(), idx.end(),
        [&spec](size_t i1, size_t i2) { return spec[i1].getIntensity() > spec[i2].getIntensity(); });
            
      // calculate 1.0/rank intensity
      for (Size i = 0; i != spec.size(); ++i) { spec[i].setIntensity(1.0 / static_cast<double>(idx[i])); }
    }


    // TODO: depends on terminal mods
    double cTermMass = Residue::getInternalToCTerm().getMonoWeight();
    double nTermMass = Residue::getInternalToNTerm().getMonoWeight();
    OPENMS_LOG_DEBUG << "nTermMass: " << nTermMass << endl;
    OPENMS_LOG_DEBUG << "cTermMass: " << cTermMass << endl;

    // determine which bin each amino acid mass is in
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural20");
    vector<int> aaMassBin;
    vector<double> aaMass;
    for (const auto& r : aas)
    {
      const float mass = r->getMonoWeight(Residue::Internal); // TODO: check if internal is correct
      aaMass.push_back(mass);
      int binMass = (int)floor(mass2bin_(mass));  
      aaMassBin.push_back(binMass);
    }

    for (size_t spectrum_index = 0; spectrum_index != exp.size(); ++spectrum_index)
    {
      MSSpectrum& spectrum = exp[spectrum_index];
      const double precurMz = spectrum.getPrecursors()[0].getMZ(); 
      const int precurCharge = spectrum.getPrecursors()[0].getCharge();
      const double precursorMass = precurMz * precurCharge - precurCharge * Constants::PROTON_MASS_U;

      // calculate residue evidence matrix rem
      // rem[row][column]: rows are amino acids (or modified amino acids), columns are mass bins 
      ResidueEvidenceMatrix rem = ResidueEvidenceMatrix(ResidueDB::getInstance()->getResidues("Natural20").size(), vector<double>(maxPrecurMassBin, 0));
      createResidueEvidenceMatrix(spectrum, 0.02, maxPrecurMassBin, rem); // 0.02 = default fragment mass tolerance for high-res
      
      // calculate s_max: an upper bound for the maximum score to limit the number of rows of the count matrix
      // 1. maximum number of amino acids q is bounded by ceil(precursor_mass / lightest amino acid)
      // 2. maximum possible score is then the sum of the q largest column maxima
      size_t s_max(0);
      {
        vector<double> column_maxima(maxPrecurMassBin, 0);
        for (Size c = 0; c != maxPrecurMassBin; ++c) // for all mass bins (columns)
        {
          double column_maximum(0);
          for (Size r = 0; r != rem.size(); ++r) // for all entries in current mass bin column c
          {
            if (rem[r][c] > column_maximum) column_maximum = rem[r][c];
          }
          column_maxima[c] = column_maximum;
        }
        sort(column_maxima.begin(), column_maxima.end(), greater<double>()); // sort decreasing

        double q = ceil(precursorMass / 71.03712); // Alanina: TODO: potentially there is a lighter amino acid if there is a modification with negative delta mass
        s_max = ceil(accumulate(column_maxima.begin(), column_maxima.begin() + q, 0.0));
        OPENMS_LOG_DEBUG << "Maximum evidence score theoretically achievable for this spectrum: " << s_max << endl;
      }

      // precalculate histogram cumsumsC_sm of peptide counts (y-axis) with equal or better score than s (x-axis) that match to the precursor mass
      vector<double> cumsumsC_sm(s_max + 1, 0); // index = 0..max score, value = peptide count
      {
        // calculate count matrix C using dynamic programming
        // C[row][column]: rows are discretized scores s, columns are mass bins
        vector<vector<double>> C(s_max + 1, vector<double>(maxPrecurMassBin, 0));
        C[0][mass2bin_(nTermMass)] = 1; // count the N-terminus once
        for (size_t m = mass2bin_(nTermMass); m <= mass2bin_(precursorMass); ++m)
        {
          for (size_t s = 0; s <= s_max; ++s)
          {
            // to determine the number of peptides with score s at mass m 
            // we sum up all counts corresponding to cells that lead to C[s][m] by adding one of the amino acids
            // while ensuring that the score adds up to s.
            // Counts are weighted by the amino acid frequency aaP[a]
            size_t left_sum(0);
            for (size_t a = 0; a != aaMassBin.size(); ++a)
            {
              size_t a_bin = aaMassBin[a];
              size_t R_am = rem[a][m];  // residue evidence that a peptide(-prefix) with mass m ended with a
              if ((int)s - (int)R_am < 0 || (int)m - (int)a_bin < 0) continue;
              left_sum += C[s - R_am][m - a_bin] /* * aaP[a]*/; // TODO aa frequencies
            }
            C[s][m] += left_sum; // += here because we don't want to overwrite the only non-zero entry C[0][mass2bin_(nTermMass)] right from the start
          }
        }
      
        /* 
        for (auto & r : C)
        {
          String line;
          bool only_zeros = true;
          for (auto & c : r)
          {
            if (c != 0) only_zeros = false;
            line += String(c) + " ";
          }
          if (!only_zeros) cout << line << endl;
        }
        */

        // calculate (cumulative) sum of counts for precursor mass with score larger than s once
        size_t m = mass2bin_(precursorMass);
        OPENMS_LOG_DEBUG << "Precursor mass (neutral) bin: " << m << endl;
        double cumsum(0);
        for (int s = s_max; s >= 0; --s) 
        { 
          cumsum += C[s][m]; 
          cumsumsC_sm[s] = cumsum;
        }

        // output cumulative count distribution
        // for (int s = 0; s <= s_max; ++s)  { cout << "Score: " << s << "\t" << cumsumsC_sm[s] << endl; }

        // cumsumsC_sm[0] holds the cumulative sum for all scores ([0] means total count for all scores > 0) - (denominator in eq. 10)
        OPENMS_LOG_DEBUG << "Sum of C_{s,m}=" << cumsumsC_sm[0] << std::endl;

        cums[spectrum_index].swap(cumsumsC_sm);
        rems[spectrum_index].swap(rem);
      }
    }    
  }

  size_t SimpleSearchEngineAlgorithm::calculateRawResEv_(const AASequence& candidate, const ResidueEvidenceMatrix& rem) const
  {

    double nTermMass = Residue::getInternalToNTerm().getMonoWeight(); //TODO: might depend on modification

    // calculate residue evidence score
    size_t res_ev_score(0);
    double m(nTermMass);
    for (Size i = 0; i != candidate.size(); ++i)
    {
      m += candidate[i].getMonoWeight(Residue::Internal);
      size_t m_bin = mass2bin_(m);
      size_t aa_index = mapResidue2Index_.at(candidate[i].getOneLetterCode()[0]); 
      res_ev_score += rem[aa_index][m_bin];
    }
    return res_ev_score;
  }

  double SimpleSearchEngineAlgorithm::calculatePValue_(
    const AASequence& candidate, 
    const ResidueEvidenceMatrix& rem, 
    const vector<double>& cumsumsC_sm) const
  {
    const size_t res_ev_score = calculateRawResEv_(candidate, rem);

    // number of peptides with score equal or better then residue evidence score of candidate / total peptide count
    const double p_value = cumsumsC_sm[res_ev_score] / cumsumsC_sm[0];
    return p_value;
  }

  unsigned int SimpleSearchEngineAlgorithm::mass2bin_(double mass, int charge /*=1*/) 
  {
    const double bin_width = 1.0005079;
    const double bin_offset = 0.68;

    return (unsigned int)((mass + (charge - 1) * Constants::PROTON_MASS_U) / (charge*bin_width) + 1.0 - bin_offset);
  }

  double SimpleSearchEngineAlgorithm::bin2mass_(int bin, int charge /*=1*/) 
  {
    const double bin_width = 1.0005079;
    const double bin_offset = 0.68;
    return (bin - 1.0 + bin_offset) * charge * bin_width + (charge - 1) * Constants::PROTON_MASS_U;
  }


// adapted from Crux/Tide Andy Lin
void SimpleSearchEngineAlgorithm::createResidueEvidenceMatrix(
    const MSSpectrum& spectrum,
    double fragment_tolerance_Da,
    size_t max_precursor_mass_bin,
    vector<vector<double> >& residueEvidenceMatrix) 
  {
    const double precurMz = spectrum.getPrecursors()[0].getMZ(); 
    const int precurCharge = spectrum.getPrecursors()[0].getCharge();
    const double precursorMass = precurMz * precurCharge - precurCharge * Constants::PROTON_MASS_U;
    const int granularityScale = 25;

    // TODO: these are different for fixed modified termini
    double cTermMass = Residue::getInternalToCTerm().getMonoWeight();
    double nTermMass = Residue::getInternalToNTerm().getMonoWeight();

    // determine which bin each amino acid mass is in
    vector<int> aaMassBin;
    vector<double> aaMass;
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural20");
    for (const auto& r : aas)
    {
      const float mass = r->getMonoWeight(Residue::Internal); // TODO: check if internal is correct
      aaMass.push_back(mass);
      int binMass = (int)floor(mass2bin_(mass));  
      aaMassBin.push_back(binMass); 
    }

  // Need to add lines for matlab line 56?
  // Basically bounds the ion masses and ion mass bins so that
  // 1: ionMassBinB >=1
  // 2: ionMassBinB <= maxMassBin - max(aaMassBin)
  // Also need to figure out maxMassBin (currently hard coded to 8500)
  // need to do 4 times (+1/+2/B ion/ yion)

  // Populate ionMass and ionMassBin for b ions in 1+ charge state
  vector<double> ionMass;
  vector<int> ionMassBin;
  vector<double> ionIntens;

  ionMass.push_back(nTermMass);
  ionMassBin.push_back(mass2bin_(nTermMass));
  ionIntens.push_back(0.0);

  for (size_t ion = 0; ion < spectrum.size(); ++ion) 
  {
    double tmpIonMass = spectrum[ion].getMZ();
    int binTmpIonMass = (int)floor(mass2bin_(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
    ionIntens.push_back(spectrum[ion].getIntensity());
  }
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  addEvidToResEvMatrix(ionMass, ionMassBin, ionIntens,
                       max_precursor_mass_bin, aas.size(), aaMass, aaMassBin,
                       fragment_tolerance_Da, residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();


  // Find pairs of y ions in 1+ charge state
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  for (size_t ion = 0; ion < spectrum.size(); ++ion) 
  {
    // Convert to equivalent b ion masses for ease of processing
    double tmpIonMass = precursorMass - spectrum[ion].getMZ() + (2.0 * Constants::PROTON_MASS_U);
    int binTmpIonMass = (int)floor(mass2bin_(tmpIonMass));
    
    if (tmpIonMass > 0) 
    {
      ionMass.push_back(tmpIonMass);
      ionMassBin.push_back(binTmpIonMass);
      ionIntens.push_back(spectrum[ion].getIntensity());
    }
  }
  ionMass.push_back(nTermMass);
  ionMassBin.push_back(mass2bin_(nTermMass));
  ionIntens.push_back(0.0);

  reverse(ionMass.begin(), ionMass.end());
  reverse(ionMassBin.begin(), ionMassBin.end());
  reverse(ionIntens.begin(), ionIntens.end());

  addEvidToResEvMatrix(ionMass, ionMassBin, ionIntens,
                       max_precursor_mass_bin, aas.size(), aaMass, aaMassBin,
                       fragment_tolerance_Da, residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();

  // Assuming fragment ion peaks are 2+ charge only works if precusor mass charge is greater than 1.
  if (precurCharge != 1) 
  {
    // Find pairs of b ions in 2+ charge state
    ionMass.push_back(nTermMass);
    ionMassBin.push_back(mass2bin_(nTermMass));
    ionIntens.push_back(0.0);

    for (size_t ion = 0; ion < spectrum.size(); ++ion) 
    {
      double tmpIonMass = 2.0 * spectrum[ion].getMZ() - Constants::PROTON_MASS_U;
      int binTmpIonMass = (int)floor(mass2bin_(tmpIonMass));
    
      ionMass.push_back(tmpIonMass);
      ionMassBin.push_back(binTmpIonMass);
      ionIntens.push_back(spectrum[ion].getIntensity());
    }
    ionMass.push_back(precursorMass - cTermMass);
    ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
    ionIntens.push_back(0.0);

    addEvidToResEvMatrix(ionMass, ionMassBin, ionIntens,
                         max_precursor_mass_bin, aas.size(), aaMass, aaMassBin,
                         fragment_tolerance_Da, residueEvidenceMatrix);
    ionMass.clear();
    ionMassBin.clear();
    ionIntens.clear();

    // Find pairs of y ions in 2+ charge state
    ionMass.push_back(precursorMass - cTermMass);
    ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
    ionIntens.push_back(0.0);

    for (size_t ion = 0; ion < spectrum.size(); ++ion) 
    {
      double tmpIonMass = precursorMass - (2.0 * spectrum[ion].getMZ() - Constants::PROTON_MASS_U) + (2.0 * Constants::PROTON_MASS_U);
      int binTmpIonMass = (int)floor(mass2bin_(tmpIonMass));

      if (tmpIonMass > 0.0) 
      {
        ionMass.push_back(tmpIonMass);
        ionMassBin.push_back(binTmpIonMass);
        ionIntens.push_back(spectrum[ion].getIntensity());
      }
    }
    ionMass.push_back(nTermMass);
    ionMassBin.push_back(mass2bin_(nTermMass));
    ionIntens.push_back(0.0);

    reverse(ionMass.begin(), ionMass.end());
    reverse(ionMassBin.begin(),ionMassBin.end());
    reverse(ionIntens.begin(), ionIntens.end());

    addEvidToResEvMatrix(ionMass, ionMassBin, ionIntens,
                         max_precursor_mass_bin, aas.size(), aaMass, aaMassBin,
                         fragment_tolerance_Da, residueEvidenceMatrix);
    ionMass.clear();
    ionMassBin.clear();
    ionIntens.clear();
  }

  // Get maxEvidence value
  double maxEvidence = -1.0;
  for (size_t i = 0; i < max_precursor_mass_bin; ++i) 
  {
    for (size_t curAaMass = 0; curAaMass < aas.size(); ++curAaMass) 
    {
      if (residueEvidenceMatrix[curAaMass][i] > maxEvidence) { maxEvidence = residueEvidenceMatrix[curAaMass][i]; }
    }
  }

  // Discretize residue evidence so largest value is residueEvidenceIntScale
  double residueEvidenceIntScale = (double)granularityScale;
  for (int i = 0; i < max_precursor_mass_bin; i++) 
  {
    for (size_t curAaMass = 0; curAaMass < aas.size(); ++curAaMass) 
    {
      if (residueEvidenceMatrix[curAaMass][i] > 0) 
      {
        double residueEvidence = residueEvidenceMatrix[curAaMass][i];
        residueEvidenceMatrix[curAaMass][i] = round(residueEvidenceIntScale * residueEvidence / maxEvidence);
      }
    }
  }
}

// Written by Andy Lin in Feb 2018
// Helper function for CreateResidueEvidenceMatrix
void SimpleSearchEngineAlgorithm::addEvidToResEvMatrix(
  vector<double>& ionMass,
  vector<int>& ionMassBin,
  vector<double>& ionIntens,
  int maxPrecurMassBin,
  int nAA,
  const vector<double>& aaMass,
  const vector<int>& aaMassBin,
  const double fragment_tolerance_Da,
  vector<vector<double> >& residueEvidenceMatrix
  ) 
  {
    for (size_t ion = 0; ion < ionMass.size(); ++ion) 
    {
      // bIonMass is named correctly because we assume that all y ions have been converted to their b ion equivalent
      // before ionMass is passed in
      double bIonMass1 = ionMass[ion];
      int bIonMassBin1 = ionMassBin[ion];
      double bIonIntensity1 =  ionIntens[ion];
    
      for (size_t curAaMass = 0; curAaMass < nAA; ++curAaMass) 
      {
        int newResMassBin = bIonMassBin1 + aaMassBin[curAaMass];
        
        // Find all ion mass bins that match newResMassBin
        int index = find(ionMassBin.begin(), ionMassBin.end(), newResMassBin) - ionMassBin.begin();
        double score = 0.0;
        for (int i = index; i < ionMassBin.size(); ++i) // TODO: check if this line is also valid
        {
          if (newResMassBin != ionMassBin[i]) { continue; } 
          const double bIonMass2 = ionMass[i];
          const double ionMassDiff =  bIonMass2 - bIonMass1;
          const double bIonIntensity2 = ionIntens[i];

          double aaTolScore = 1.0 - (std::abs(ionMassDiff - aaMass[curAaMass]) / fragment_tolerance_Da);

          if (aaTolScore > 0.0) // in tolerance window?
          {         
            const double tmpScore = aaTolScore * (bIonIntensity1 + bIonIntensity2);  // rank score
            if (tmpScore > score) { score = tmpScore; }
          }
        }

        // Add evidence to matrix
        // Use -1 since all mass bins are index 1 instead of index 0
        // Bounds checks. Only add score if smaller than precursor mass bin.
        // When assuming each fragment peak is a 2+ charge, it is possible
        // to have a fragment peak larger than precursor mass (ie why
        // bounds check is needed).
        if (newResMassBin <= maxPrecurMassBin) 
        {
          residueEvidenceMatrix[curAaMass][newResMassBin-1] += score;
        }
      }
    }
  }

  // static
  void SimpleSearchEngineAlgorithm::preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // sort by mz
      exp[exp_index].sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(exp[exp_index], 
        fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
        1, 3,   // min / max charge 
        false,  // keep only deisotoped
        3, 10,  // min / max isopeaks 
        true);  // convert fragment m/z to mono-charge

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
    }
  }


  // static
void SimpleSearchEngineAlgorithm::postProcessHits_(const PeakMap& exp, 
      std::vector<std::vector<SimpleSearchEngineAlgorithm::AnnotatedHit_> >& annotated_hits, 
      std::vector<ProteinIdentification>& protein_ids, 
      std::vector<PeptideIdentification>& peptide_ids, 
      Size top_hits,
      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
      Size max_variable_mods_per_peptide,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const String& precursor_mass_tolerance_unit_ppm,
      const String& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const String& enzyme,
      const String& database_name)
  {
    // remove all but top n scoring
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
      annotated_hits.shrink_to_fit();
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        pi.setScoreType("hyperscore");
        pi.setHigherScoreBetter(true);
        pi.setRT(exp[scan_index].getRT());
        pi.setMZ(exp[scan_index].getPrecursors()[0].getMZ());
        Size charge = exp[scan_index].getPrecursors()[0].getCharge();

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        for (vector<AnnotatedHit_>::const_iterator a_it = annotated_hits[scan_index].begin(); a_it != annotated_hits[scan_index].end(); ++a_it)
        {
          PeptideHit ph;
          ph.setCharge(charge);

          // get unmodified string
          AASequence aas = AASequence::fromString(a_it->sequence.getString());

          // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
          vector<AASequence> all_modified_peptides;
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

          // reannotate much more memory heavy AASequence object
          AASequence fixed_and_variable_modified_peptide = all_modified_peptides[a_it->peptide_mod_index]; 
          ph.setScore(a_it->score);
          ph.setSequence(fixed_and_variable_modified_peptide);
          phs.push_back(ph);
        }
        pi.setHits(phs);
        pi.assignRanks();

#ifdef _OPENMP
#pragma omp critical (peptide_ids_access)
#endif
        {
          //clang-tidy: seems to be a false-positive in combination with omp
          peptide_ids.push_back(std::move(pi));
        }
      }
    }

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("SimpleSearchEngine");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges = String(precursor_min_charge) + ":" + String(precursor_max_charge);

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = modifications_fixed;
    search_parameters.variable_modifications = modifications_variable;
    search_parameters.missed_cleavages = peptide_missed_cleavages;
    search_parameters.fragment_mass_tolerance = fragment_mass_tolerance;
    search_parameters.precursor_mass_tolerance = precursor_mass_tolerance;
    search_parameters.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);
    protein_ids[0].setSearchParameters(std::move(search_parameters));
  }

  SimpleSearchEngineAlgorithm::ExitCodes SimpleSearchEngineAlgorithm::search(const String& in_mzML, const String& in_db, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids) const
  {
    boost::regex peptide_motif_regex(peptide_motif_);

    bool precursor_mass_tolerance_unit_ppm = (precursor_mass_tolerance_unit_ == "ppm");
    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    set<String> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());

    if (fixed_unique.size() != modifications_fixed_.size())
    {
      cout << "duplicate fixed modification provided." << endl;
      return ExitCodes::ILLEGAL_PARAMETERS;
    }

    set<String> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      cout << "duplicate variable modification provided." << endl;
      return ExitCodes::ILLEGAL_PARAMETERS;
    }

    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
    
    // load MS2 map
    PeakMap spectra;
    MzMLFile f;
    //f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzML, spectra);
    spectra.sortSpectra(true);

    // residue evidence matrix and cumulative peptide vs score counts
    std::vector<ResidueEvidenceMatrix> rems;
    std::vector<CumScoreHistogram> cums;
    preprocessResidueEvidence_(spectra, rems, cums);
/* 
    startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
    endProgress();
*/

    // build multimap of precursor mass to scan index
    multimap<double, Size> multimap_mass_2_scan_index;
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size_)
      {
        Size precursor_charge = precursor[0].getCharge();

        if (precursor_charge < precursor_min_charge_ 
         || precursor_charge > precursor_max_charge_)
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();

        // calculate precursor mass (optionally corrected for misassignment) and map it to MS scan index
        for (int isotope_number : precursor_isotopes_)
        {
          double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

          // correct for monoisotopic misassignments of the precursor annotation
          if (isotope_number != 0) { precursor_mass -= isotope_number * Constants::C13C12_MASSDIFF_U; }

          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        }
      }
    }

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param(spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    spectrum_generator.setParameters(param);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(2 * report_top_hits_); }

#ifdef _OPENMP
    // we want to do locking at the spectrum level so we get good parallelisation 
    vector<omp_lock_t> annotated_hits_lock(annotated_hits.size());
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_init_lock(&(annotated_hits_lock[i])); }
#endif

    startProgress(0, 1, "Load database from FASTA file...");
    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile::load(in_db, fasta_db);
    endProgress();

    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme_);
    digestor.setMissedCleavages(peptide_missed_cleavages_);

    startProgress(0, fasta_db.size(), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    set<StringView> processed_petides;

    std::atomic<Size> count_proteins(0), count_peptides(0);

    #pragma omp parallel for schedule(static)
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
      ++count_proteins;

      IF_MASTERTHREAD
      {
        setProgress(count_proteins);
      }

      vector<StringView> current_digest;
      digestor.digestUnmodified(fasta_db[fasta_index].sequence, current_digest, peptide_min_size_, peptide_max_size_);

      for (auto const & c : current_digest)
      { 
        const String current_peptide = c.getString();
        if (current_peptide.find_first_of("XBZ") != std::string::npos) { continue; }

        // if a peptide motif is provided skip all peptides without match
        if (!peptide_motif_.empty() && !boost::regex_match(current_peptide, peptide_motif_regex)) { continue; }          
      
        bool already_processed = false;
        #pragma omp critical (processed_peptides_access)
        {
          // peptide (and all modified variants) already processed so skip it
          if (processed_petides.find(c) != processed_petides.end())
          {
            already_processed = true;
          }
          else
          {
            processed_petides.insert(c);
          }
        }

        // skip peptides that have already been processed
        if (already_processed) { continue; }

        ++count_peptides;

        vector<AASequence> all_modified_peptides;

        // this critial section is because ResidueDB is not thread safe and new residues are created based on the PTMs
        #pragma omp critical (residuedb_access)
        {
          AASequence aas = AASequence::fromString(current_peptide);
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, modifications_max_variable_mods_per_peptide_, all_modified_peptides);
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass = candidate.getMonoWeight();

          // determine MS2 precursors that match to the current peptide mass
          multimap<double, Size>::const_iterator low_it;
          multimap<double, Size>::const_iterator up_it;

          if (precursor_mass_tolerance_unit_ppm) // ppm
          {
            low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * current_peptide_mass * precursor_mass_tolerance_ * 1e-6);
            up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * current_peptide_mass * precursor_mass_tolerance_ * 1e-6);
          }
          else // Dalton
          {
            low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * precursor_mass_tolerance_);
            up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * precursor_mass_tolerance_);
          }

          // no matching precursor in data
          if (low_it == up_it) { continue; }

          // create theoretical spectrum
          PeakSpectrum theo_spectrum;

          // add peaks for b and y ions with charge 1
          spectrum_generator.getSpectrum(theo_spectrum, candidate, 1, 1);

          // sort by mz
          theo_spectrum.sortByPosition();

          for (; low_it != up_it; ++low_it)
          {
            const Size& scan_index = low_it->second;
            const PeakSpectrum& exp_spectrum = spectra[scan_index];
            // const int& charge = exp_spectrum.getPrecursors()[0].getCharge();
            //const double& score = HyperScore::compute(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);
            const double& score = -log(calculatePValue_(candidate, rems[scan_index], cums[scan_index]));

            //if (score == 0) { continue; } // no hit?

            // add peptide hit
            AnnotatedHit_ ah;
            ah.sequence = c;
            ah.peptide_mod_index = mod_pep_idx;
            ah.score = score;

#ifdef _OPENMP
            omp_set_lock(&(annotated_hits_lock[scan_index]));
            {
#endif
              annotated_hits[scan_index].push_back(ah);

              // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
              if (annotated_hits[scan_index].size() >= 2 * report_top_hits_)
              {
                std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits_, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
                annotated_hits[scan_index].resize(report_top_hits_); 
              }
#ifdef _OPENMP
            }
            omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
          }
        }
      }
    }
    endProgress();

    OPENMS_LOG_INFO << "Proteins: " << count_proteins << endl;
    OPENMS_LOG_INFO << "Peptides: " << count_peptides << endl;
    OPENMS_LOG_INFO << "Processed peptides: " << processed_petides.size() << endl;

    startProgress(0, 1, "Post-processing PSMs...");
    SimpleSearchEngineAlgorithm::postProcessHits_(spectra, 
      annotated_hits, 
      protein_ids, 
      peptide_ids, 
      report_top_hits_,
      fixed_modifications, 
      variable_modifications, 
      modifications_max_variable_mods_per_peptide_,
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      precursor_mass_tolerance_,
      fragment_mass_tolerance_,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      in_db
      );
    endProgress();

    // add meta data on spectra file
    StringList ms_runs;
    spectra.getPrimaryMSRunPath(ms_runs);
    protein_ids[0].setPrimaryMSRunPath(ms_runs);

    // reindex peptides to proteins
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string", "DECOY_");
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", enzyme_);
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
      {
        return ExitCodes::INPUT_FILE_EMPTY;       
      }
      else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
      {
        return ExitCodes::UNEXPECTED_RESULT;
      }
      else
      {
        return ExitCodes::UNKNOWN_ERROR;
      }
    } 

#ifdef _OPENMP
    // free locks
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_destroy_lock(&(annotated_hits_lock[i])); }
#endif

    return ExitCodes::EXECUTION_OK;
  }

} // namespace

