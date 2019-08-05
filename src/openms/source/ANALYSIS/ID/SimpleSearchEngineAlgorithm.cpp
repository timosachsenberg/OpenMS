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


//#define RES_EV_DEBUG 1

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
  void SimpleSearchEngineAlgorithm::preprocessResidueEvidence_(
    PeakMap& exp,
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm,
    const std::set<const Residue*>& aas,
    vector<ResidueEvidenceMatrix>& rems, 
    vector<CumScoreHistogram>& cums)
  {
    rems.clear();
    cums.clear();
    rems = vector<ResidueEvidenceMatrix>(exp.size(), ResidueEvidenceMatrix());
    cums = vector<CumScoreHistogram>(exp.size(), CumScoreHistogram());

    // sort
    exp.sortSpectra(true);

    // remove 0 or negative intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    // deisotope (not part of original implementation)
    for (auto & spectrum : exp)
    {
      Deisotoper::deisotopeAndSingleCharge(spectrum, 
        fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
        1, 3,   // min / max charge 
        false,  // keep all peaks
        2, 5,  // min / max isopeaks 
        false, false);  // don't convert fragment m/z to mono charge and don't annotate charge
        
      // remove precursor peaks?  disabled by default in tide but in publication on residue evidence score
      const double precurMz = spectrum.getPrecursors()[0].getMZ();
      auto left = spectrum.MZBegin(precurMz - 0.02);
      auto right = spectrum.MZEnd(precurMz + 0.02); 
      std::transform(left, right, left, [](Peak1D& p) { p.setIntensity(0); return p; }); // set to zero - removed in threshold mower
    }

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
        auto right = spec.MZEnd(min_mz + (i + 1) * segment_mz + 0.001); // TODO: check if 0.001 required

        auto max_intensity_it = std::max_element(left, right, 
          [](const Peak1D& a, const Peak1D& b) 
          { 
            return a.getIntensity() < b.getIntensity(); 
          });
        double max_intensity = max_intensity_it->getIntensity();
        std::transform(left, right, left, [max_intensity](Peak1D& p) { p.setIntensity(p.getIntensity() / max_intensity * 50.0); return p; });        
      }    

      // calculate ranks and assign intensity/ranks as new intensities. 
      // Identical intensities (e.g, the ten 50.0 the segment maxima) get the same rank - seems to be critical for score perfomance
      vector<double> intensities(spec.size());
      for (Size i = 0; i != spec.size(); ++i) { intensities[i] = spec[i].getIntensity(); }
      std::sort(intensities.begin(), intensities.end(), greater<double>()); // sort decreasing
      OPENMS_LOG_DEBUG << "before unique: " <<  intensities.size() << endl;
      intensities.resize(std::distance(intensities.begin(), std::unique(intensities.begin(), intensities.end(), 
          [](double l, double r) { return std::abs(l - r) < 1e-6; })
         ));
      OPENMS_LOG_DEBUG << "after unique: " <<  intensities.size() << endl;
      // calculate 1.0 / rank intensity
      for (Size i = 0; i != spec.size(); ++i) 
      {
        int rank = 1 + lower_bound(intensities.begin(), intensities.end(), spec[i].getIntensity() + 1e-3, greater<double>()) - intensities.begin();
        spec[i].setIntensity(1.0 / (double)rank);
      }
       
    } 

    // MzMLFile().store("testX.mzML", exp);

    // TODO: depends on terminal mods
    double cTermMass = Residue::getInternalToCTerm().getMonoWeight();
    double nTermMass = Residue::getInternalToNTerm().getMonoWeight();
    OPENMS_LOG_DEBUG << "nTermMass: " << nTermMass << endl;
    OPENMS_LOG_DEBUG << "cTermMass: " << cTermMass << endl;

    // determine which bin each amino acid mass is in
    vector<int> aaMassBin;
    vector<double> aaMass;
    for (const auto& r : aas)
    {
      const float mass = r->getMonoWeight(Residue::Internal); // TODO: check if internal is correct
      aaMass.push_back(mass);
      int binMass = mass2bin_(mass);  
      aaMassBin.push_back(binMass);
    }

    #pragma omp parallel for
    for (SignedSize spectrum_index = 0; spectrum_index < exp.size(); ++spectrum_index)
    {
      MSSpectrum& spectrum = exp[spectrum_index];
      const double precurMz = spectrum.getPrecursors()[0].getMZ(); 
      const int precurCharge = spectrum.getPrecursors()[0].getCharge();
      const double precursorMass = precurMz * precurCharge - precurCharge * Constants::PROTON_MASS_U;

      // calculate residue evidence matrix rem
      // rem[row][column]: rows are amino acids (or modified amino acids), columns are mass bins 
      ResidueEvidenceMatrix rem = ResidueEvidenceMatrix(aas.size(), vector<double>(maxPrecurMassBin, 0));  
//TODO: changed TOOOOOOOOOOOOOOOOOOOOOOOOOODOOOOOOOOOOOODDDDOOOOOOO but should be correct
      createResidueEvidenceMatrix(spectrum, aas, 0.02,  mass2bin_(precursorMass), rem); // 0.02 = default fragment mass tolerance for high-res
      
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

        double q = ceil(precursorMass / 57.02147); // lightest amino acid
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
            double left_sum(0);
            for (size_t a = 0; a != aaMassBin.size(); ++a)
            {
              size_t a_bin = aaMassBin[a];
              size_t R_am = rem[a][m];  // residue evidence that a peptide(-prefix) with mass m (z=0) ended with a
              if ((int)s - (int)R_am < 0 || (int)m - (int)a_bin < 0) continue;
             
              const size_t left_bin = m - a_bin;
              const size_t aa_mass = aaMass[a] * 1000.0; // int scaling
              double pAA(1);
              if (left_bin == mass2bin_(nTermMass)) pAA = dAAFreqN.at(aa_mass);
              else if (m == mass2bin_(precursorMass)) pAA = dAAFreqC.at(aa_mass);
              else pAA = dAAFreqI.at(aa_mass);

              left_sum += C[s - R_am][m - a_bin] * pAA;
            }
            C[s][m] += left_sum; // += here because we don't want to overwrite the only non-zero entry C[0][mass2bin_(nTermMass)] right from the start
          }
        }
      
        // calculate (cumulative) sum of counts for precursor mass with score larger than s once
        size_t m = mass2bin_(precursorMass);
        OPENMS_LOG_DEBUG << "Precursor mass (neutral) bin: " << m << endl;
        double cumsum(0);
        for (int s = s_max; s >= 0; --s) 
        { 
          cumsum += C[s][m]; 
          cumsumsC_sm[s] = cumsum;
        }

#ifdef RES_EV_DEBUG
	int r_index(0);
        for (auto & r : C)
        {
          cout << "Counts for score: " << r_index << endl;
          String line;
          bool only_zeros = true;
          int m_bin(0);
          for (auto & c : r)
          {
            if (c != 0) only_zeros = false;
            if (m_bin == m)
            {
              line +="(" + String((int)log10(c+1)) + ") ";
              break; // no entries after this one relevant 
            }
            else 
              line +=String((int)log10(c+1)) + " ";
            ++m_bin;
          }
          if (!only_zeros) cout << line << endl; else cout << "zeros only" << endl;
          ++r_index;
        }
#endif
        // output cumulative count distribution
        /// for (int s = 0; s <= s_max; ++s)  { if (s != 0) cout << "Score: " << s << "\t" << cumsumsC_sm[s] << endl; }

        // cumsumsC_sm[0] holds the cumulative sum for all scores ([0] means total count for all scores > 0) - (denominator in eq. 10)
        OPENMS_LOG_DEBUG << "Sum of C_{s,m}=" << cumsumsC_sm[0] << std::endl;

        cums[spectrum_index].swap(cumsumsC_sm);
        rems[spectrum_index].swap(rem);
      }
    }    
  }

  size_t SimpleSearchEngineAlgorithm::calculateRawResEv_(const AASequence& candidate, const ResidueEvidenceMatrix& rem) const
  {
    // calculate residue evidence score
    size_t res_ev_score(0);
    double m(Residue::getInternalToBIon().getMonoWeight());

    // TODO: add modified terminus
    for (Size i = 0; i != candidate.size(); ++i)
    {
      const double aa_mass = candidate[i].getMonoWeight(Residue::Internal);
      m += aa_mass;
      size_t m_bin = mass2bin_(m);
      
      auto it = mapResidue2Index_.find(&(candidate[i])); // TODO: handle modified residues / I vs. L
      if (it != mapResidue2Index_.end())
      {
#ifdef RES_EV_DEBUG
        cout << "candidate residue: '" << candidate[i].getOneLetterCode() << "' mass_bin: " << m_bin << " index: " << it->second << " ev:" << rem[it->second][m_bin] << " sum:" << res_ev_score << endl; 
#endif
        res_ev_score += rem[it->second][m_bin];  
      }
      else
      {
         OPENMS_LOG_DEBUG << "Warning: residue '" << candidate[i].getOneLetterCode() << "' not found." << endl;
      }
    }
    // OPENMS_LOG_DEBUG << "ResEv: Raw residue evidence score of " << candidate.toString() << ": " << res_ev_score << std::endl;
    return res_ev_score;
  }

  // probability of other peptide of same mass with higher score
  double SimpleSearchEngineAlgorithm::calculatePValue_(
    const double res_ev_score, 
    const vector<double>& cumsumsC_sm) const
  {
    if (res_ev_score <= 0) return 1.0;

    OPENMS_PRECONDITION(cumsumsC_sm[0] > 0, "No DP-amino acid sequence exist matching to mass bin of precursor.")
    OPENMS_PRECONDITION(cumsumsC_sm[res_ev_score] > 0, "At least one peptide needs to match score.")

    if (cumsumsC_sm[res_ev_score] == 0) return 1e-200; // TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO fix

    // number of peptides with score equal or better then residue evidence score of candidate / total peptide count
    //cout << "(res-ev: p count>=res-ev count>=0) " << res_ev_score << ": " << cumsumsC_sm[res_ev_score]  / cumsumsC_sm[0] << " " << cumsumsC_sm[res_ev_score] << " " << cumsumsC_sm[0] << std::endl;
    const double p_value = cumsumsC_sm[res_ev_score]  / cumsumsC_sm[0]; // TODO: check why +1 is needed here? should be at least 1 but there are zeros?
    return p_value;
  }

  unsigned int SimpleSearchEngineAlgorithm::mass2bin_(double mass, int charge /*=1*/) 
  {
    const double bin_width = 1.0005079;
    const double bin_offset = 0.4; // 0.4 tide & comet default - paper default: 0.68
    return (unsigned int)((mass + (charge - 1) * Constants::PROTON_MASS_U) / (charge * bin_width) + 1.0 - bin_offset);
  }

  double SimpleSearchEngineAlgorithm::bin2mass_(int bin, int charge /*=1*/) 
  {
    const double bin_width = 1.0005079;
    const double bin_offset = 0.4; // tide & comet default - paper default: 0.68
    return (bin - 1.0 + bin_offset) * charge * bin_width + (charge - 1) * Constants::PROTON_MASS_U;
  }


// adapted from Crux/Tide Andy Lin
void SimpleSearchEngineAlgorithm::createResidueEvidenceMatrix(
    const MSSpectrum& spectrum,
    const std::set<const Residue*>& aas,
    double fragment_tolerance_Da,
    size_t max_precursor_mass_bin,
    vector<vector<double> >& residueEvidenceMatrix) 
  {
    const double precurMz = spectrum.getPrecursors()[0].getMZ(); 
    const int precurCharge = spectrum.getPrecursors()[0].getCharge();
    const double precursorMass = precurMz * precurCharge - precurCharge * Constants::PROTON_MASS_U; // neutral mass
    const int granularityScale = 25;

    // TODO: these are different for fixed modified termini
    double cTermMass = Residue::getInternalToCTerm().getMonoWeight();
    double nTermMass = Residue::getInternalToNTerm().getMonoWeight();

    // determine which bin each amino acid mass is in
    vector<int> aaMassBin;
    vector<double> aaMass;
    for (const auto& r : aas)
    {
      const float mass = r->getMonoWeight(Residue::Internal); 
      aaMass.push_back(mass);
      int binMass = mass2bin_(mass);  
      aaMassBin.push_back(binMass);
#ifdef RES_EV_DEBUG 
      OPENMS_LOG_DEBUG << "ResEv: Amino acid: " << r->getName() << "(" << r->getOneLetterCode() << "):"  << binMass << "\t" << endl;
#endif
    }

  ////////////////////////////////////////////////////////////////////////////
  // Populate ionMass and ionMassBin for b ions in 1+ charge state
  vector<double> ionMass;
  vector<int> ionMassBin;
  vector<double> ionIntens;

  // add nTerm peak (no proton charge)
  ionMass.push_back(nTermMass);
  ionMassBin.push_back(mass2bin_(nTermMass));
  ionIntens.push_back(0.0);
  
#ifdef RES_EV_DEBUG
  cout << "nTerm: " << mass2bin_(nTermMass) << "\t" << nTermMass << endl;
#endif

  // 2. add residues
  for (size_t ion = 0; ion < spectrum.size(); ++ion) 
  {
    // add peak m/z (contains proton charge - assuming charge +1)
    double tmpIonMass = spectrum[ion].getMZ();
    int binTmpIonMass = mass2bin_(tmpIonMass);

//    cout << "b+: " << tmpIonMass << endl;

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
    ionIntens.push_back(spectrum[ion].getIntensity());
  }
  // 3. add peak for last amino acid without cTerm (no proton charge)
  ionMass.push_back(precursorMass - cTermMass); // TODO: should this get a proton charge? otherwise distance calculation might not be correct
  ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

#ifdef RES_EV_DEBUG
  cout << "b+: ";
  for (Size i = 0; i != ionMassBin.size(); ++i) cout << " " << ionMassBin[i];
  cout << endl;
#endif

  // add evidence fo
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
    double tmpIonMass = precursorMass - spectrum[ion].getMZ() + (2.0 * Constants::PROTON_MASS_U); //TODO: should this be H instead of H+
    int binTmpIonMass = mass2bin_(tmpIonMass);

    
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

  // why reverse? to find first all mass bins matching a mass?
  reverse(ionMass.begin(), ionMass.end());
  reverse(ionMassBin.begin(), ionMassBin.end());
  reverse(ionIntens.begin(), ionIntens.end());

#ifdef RES_EV_DEBUG
  cout << "y+: ";
  for (Size i = 0; i != ionMassBin.size(); ++i) cout << " " << ionMassBin[i];
  cout << endl;
#endif

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
      int binTmpIonMass = mass2bin_(tmpIonMass);
    
      ionMass.push_back(tmpIonMass);
      ionMassBin.push_back(binTmpIonMass);
      ionIntens.push_back(spectrum[ion].getIntensity());
    }
    ionMass.push_back(precursorMass - cTermMass);
    ionMassBin.push_back(mass2bin_(precursorMass - cTermMass));
    ionIntens.push_back(0.0);

#ifdef RES_EV_DEBUG
  cout << "b++: ";
  for (Size i = 0; i != ionMassBin.size(); ++i) cout << " " << ionMassBin[i];
  cout << endl;
#endif

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
      int binTmpIonMass = mass2bin_(tmpIonMass);

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

#ifdef RES_EV_DEBUG
  cout << "y++: ";
  for (Size i = 0; i != ionMassBin.size(); ++i) cout << " " << ionMassBin[i];
  cout << endl;
#endif
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
    for (size_t aa_index = 0; aa_index < aas.size(); ++aa_index) 
    {
      if (residueEvidenceMatrix[aa_index][i] > maxEvidence) { maxEvidence = residueEvidenceMatrix[aa_index][i]; }
    }
  }

#ifdef RES_EV_DEBUG
  cout << "max evidence before scaling: " << maxEvidence << endl;
#endif

  // Discretize residue evidence so largest value is residueEvidenceIntScale
  double residueEvidenceIntScale = (double)granularityScale;
  size_t truncated_evidence(0);
  for (int i = 0; i < max_precursor_mass_bin; i++) 
  {
    for (size_t aa_index = 0; aa_index < aas.size(); ++aa_index) 
    {
      double residueEvidence = residueEvidenceMatrix[aa_index][i];
      if (residueEvidence > 0) 
      {
        residueEvidenceMatrix[aa_index][i] = round(residueEvidenceIntScale * residueEvidence / maxEvidence);  
        //if (residueEvidenceMatrix[aa_index][i] == 0) { ++truncated_evidence; residueEvidenceMatrix[aa_index][i] = 1; } // prevent truncation TODOOOOOOOOOOOOOOOO makes this sense
      }
    }
  }
  //cout << "Truncated evidence: " << truncated_evidence << endl; // these evidences will not contribute anymore (Improvements?)

#ifdef RES_EV_DEBUG
  cout << "max evidence after scaling: " << residueEvidenceIntScale << endl;
#endif

/*
  size_t aa_index(0);
  for (auto & r : residueEvidenceMatrix)
  {
    String line;
    bool only_zeros = true;
    size_t index(0);
    for (auto & c : r)
    {
      if (c > 0) only_zeros = false;
      line += String(index)+":"+String(c) + " ";
      index++;
    }
    auto it = aas.begin();
    std::advance(it, aa_index);
    cout << (*it)->getName() << endl;
    if (!only_zeros) cout << line << endl;
    ++aa_index;
  }
*/
}

// adapted from code written by Andy Lin in Feb 2018
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
      // bIonMass is named correctly because we assume that all y ions have been converted to their b+ ion equivalent
      // before ionMass is passed in
      double bIonMass1 = ionMass[ion];
      int bIonMassBin1 = ionMassBin[ion];
      double bIonIntensity1 =  ionIntens[ion];
   
      // for all amino acids:
      for (size_t aa_index = 0; aa_index < nAA; ++aa_index) 
      { 
        // calculate bin with amino acid distance to current peak
        int newResMassBin = bIonMassBin1 + aaMassBin[aa_index];
        
        // Find all ion mass bins that match newResMassBin
        int index = find(ionMassBin.begin(), ionMassBin.end(), newResMassBin) - ionMassBin.begin();
        double score = 0.0;
        for (int i = index; i < ionMassBin.size(); ++i) // there can be several spectrum peaks in the newResMassBin, process all
        {
          if (newResMassBin != ionMassBin[i]) { break; } 
          const double bIonMass2 = ionMass[i];
          const double ionMassDiff =  bIonMass2 - bIonMass1;
          const double bIonIntensity2 = ionIntens[i];

          const double aa_mass = aaMass[aa_index];
          // linearly weight deviation from expected amino acid distance
          const double aaTolScore = 1.0 - (std::fabs(ionMassDiff - aa_mass) / fragment_tolerance_Da);

          if (aaTolScore > 0.0) // in tolerance window?
          {         
            const double tmpScore = aaTolScore * (bIonIntensity1 + bIonIntensity2);  // weighted rank score
            
            if (tmpScore > score) { score = tmpScore; }
          }
        }

        // Add evidence to matrix
        // Use -1 since all mass bins are index 1 instead of index 0
        // Bounds checks. Only add score if smaller than precursor mass bin.
        // When assuming each fragment peak is a 2+ charge, it is possible
        // to have a fragment peak larger than precursor mass (ie why
        // bounds check is needed).
        if (score != 0.0 && newResMassBin <= maxPrecurMassBin) 
        {
          residueEvidenceMatrix[aa_index][newResMassBin - 1] += score; // add evidence for neutral mass (hence -1 proton)
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

  void SimpleSearchEngineAlgorithm::getAminoAcidFrequencies_(
      const vector<FASTAFile::FASTAEntry>& fasta_db, 
      const ProteaseDigestion& digestor,   
      map<double, double>& dAAFreqN,
      map<double, double>& dAAFreqI,
      map<double, double>& dAAFreqC,
      map<double, double>& dAAMass 
    )
  {
    map<double, double> cN;
    map<double, double> cI;
    map<double, double> cC;
    map<double, double> cAll;
    size_t cN_total(0);
    size_t cI_total(0);
    size_t cC_total(0);

    for (int fasta_index = 0; fasta_index != fasta_db.size(); ++fasta_index)
    {
      vector<StringView> current_digest;
      digestor.digestUnmodified(fasta_db[fasta_index].sequence, current_digest, peptide_min_size_, peptide_max_size_);
      for (auto const & c : current_digest)
      { 
        const String current_peptide = c.getString();
        if (current_peptide.find_first_of("XBZ") != std::string::npos) { continue; }
        const AASequence aas = AASequence::fromString(current_peptide);
        // count how often a (unique) AA mass occurs at N-,C-term or internal 
        for (Size i = 0; i != aas.size(); ++i)
        {
          // TODO: mods
          size_t m = aas[i].getMonoWeight(Residue::Internal) * 1000.0; // int scale
          if ( i == 0 ) 
          { 
            cN[m] += 1; ++cN_total; 
          }
          else if ( i > 0 && i < aas.size()-1) 
          { 
            cI[m] += 1; ++cI_total; 
          }
          else if ( i == aas.size()-1 ) { cC[m] += 1; ++cC_total;}
          cAll[m] += 1;
        } 
      }
    }

    // add one pseudo count so we fill all amino acids
    for (auto & m : cAll)
    {
      ++cN[m.first]; ++cN_total; 
      ++cI[m.first]; ++cI_total; 
      ++cC[m.first]; ++cC_total;
    } 

    for (const auto& m : cN) dAAFreqN[m.first] = m.second / (double)cN_total; 
    for (const auto& m : cI) dAAFreqI[m.first] = m.second / (double)cI_total; 
    for (const auto& m : cC) dAAFreqC[m.first] = m.second / (double)cC_total; 
    cout << "AAFrequencies: N,C,internal:" << endl;
    for (const auto& m : dAAFreqN) cout << m.first << " N-term " << m.second << endl; 
    for (const auto& m : dAAFreqI) cout << m.first << " Intern " << m.second << endl; 
    for (const auto& m : dAAFreqC) cout << m.first << " C-term " << m.second << endl; 
  }

  SimpleSearchEngineAlgorithm::ExitCodes SimpleSearchEngineAlgorithm::search(const String& in_mzML, const String& in_db, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids) 
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

/* 
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
*/

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

    getAminoAcidFrequencies_(fasta_db, 
      digestor,   
      dAAFreqN,
      dAAFreqI,
      dAAFreqC,
      dAAMass 
    );

#ifdef RES_EV_DEBUG
    TheoreticalSpectrumGenerator test_spectrum_generator;
    Param test_param(test_spectrum_generator.getParameters());
    test_param.setValue("add_first_prefix_ion", "true");
    test_param.setValue("add_b_ions", "false"); 
    test_param.setValue("add_y_ions", "true");
    test_spectrum_generator.setParameters(test_param);
    fasta_db[0].sequence = "SSSSSSSSSS";
//    fasta_db[1].sequence = "AAAAAAAAAA";
//    fasta_db[2].sequence = "GGGGGGGGGG";
    fasta_db.resize(1);
    for (Size i = 0; i != fasta_db.size(); ++i)
    {
      spectra[i].clear(true);
      test_spectrum_generator.getSpectrum(spectra[i], AASequence::fromString(fasta_db[i].sequence), 1, 1);
      for (const Peak1D & p : spectra[i]) { cout << p.getMZ() << "\t" << p.getIntensity() << endl; }
      vector<Precursor> pcs;
      pcs.push_back(Precursor());
      pcs[0].setCharge(2);
      pcs[0].setMZ(AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, pcs[0].getCharge())/pcs[0].getCharge());
      cout << fasta_db[i].sequence << "\t" << pcs[0].getMZ() << "\t neutral mass:" << AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, 0)  
           << "\t z1: " << AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, 1)  << endl; 
      spectra[i].setPrecursors(pcs);
    }

    test_param.setValue("add_b_ions", "true"); 
    test_param.setValue("add_y_ions", "false");
    test_spectrum_generator.setParameters(test_param);
    for (Size i = 0; i != fasta_db.size(); ++i)
    {
      spectra[fasta_db.size() + i].clear(true);
      test_spectrum_generator.getSpectrum(spectra[fasta_db.size() + i], AASequence::fromString(fasta_db[i].sequence), 1, 1);
      for (const Peak1D & p : spectra[fasta_db.size() + i]) { cout << p.getMZ() << "\t" << p.getIntensity() << endl; }
      vector<Precursor> pcs;
      pcs.push_back(Precursor());
      pcs[0].setCharge(2);
      pcs[0].setMZ(AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, pcs[0].getCharge())/pcs[0].getCharge());
      cout << fasta_db[i].sequence << "\t" << pcs[0].getMZ() << "\t neutral mass:" << AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, 0)  
           << "\t z1: " << AASequence::fromString(fasta_db[i].sequence).getMonoWeight(Residue::Full, 1)  << endl; 
      spectra[fasta_db.size() + i].setPrecursors(pcs);
    }
#endif


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

    // residue evidence matrix and cumulative peptide vs score counts
    startProgress(0, 1, "Filtering spectra...");
    std::vector<ResidueEvidenceMatrix> rems;
    std::vector<CumScoreHistogram> cums;

    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural20");
    size_t index(0);
    for (auto & r : aas)
    {
      mapResidue2Index_[r] = index;
      ++index;
    }
    preprocessResidueEvidence_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, aas, rems, cums);

#ifdef RES_EV_DEBUG
     double test_score = -log(calculatePValue_(AASequence::fromString(fasta_db[0].sequence), rems[0], cums[0]));
     cout << "test_score: " << test_score << endl;
     test_score = -log(calculatePValue_(AASequence::fromString(fasta_db[0].sequence), rems[1], cums[1]));
     cout << "test_score: " << test_score << endl;
#endif

    endProgress();

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

          // cout << "ResEv: Scoring: " << candidate.toString() << endl;

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


          for (; low_it != up_it; ++low_it)
          {
            const Size& scan_index = low_it->second;
            // cout << "ResEv: Matching scan_index: " << scan_index << endl;
            // const PeakSpectrum& exp_spectrum = spectra[scan_index];

            double score;
            score = calculateRawResEv_(candidate, rems[scan_index]);
#ifdef RES_EV_DEBUG
            cout << candidate.toString() << " res-ev: " << score << endl;
#endif
            score = -log10(calculatePValue_(score, cums[scan_index]));
#ifdef RES_EV_DEBUG
            cout << candidate.toString() << " p-value: " << score << endl;
#endif

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

