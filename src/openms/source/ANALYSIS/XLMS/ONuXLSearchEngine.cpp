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

#include <OpenMS/ANALYSIS/XLMS/ONuXLSearchEngine.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLSpectrumProcessing.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <map>

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

using namespace std;
using namespace OpenMS;
using namespace OpenNuXL;

 ONuXLSearchEngine::ONuXLSearchEngine()
    : DefaultParamHandler("ONuXLSearchEngine")
  {
    defaults_.setValue("precursor:mass_tolerance", 10.0, "Precursor mass tolerance (+/- around precursor m/z)");
    defaults_.setSectionDescription("precursor", "Precursor (Parent Ion) Options");

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered.");

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    defaults_.setValue("precursor:isotopes", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)");

    defaults_.setValue("fragment:mass_tolerance", 10.0, "Fragment mass tolerance (+/- around fragment m/z)");

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment mass");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);
    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");


    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("modifications:fixed", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", all_mods);
    defaults_.setValue("modifications:variable", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", all_mods);
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    defaults_.setValue("peptide:min_size", 6, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 1000, "Maximum size a peptide may have after digestion to be considered in the search.");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");

    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    defaults_.setValue("peptide:enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("peptide:enzyme", all_enzymes);
    defaults_.setSectionDescription("peptide", "Peptide Options");

    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setSectionDescription("report", "Reporting Options");

    // RNPxl specific
    defaults_.setValue("RNPxl:length", 2, "Oligonucleotide maximum length. 0 = disable search for RNA variants.");

    defaults_.setValue("RNPxl:sequence", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)");

    defaults_.setValue("RNPxl:target_nucleotides", 
                        StringList{"A=C10H14N5O7P", "C=C9H14N3O8P", "G=C10H14N5O8P", "U=C9H13N2O9P"}, 
                        "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG"
                        );

    defaults_.setValue("RNPxl:mapping", StringList{"A->A", "C->C", "G->G", "U->U"}, "format: source->target e.g. A->A, ..., U->U, U->X");

    // define if nucleotide can cross-link (produce y,b,a,immonium-ion shifts) in addition to marker ions
    defaults_.setValue("RNPxl:can_cross_link",                        
                        "U", 
                        "format: 'U' if only U forms cross-links. 'CATG' if C, A, G, and T form cross-links.");

    StringList modifications;
    modifications.emplace_back("U:");
    modifications.emplace_back("U:-H2O");
    modifications.emplace_back("U:-H2O-HPO3");
    modifications.emplace_back("U:-HPO3");

    // fragment adducts that may occur for every precursor adduct (if chemically feasible in terms of elements may not be negative)
    StringList fragment_adducts = {"U:C9H10N2O5;U-H3PO4", 
                                   "U:C4H4N2O2;U'", 
                                   "U:C4H2N2O1;U'-H2O",
                                   "U:C3O;C3O",
                                   "U:C9H13N2O9P1;U",
                                   "U:C9H11N2O8P1;U-H2O",
                                   "U:C9H12N2O6;U-HPO3"
                                  };    

    defaults_.setValue("RNPxl:fragment_adducts",                         
                        fragment_adducts, 
                        "format: [target nucleotide]:[formula] or [precursor adduct]->[fragment adduct formula];[name]: e.g., 'U:C9H10N2O5;U-H3PO4' or 'U:U-H2O->C9H11N2O8P1;U-H2O'," 
                        );

    defaults_.setValue("RNPxl:modifications", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3");

    defaults_.setValue("RNPxl:scoring", "fast", "Scoring algorithm used in prescoring (fast: total-loss, slow: all losses).");
    defaults_.setValue("RNPxl:scoring", StringList{"fast", "slow"});

    defaults_.setValue("RNPxl:decoys", "false", "Generate decoy sequences and spectra.");
    defaults_.setValidStrings("RNPxl:decoys", ListUtils::create<String>("true,false"));
    defaults_.setValue("RNPxl:CysteineAdduct", "false", "Use this flag if the +152 adduct is expected.");
    defaults_.setValidStrings("RNPxl:CysteineAdduct", ListUtils::create<String>("true,false"));
    defaults_.setValue("RNPxl:filter_fractional_mass", "false", "Use this flag to filter non-crosslinks by fractional mass.");
    defaults_.setValidStrings("RNPxl:filter_fractional_mass", ListUtils::create<String>("true,false"));
    defaults_.setValue("RNPxl:localization", "true", "Use this flag to perform crosslink localization by partial loss scoring as post-analysis.");
    defaults_.setValidStrings("RNPxl:localization", ListUtils::create<String>("true,false"));
    defaults_.setValue("RNPxl:carbon_labeled_fragments", "false", "Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).");
    defaults_.setValidStrings("RNPxl:carbon_labeled_fragments", ListUtils::create<String>("true,false"));

    defaults_.setValue("RNPxl:filter_small_peptide_mass", 600.0, "Filter precursor that can only correspond to non-crosslinks by mass.");
    defaults_.setValue("RNPxl:marker_ions_tolerance", 0.05, "Tolerance used to determine marker ions (Da).");
    defaults_.setSectionDescription("RNPxl", "RNPxl Options");

    defaultsToParam_();
  }

  void ONuXLSearchEngine::updateMembers_()
  {    
  }

  ONuXLSearchEngine::~ONuXLSearchEngine()
  {
  }

  ONuXLSearchEngine::ExitCodes ONuXLSearchEngine::run(
    const MSExperiment& in_spectra, 
    vector<FASTAFile::FASTAEntry>& fasta_db,
    const String& database_filename,
    vector<ProteinIdentification>& protein_ids,
    vector<PeptideIdentification>& peptide_ids,
    TextFile& csv_file)
  {
    if (!in_spectra.isSorted())
    {
      LOG_ERROR << "Spectra need to be sorted by m/z and RT!" << endl;
      return ExitCodes::ILLEGAL_PARAMETERS;
    }

    for (auto const & s : in_spectra)
    {
      if (s.getMSLevel() != 2)
      {
        LOG_ERROR << "Only MS2-level spectra must be provided!" << endl;
        return ExitCodes::ILLEGAL_PARAMETERS;
      }
    }

    MSExperiment spectra = in_spectra;

    // force initialization of residue db
    AASequence::fromString("ACDEFGHIJKLMNQRSTVWY");

    fast_scoring_ = param_.getValue("RNPxl:scoring") == "fast" ? true : false;

    bool generate_decoys = param_.getValue("RNPxl:decoys") == "true" ? true : false;

    Int min_precursor_charge = (int)param_.getValue("precursor:min_charge");
    Int max_precursor_charge = (int)param_.getValue("precursor:max_charge");
    double precursor_mass_tolerance = (double)param_.getValue("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = ((String)param_.getValue("precursor:mass_tolerance_unit") == "ppm");
    IntList precursor_isotopes = param_.getValue("precursor:isotopes").toIntList();

    double fragment_mass_tolerance = (double)param_.getValue("fragment:mass_tolerance");
    bool fragment_mass_tolerance_unit_ppm = ((String)param_.getValue("fragment:mass_tolerance_unit") == "ppm");

    double marker_ions_tolerance = (double)param_.getValue("RNPxl:marker_ions_tolerance");

    double small_peptide_mass_filter_threshold = (double)param_.getValue("RNPxl:filter_small_peptide_mass");

    StringList fixedModNames = param_.getValue("modifications:fixed").toStringList();
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = (int)param_.getValue("peptide:min_size");

    if (fixed_unique.size() != fixedModNames.size())
    {
      LOG_WARN << "duplicate fixed modification provided." << endl;
      return ONuXLSearchEngine::ExitCodes::ILLEGAL_PARAMETERS;
    }

    StringList varModNames = param_.getValue("modifications:variable").toStringList();
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      LOG_WARN << "duplicate variable modification provided." << endl;
      return ONuXLSearchEngine::ExitCodes::ILLEGAL_PARAMETERS;
    }

    vector<ResidueModification> fixed_modifications = ONuXLParameterParsing::getModifications(fixedModNames);
    vector<ResidueModification> variable_modifications = ONuXLParameterParsing::getModifications(varModNames);
    Size max_variable_mods_per_peptide = (int)param_.getValue("modifications:variable_max_per_peptide");

    size_t report_top_hits = (size_t)((int)param_.getValue("report:top_hits"));

    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides = param_.getValue("RNPxl:target_nucleotides").toStringList();

    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings = param_.getValue("RNPxl:mapping").toStringList();

    // read list of nucleotides that can directly cross-link
    // these are responsible for shifted fragment ions. Their fragment adducts thus determine which shifts will be observed on b-,a-,y-ions
    String can_cross_link = (String)param_.getValue("RNPxl:can_cross_link");
    for (auto c : can_cross_link) { can_xl_.insert(c); }

    StringList modifications = param_.getValue("RNPxl:modifications").toStringList();

    String sequence_restriction = (String)param_.getValue("RNPxl:sequence");

    Int max_nucleotide_length = (int)param_.getValue("RNPxl:length");

    bool cysteine_adduct = (String)param_.getValue("RNPxl:CysteineAdduct") == "true" ? true : false;

    localization_ = (String)param_.getValue("RNPxl:localization") == "true" ? true : false;

    // generate all precursor adducts
    ONuXLModificationMassesResult mm;
    if (max_nucleotide_length != 0)
    {
      mm = ONuXLModificationsGenerator::initModificationMassesRNA(
            target_nucleotides,
            can_xl_,
            mappings,
            modifications, 
            sequence_restriction, 
            cysteine_adduct, 
            max_nucleotide_length);
    }

    mm.mod_masses[""] = 0; // insert "null" modification otherwise peptides without RNA will not be searched
    mm.mod_combinations[""].insert("none");

    // parse tool parameter and generate all fragment adducts

    // first, we determine which fragments adducts can be generated from a single nucleotide (that has no losses)
    ONuXLParameterParsing::NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts = 
        ONuXLParameterParsing::getTargetNucleotideToFragmentAdducts(param_.getValue("RNPxl:fragment_adducts").toStringList());

    // calculate all feasible fragment adducts from all possible precursor adducts
    ONuXLParameterParsing::PrecursorsToMS2Adducts all_feasible_fragment_adducts = 
        ONuXLParameterParsing::getAllFeasibleFragmentAdducts(mm, nucleotide_to_fragment_adducts, can_xl_);

    startProgress(0, 1, "Filtering spectra...");
    const bool convert_to_single_charge = false;  // whether to convert fragment peaks with isotopic patterns to single charge
    const bool annotate_charge = false;  // whether the charge and type is annotated
    ONuXLSpectrumProcessing::preprocess(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, convert_to_single_charge, annotate_charge);
    endProgress();

    // build multimap of precursor mass to scan index (and perform some mass and length based filtering)
    using MassToScanMultiMap = multimap<double, pair<Size, int>>;
    MassToScanMultiMap multimap_mass_2_scan_index;  // map precursor mass to scan index and (potential) isotopic missassignment
    mapPrecursorMassesToScans_(min_precursor_charge,
                              max_precursor_charge,
                              precursor_isotopes,
                              small_peptide_mass_filter_threshold,
                              peptide_min_size,
                              spectra,
                              multimap_mass_2_scan_index);

    // initialize spectrum generators (generated ions, etc.)
    TheoreticalSpectrumGenerator total_loss_spectrum_generator;
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    TheoreticalSpectrumGenerator a_ion_sub_score_spectrum_generator;
    TheoreticalSpectrumGenerator immonium_ion_sub_score_spectrum_generator;
    TheoreticalSpectrumGenerator precursor_ion_sub_score_spectrum_generator;
    initializeSpectrumGenerators(total_loss_spectrum_generator,
                                 partial_loss_spectrum_generator,
                                 a_ion_sub_score_spectrum_generator,
                                 immonium_ion_sub_score_spectrum_generator,
                                 precursor_ion_sub_score_spectrum_generator);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(2 * report_top_hits); }

#ifdef _OPENMP     
    // we want to do locking at the spectrum level so we get good parallelisation 
    vector<omp_lock_t> annotated_hits_lock(annotated_hits.size());
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_init_lock(&(annotated_hits_lock[i])); }
#endif

    // generate decoy protein sequences by reversing them
    if (generate_decoys)
    {
      startProgress(0, 1, "Generate decoys...");

      // append decoy proteins
      const size_t old_size = fasta_db.size();
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = fasta_db[i];
        e.sequence.reverse();
        e.identifier = "DECOY_" + e.identifier;
        fasta_db.emplace_back(e);
      }
      endProgress();
    }

    const Size missed_cleavages = (int)param_.getValue("peptide:missed_cleavages");
    ProteaseDigestion digestor;
    digestor.setEnzyme((String)param_.getValue("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

    startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    set<StringView> processed_petides;

    // set minimum size of peptide after digestion
    Size min_peptide_length = (Size)((int)param_.getValue("peptide:min_size"));
    Size max_peptide_length = (Size)((int)param_.getValue("peptide:max_size"));

    Size count_proteins(0), count_peptides(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 100)
#endif
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif

      ++count_proteins;

      IF_MASTERTHREAD
      {
        setProgress((SignedSize)count_proteins);
      }

      vector<StringView> current_digest;

      auto const & current_fasta_entry = fasta_db[fasta_index];

      digestor.digestUnmodified(current_fasta_entry.sequence, current_digest, min_peptide_length, max_peptide_length);

      for (auto cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          // skip peptide (and all modified variants) if already processed
          if (processed_petides.find(*cit) != processed_petides.end())
          {
            already_processed = true;
          }
        }

        if (already_processed) { continue; }

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          processed_petides.insert(*cit);
        }

#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;
        vector<AASequence> all_modified_peptides;

        const String unmodified_sequence = cit->getString();

#ifdef _OPENMP
#pragma omp critical (residuedb_access)
#endif
        {
           // only process peptides without ambiguous amino acids (placeholder / any amino acid)
          if (unmodified_sequence.find_first_of("XBZ") == std::string::npos)
          {
            AASequence aas = AASequence::fromString(unmodified_sequence);
            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
          }
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass_without_RNA = fixed_and_variable_modified_peptide.getMonoWeight();

          //create empty theoretical spectrum.  total_loss_spectrum_z2 contains both charge 1 and charge 2 peaks
          PeakSpectrum total_loss_spectrum_z1, total_loss_spectrum_z2;

          // spectrum containing additional peaks for sub scoring
          PeakSpectrum immonium_sub_score_spectrum, 
                       a_ion_sub_score_spectrum, 
                       precursor_sub_score_spectrum,
                       marker_ions_sub_score_spectrum;

          // iterate over all RNA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          Size rna_mod_index = 0;

          // TODO: track the XL-able nt here
          for (std::map<String, double>::const_iterator rna_mod_it = mm.mod_masses.begin(); rna_mod_it != mm.mod_masses.end(); ++rna_mod_it, ++rna_mod_index)
          {            
            const double precursor_rna_weight = rna_mod_it->second;
            const double current_peptide_mass = current_peptide_mass_without_RNA + precursor_rna_weight; // add RNA mass
            // TODO: const char xl_nucleotide; // can be none

            // determine MS2 precursors that match to the current peptide mass
            MassToScanMultiMap::const_iterator low_it, up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + current_peptide_mass * precursor_mass_tolerance * 1e-6);
            }
            else // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + precursor_mass_tolerance);
            }

            if (low_it == up_it) { continue; } // no matching precursor in data

            // add peaks for b- and y- ions with charge 1 (sorted by m/z)

            // total / complete loss spectra are generated for fast and (slow) full scoring
            if (total_loss_spectrum_z1.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              total_loss_spectrum_generator.getSpectrum(total_loss_spectrum_z1, fixed_and_variable_modified_peptide, 1, 1);
              total_loss_spectrum_generator.getSpectrum(total_loss_spectrum_z2, fixed_and_variable_modified_peptide, 1, 2);
              immonium_ion_sub_score_spectrum_generator.getSpectrum(immonium_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
              ONuXLFragmentIonGenerator::addSpecialLysImmonumIons(
                unmodified_sequence, 
                immonium_sub_score_spectrum, 
                immonium_sub_score_spectrum.getIntegerDataArrays()[0], 
                immonium_sub_score_spectrum.getStringDataArrays()[0]);
              immonium_sub_score_spectrum.sortByPosition();
              precursor_ion_sub_score_spectrum_generator.getSpectrum(precursor_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
              a_ion_sub_score_spectrum_generator.getSpectrum(a_ion_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
            }

            if (!fast_scoring_)
            {
              PeakSpectrum marker_ions_sub_score_spectrum_z1;
              //shifted_immonium_ions_sub_score_spectrum;
              PeakSpectrum partial_loss_spectrum_z1, partial_loss_spectrum_z2;

              // retrieve RNA adduct name
              auto mod_combinations_it = mm.mod_combinations.begin();
              std::advance(mod_combinations_it, rna_mod_index);
              const String& precursor_rna_adduct = *mod_combinations_it->second.begin();

              if (precursor_rna_adduct == "none")
              {
                // score peptide without RNA (same method as fast scoring)
                for (auto l = low_it; l != up_it; ++l)
                {
                  //const double exp_pc_mass = l->first;
                  const Size & scan_index = l->second.first;
                  const int & isotope_error = l->second.second;
                  const PeakSpectrum & exp_spectrum = spectra[scan_index];
                  const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                  PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                  float total_loss_score(0), 
                        immonium_sub_score(0), 
                        precursor_sub_score(0), 
                        a_ion_sub_score(0), 
                        tlss_MIC(0),
                        tlss_err(0), 
                        tlss_Morph(0);

                  scoreTotalLossFragments_(exp_spectrum,
                                         total_loss_spectrum,
                                         fragment_mass_tolerance,
                                         fragment_mass_tolerance_unit_ppm,
                                         a_ion_sub_score_spectrum,
                                         precursor_sub_score_spectrum,
                                         immonium_sub_score_spectrum,
                                         total_loss_score,
                                         tlss_MIC,
                                         tlss_err,
                                         tlss_Morph,
                                         immonium_sub_score,
                                         precursor_sub_score,
                                         a_ion_sub_score);


                  // no good hit
                  if (total_loss_score < 0.001) { continue; }

                  // add peptide hit
                  AnnotatedHit_ ah;
                  ah.sequence = *cit; // copy StringView
                  ah.peptide_mod_index = mod_pep_idx;
                  ah.MIC = tlss_MIC;
                  ah.err = tlss_err;
                  ah.Morph = tlss_Morph;
                  ah.immonium_score = immonium_sub_score;
                  ah.precursor_score = precursor_sub_score;
                  ah.a_ion_score = a_ion_sub_score;
                  ah.total_MIC = tlss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                  ah.rna_mod_index = rna_mod_index;
                  ah.isotope_error = isotope_error;

                  ah.score = total_loss_score + ah.total_MIC; 

#ifdef DEBUG_ONUXLSEARCH
                  LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP 
                  omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                  {
                    annotated_hits[scan_index].push_back(ah);

                    // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                    if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                    {
                      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
                      annotated_hits[scan_index].resize(report_top_hits); 
                    }
                  }
#ifdef _OPENMP 
                  omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
                }
              }
              else  // score peptide with RNA adduct
              {
                // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
                // get RNA fragment shifts in the MS2 (based on the precursor RNA/DNA)
                const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_fragment_adducts.at(precursor_rna_adduct).feasible_adducts;

                // get marker ions
                const vector<ONuXLFragmentAdduct>& marker_ions = all_feasible_fragment_adducts.at(precursor_rna_adduct).marker_ions;

                //cout << "'" << precursor_rna_adduct << "'" << endl;
                //OPENMS_POSTCONDITION(!feasible_MS2_adducts.empty(),
                //                String("FATAL: No feasible adducts for " + precursor_rna_adduct).c_str());


                // Do we have (nucleotide) specific fragmentation adducts? for the current RNA adduct on the precursor?
                // If so, generate spectra for shifted ion series

                // score individually for every nucleotide
                for (auto const & nuc_2_adducts : feasible_MS2_adducts)
                {
                  char cross_linked_nucleotide = nuc_2_adducts.first;
                  const vector<ONuXLFragmentAdduct>& partial_loss_modification = nuc_2_adducts.second;

                  if (!partial_loss_modification.empty())
                  {
                    // shifted b- / y- / a-ions
                    // generate shifted_immonium_ions_sub_score_spectrum.empty
                    ONuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                                fixed_and_variable_modified_peptide,
                                                current_peptide_mass_without_RNA,
                                                precursor_rna_adduct,
                                                precursor_rna_weight,
                                                1,
                                                partial_loss_modification,
                                                partial_loss_spectrum_generator,
                                                partial_loss_spectrum_z1);
                    for (auto& n : partial_loss_spectrum_z1.getStringDataArrays()[0]) { n[0] = 'y'; } // hyperscore hack

                    ONuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                                fixed_and_variable_modified_peptide,
                                                current_peptide_mass_without_RNA,
                                                precursor_rna_adduct,
                                                precursor_rna_weight,
                                                2, // don't know the charge of the precursor at that point
                                                partial_loss_modification,
                                                partial_loss_spectrum_generator,
                                                partial_loss_spectrum_z2);
                    for (auto& n : partial_loss_spectrum_z2.getStringDataArrays()[0]) { n[0] = 'y'; } // hyperscore hack
                  }

                  // add shifted marker ions
                  marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
                  marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation
                  ONuXLFragmentIonGenerator::addMS2MarkerIons(
                    marker_ions,
                    marker_ions_sub_score_spectrum_z1,
                    marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[0],
                    marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

                  for (auto l = low_it; l != up_it; ++l)
                  {
                    //const double exp_pc_mass = l->first;
                    const Size& scan_index = l->second.first;
                    const int& isotope_error = l->second.second;
                    const PeakSpectrum& exp_spectrum = spectra[scan_index];
                    float tlss_MIC(0), tlss_err(0), tlss_Morph(0),
                      immonium_sub_score(0), precursor_sub_score(0),
                      a_ion_sub_score(0), partial_loss_sub_score(0), marker_ions_sub_score(0),
                      plss_MIC(0), plss_err(0), plss_Morph(0), score;

                    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                    PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                    scoreTotalLossFragments_(exp_spectrum,
                                             total_loss_spectrum,
                                             fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                             a_ion_sub_score_spectrum,
                                             precursor_sub_score_spectrum,
                                             immonium_sub_score_spectrum,
                                             score,
                                             tlss_MIC,
                                             tlss_err,
                                             tlss_Morph,
                                             immonium_sub_score,
                                             precursor_sub_score,
                                             a_ion_sub_score);

                    scorePartialLossFragments_(exp_spectrum,
                                               fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                               partial_loss_spectrum_z1, partial_loss_spectrum_z2,
                                               marker_ions_sub_score_spectrum_z1,
                                               partial_loss_sub_score,
                                               marker_ions_sub_score,
                                               plss_MIC, plss_err, plss_Morph);

                    // no good hit
                    if (score < 0.001) { continue; }

                    // add peptide hit
                    AnnotatedHit_ ah;
                    ah.sequence = *cit; // copy StringView
                    ah.peptide_mod_index = mod_pep_idx;
                    ah.MIC = tlss_MIC;
                    ah.err = tlss_err;
                    ah.Morph = tlss_Morph;
                    ah.pl_MIC = plss_MIC;
                    ah.pl_err = plss_err;
                    ah.pl_Morph = plss_Morph;
                    ah.immonium_score = immonium_sub_score;
                    ah.precursor_score = precursor_sub_score;
                    ah.a_ion_score = a_ion_sub_score;
                    ah.cross_linked_nucleotide = cross_linked_nucleotide;
                    ah.total_MIC = tlss_MIC + plss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                    // scores from shifted peaks
                    ah.marker_ions_score = marker_ions_sub_score;
                    ah.partial_loss_score = partial_loss_sub_score;

                    ah.rna_mod_index = rna_mod_index;
                    ah.isotope_error = isotope_error;

// TODO: currently mainly a tie-breaker!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ah.score = score + ah.total_MIC; 
#ifdef DEBUG_ONUXLSEARCH
                    LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
                    omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                    {
                      annotated_hits[scan_index].push_back(ah);

                      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                      if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                      {
                        std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
                        annotated_hits[scan_index].resize(report_top_hits); 
                      }
                    }
#ifdef _OPENMP
                    omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
                  }
                } // for every nucleotide in the precursor
              }
            }
            else // fast scoring
            {
              for (auto l = low_it; l != up_it; ++l)
              {
                //const double exp_pc_mass = l->first;
                const Size &scan_index = l->second.first;
                const int &isotope_error = l->second.second;
                const PeakSpectrum &exp_spectrum = spectra[scan_index];
                float total_loss_score;
                float immonium_sub_score;
                float precursor_sub_score;
                float a_ion_sub_score;
                float tlss_MIC;
                float tlss_err;
                float tlss_Morph;

                const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                scoreTotalLossFragments_(exp_spectrum, 
                                         total_loss_spectrum, 
                                         fragment_mass_tolerance,
                                         fragment_mass_tolerance_unit_ppm, 
                                         a_ion_sub_score_spectrum,
                                         precursor_sub_score_spectrum, 
                                         immonium_sub_score_spectrum, 
                                         total_loss_score, 
                                         tlss_MIC,
                                         tlss_err,
                                         tlss_Morph,
                                         immonium_sub_score,
                                         precursor_sub_score,
                                         a_ion_sub_score);

                // no good hit
                if (total_loss_score < 0.001) { continue; }

                // add peptide hit
                AnnotatedHit_ ah;
                ah.sequence = *cit; // copy StringView
                ah.peptide_mod_index = mod_pep_idx;
                ah.MIC = tlss_MIC;
                ah.err = tlss_err;
                ah.Morph = tlss_Morph;
                ah.immonium_score = immonium_sub_score;
                ah.precursor_score = precursor_sub_score;
                ah.a_ion_score = a_ion_sub_score;

                ah.total_MIC = tlss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                ah.rna_mod_index = rna_mod_index;
                ah.isotope_error = isotope_error;

// TODO: currently mainly a tie-breaker!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ah.score = total_loss_score + ah.total_MIC; 

#ifdef DEBUG_ONUXLSEARCH
                LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
                omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                {
                  annotated_hits[scan_index].push_back(ah);

                  // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                  if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                  {
                    std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
                    annotated_hits[scan_index].resize(report_top_hits); 
                  }
                }
#ifdef _OPENMP
                omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
              }
            }
          }
        }
      }
    }
    endProgress();

    LOG_INFO << "Proteins: " << count_proteins << endl;
    LOG_INFO << "Peptides: " << count_peptides << endl;
    LOG_INFO << "Processed peptides: " << processed_petides.size() << endl;

    startProgress(0, 1, "Post-processing PSMs...");

    if (localization_)
    {
      // reload spectra from disc with same settings as before (important to keep same spectrum indices)
      spectra = in_spectra;

      // for post scoring don't convert fragments to single charge. Annotate charge instead to every peak.
      ONuXLSpectrumProcessing::preprocess(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, true); // no single charge (false), annotate charge (true)

      startProgress(0, 1, "localization...");

      // create spectrum generator. For convenience we add more peak types here.
      Param param(total_loss_spectrum_generator.getParameters());
      param.setValue("add_first_prefix_ion", "true");
      param.setValue("add_abundant_immonium_ions", "true");
      param.setValue("add_precursor_peaks", "true");
      param.setValue("add_metainfo", "true");
      param.setValue("add_a_ions", "false");
      param.setValue("add_b_ions", "true");
      param.setValue("add_c_ions", "false");
      param.setValue("add_x_ions", "false");
      param.setValue("add_y_ions", "true");
      param.setValue("add_z_ions", "false");
      total_loss_spectrum_generator.setParameters(param);

      postScoreHits_(spectra, 
                     annotated_hits, 
                     report_top_hits, 
                     mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, 
                     partial_loss_spectrum_generator, 
                     fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
                     all_feasible_fragment_adducts);
    }

    startProgress(0, 1, "Post-processing and annotation...");
    postProcessHits_(spectra, 
                     annotated_hits, 
                     protein_ids, peptide_ids, 
                     report_top_hits, 
                     mm, 
                     fixed_modifications, variable_modifications, 
                     max_variable_mods_per_peptide,
                     database_filename);
    endProgress();

    // annotate ONuXL related information to hits and create report
    vector<ONuXLReportRow> csv_rows = ONuXLReport::annotate(spectra, peptide_ids, marker_ions_tolerance);

    // reindex ids
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", (String)param_.getValue("peptide:enzyme"));
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
      {
        return ONuXLSearchEngine::ExitCodes::INPUT_FILE_EMPTY;       
      }
      else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
      {
        return ONuXLSearchEngine::ExitCodes::UNEXPECTED_RESULT;
      }
      else
      {
        return ONuXLSearchEngine::ExitCodes::UNKNOWN_ERROR;
      }
    } 

    if (generate_decoys)	
    {
      // calculate FDR
      FalseDiscoveryRate fdr;     	
      fdr.apply(peptide_ids);	
    }

    // store report
    csv_rows = ONuXLReport::annotate(spectra, peptide_ids, marker_ions_tolerance);
    csv_file.addLine(ONuXLReportRowHeader().getString("\t"));
    for (Size i = 0; i != csv_rows.size(); ++i)
    {
      csv_file.addLine(csv_rows[i].getString("\t"));
    }      
  
 #ifdef _OPENMP
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_destroy_lock(&(annotated_hits_lock[i])); } // free locks
 #endif
    return ExitCodes::EXECUTION_OK;
  }


  // determine main score and sub scores of peaks without shifts
  void ONuXLSearchEngine::scoreTotalLossFragments_(const PeakSpectrum &exp_spectrum,
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
                                float &a_ion_sub_score) const
  {
    total_loss_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                           exp_spectrum, total_loss_spectrum);
    immonium_sub_score = 0;
    precursor_sub_score = 0;
    a_ion_sub_score = 0;
    auto const & tl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                       fragment_mass_tolerance_unit_ppm,
                                                       exp_spectrum,
                                                       total_loss_spectrum);

    tlss_MIC = tl_sub_scores.TIC != 0 ? tl_sub_scores.MIC / tl_sub_scores.TIC : 0;
    tlss_err = tl_sub_scores.err;
    tlss_Morph = tl_sub_scores.score;

    if (!immonium_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             immonium_sub_score_spectrum);
      immonium_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    if (!precursor_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             precursor_sub_score_spectrum);
      precursor_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    if (!a_ion_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             a_ion_sub_score_spectrum);
      a_ion_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
  }

  void ONuXLSearchEngine::scorePartialLossFragments_(const PeakSpectrum &exp_spectrum,
                                  double fragment_mass_tolerance,
                                  bool fragment_mass_tolerance_unit_ppm,
                                  const PeakSpectrum &partial_loss_spectrum_z1,
                                  const PeakSpectrum &partial_loss_spectrum_z2,
                                  const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                                  float &partial_loss_sub_score,
                                  float &marker_ions_sub_score,
                                  float &plss_MIC, float &plss_err, float &plss_Morph) const
  {
    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    plss_MIC = 0;
    plss_err = 0;
    plss_Morph = 0;


    if (!marker_ions_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             marker_ions_sub_score_spectrum_z1);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    //TODO: these are currently empty
/*
              if (!shifted_immonium_ions_sub_score_spectrum.empty())
              {
                shifted_immonium_ions_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, shifted_immonium_ions_sub_score_spectrum);
                               shifted_immonium_ions_sub_score(0),
              }
*/
    if (!partial_loss_spectrum_z1.empty()) // check if we generated partial loss spectra
    {
      if (exp_pc_charge < 3)
      {
        partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                                     exp_spectrum, partial_loss_spectrum_z1);
        auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                           fragment_mass_tolerance_unit_ppm,
                                                           exp_spectrum,
                                                           partial_loss_spectrum_z1);

        plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
        plss_err = pl_sub_scores.err;
        plss_Morph = pl_sub_scores.score;
      }
      else //if (exp_pc_charge >= 3)
      {
        partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                                     exp_spectrum, partial_loss_spectrum_z2);
        auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                           fragment_mass_tolerance_unit_ppm,
                                                           exp_spectrum,
                                                           partial_loss_spectrum_z2);

        plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
        plss_err = pl_sub_scores.err;
        plss_Morph = pl_sub_scores.score;
      }
    }
#ifdef DEBUG_ONUXLSEARCH
    LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
  }


    // fast or all-shifts scoring mode
    bool fast_scoring_ = true;

    // if localization should be performed
    bool localization_ = false;

    // nucleotides can form cross-link
    set<char> can_xl_;

  /* @brief Localization step of the cross-link identification engine.
   * Given a top scoring candidate (based on total loss spectrum) it:
   *  - generates all fragment adducts based on the attached precursor adduct
   *  - annotated peaks
   *  - calculates an additive score that considers the presence or absence of evidence for a cross-linking site
   *  - the maximum score is reported
   */
  void ONuXLSearchEngine::postScoreHits_(const PeakMap& exp, 
                      vector<vector<AnnotatedHit_> >& annotated_hits, 
                      Size top_hits, 
                      const ONuXLModificationMassesResult& mm, 
                      const vector<ResidueModification>& fixed_modifications, 
                      const vector<ResidueModification>& variable_modifications, 
                      Size max_variable_mods_per_peptide, 
                      const TheoreticalSpectrumGenerator& partial_loss_spectrum_generator, 
                      double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, 
                      const ONuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    assert(exp.size() == annotated_hits.size());

    #ifdef DEBUG_ONUXLSEARCH
      LOG_DEBUG << exp.size() << " : " << annotated_hits.size() << endl;
    #endif

    SpectrumAlignment spectrum_aligner;
    Param pa = spectrum_aligner.getParameters();
    pa.setValue("tolerance", fragment_mass_tolerance, "Defines the absolute (in Da) or relative (in ppm) tolerance in the alignment");
    pa.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm ? "true" : "false");  
    spectrum_aligner.setParameters(pa);

    // remove all but top n scoring for localization (usually all but the first one)
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      const Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
      annotated_hits[scan_index].shrink_to_fit();
    }

    // If we did a (total-loss) only fast scoring, PSMs were not associated with a nucleotide.
    // To make the localization code work for both fast and slow (all-shifts) scoring,
    // we copy PSMs for every cross-linkable nucleotide present in the precursor.
    if (fast_scoring_)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
      {
        vector<AnnotatedHit_> new_hits;

        // for each PSM
        for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
        {
          // determine RNA on precursor from index in map
          auto mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, annotated_hits[scan_index][i].rna_mod_index);
          const String precursor_rna_adduct = *mod_combinations_it->second.begin();
          const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;

          // copy PSM information for each cross-linkable nucleotides
          for (auto const & c : feasible_MS2_adducts)
          {
            AnnotatedHit_ a(annotated_hits[scan_index][i]);
            a.cross_linked_nucleotide = c.first; // nucleotide
            new_hits.push_back(a);
          }
        }
        annotated_hits[scan_index].swap(new_hits);
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) { continue; }

      const PeakSpectrum & exp_spectrum = exp[scan_index];
      const Size & precursor_charge = exp_spectrum.getPrecursors()[0].getCharge();

      for (auto & a : annotated_hits[scan_index])
      {
        // get unmodified string
        const String unmodified_sequence = a.sequence.getString();

        // initialize result fields
        a.best_localization = unmodified_sequence;
        a.best_localization_score = 0;

        AASequence aas(AASequence::fromString(unmodified_sequence));

        // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);

        // sequence with modifications - note: reannotated version requires much more memory heavy AASequence object
        const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[a.peptide_mod_index];
        const double fixed_and_variable_modified_peptide_weight = fixed_and_variable_modified_peptide.getMonoWeight();

        // determine RNA on precursor from index in map
        std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, a.rna_mod_index);
        const String precursor_rna_adduct = *mod_combinations_it->second.begin();
        const double precursor_rna_weight = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

        // we don't localize on non-cross-links
        if (precursor_rna_adduct == "none") { continue; }

        // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
        // 1. get all possible RNA fragment shifts in the MS2 (based on the precursor RNA/DNA)
        LOG_DEBUG << "precursor_rna_adduct: "  << precursor_rna_adduct << endl;
        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;

        if (feasible_MS2_adducts.empty()) { continue; } // should not be the case - check case of no nucleotide but base fragment ?

        // 2. retrieve the (nucleotide specific) fragment adducts for the cross-linked nucleotide (annotated in main search)
        auto nt_to_adducts = std::find_if(feasible_MS2_adducts.begin(), feasible_MS2_adducts.end(),
          [&a](NucleotideToFeasibleFragmentAdducts const & item)
          {
            return (item.first == a.cross_linked_nucleotide);
          });

        OPENMS_POSTCONDITION(nt_to_adducts != feasible_MS2_adducts.end(), "Nucleotide not found in mapping to feasible adducts.")

        const vector<ONuXLFragmentAdduct>& partial_loss_modification = nt_to_adducts->second;

        // get marker ions (these are not specific to the cross-linked nucleotide but also depend on the whole oligo bound to the precursor)
        const vector<ONuXLFragmentAdduct>& marker_ions = all_feasible_adducts.at(precursor_rna_adduct).marker_ions;

        // generate total loss spectrum for the fixed and variable modified peptide (without RNA) (using the settings for partial loss generation)
        // but as we also add the abundant immonium ions for charge 1 and precursor ions for all charges to get a more complete annotation
        // (these have previously not been used in the scoring of the total loss spectrum)
        PeakSpectrum total_loss_spectrum;

        TheoreticalSpectrumGenerator tmp_generator;
        Param new_param(partial_loss_spectrum_generator.getParameters());
        new_param.setValue("add_all_precursor_charges", "true");
        new_param.setValue("add_abundant_immonium_ions", "true");
        tmp_generator.setParameters(new_param);
        tmp_generator.getSpectrum(total_loss_spectrum, fixed_and_variable_modified_peptide, 1, precursor_charge);

        // add special immonium ions
        ONuXLFragmentIonGenerator::addSpecialLysImmonumIons(
          unmodified_sequence,
          total_loss_spectrum, 
          total_loss_spectrum.getIntegerDataArrays()[0],
          total_loss_spectrum.getStringDataArrays()[0]);
        total_loss_spectrum.sortByPosition(); // need to resort after adding special immonium ions

        PeakSpectrum partial_loss_spectrum;
        ONuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                    fixed_and_variable_modified_peptide,
                                    fixed_and_variable_modified_peptide_weight,
                                    precursor_rna_adduct,
                                    precursor_rna_weight,
                                    precursor_charge,
                                    partial_loss_modification,
                                    partial_loss_spectrum_generator,
                                    partial_loss_spectrum);

         // add shifted marker ions
         ONuXLFragmentIonGenerator::addMS2MarkerIons(
           marker_ions,
           partial_loss_spectrum,
           partial_loss_spectrum.getIntegerDataArrays()[0],
           partial_loss_spectrum.getStringDataArrays()[0]);

        partial_loss_spectrum.sortByPosition(); // need to resort after adding marker ions

        // fill annotated spectrum information
        set<Size> peak_is_annotated;  // experimental peak index

        // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)
        using MapIonIndexToFragmentAnnotation = map<Size, vector<ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_> >;
        MapIonIndexToFragmentAnnotation unshifted_b_ions, unshifted_y_ions, unshifted_a_ions, shifted_b_ions, shifted_y_ions, shifted_a_ions;
        vector<PeptideHit::PeakAnnotation> shifted_immonium_ions;
        vector<PeptideHit::PeakAnnotation> annotated_marker_ions;
        vector<PeptideHit::PeakAnnotation> annotated_precursor_ions;
        vector<PeptideHit::PeakAnnotation> annotated_immonium_ions;

        // first annotate total loss peaks (these give no information where the actual shift occured)
        #ifdef DEBUG_ONUXLSEARCH
          LOG_DEBUG << "Annotating ion (total loss spectrum): " << fixed_and_variable_modified_peptide.toString()  << endl;
        #endif
        vector<pair<Size, Size>> alignment;
        spectrum_aligner.getSpectrumAlignment(alignment, total_loss_spectrum, exp_spectrum);

        const PeakSpectrum::StringDataArray& total_loss_annotations = total_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& total_loss_charges = total_loss_spectrum.getIntegerDataArrays()[0];

        for (auto const & aligned : alignment)
        {
          // information on the experimental fragment in the alignment
          const Size& fragment_index = aligned.second;
          const Peak1D& fragment = exp_spectrum[fragment_index];
          const double fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double fragment_mz = fragment.getMZ();
          const int fragment_charge = exp_spectrum.getIntegerDataArrays().back()[fragment_index];

          const String& ion_name = total_loss_annotations[aligned.first];
          const int charge = total_loss_charges[aligned.first];

          // define which ion names are annotated
          if (ion_name.hasPrefix("y"))
          {
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("y", "");
              ion_nr_string.substitute("+", "");
              auto ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_ONUXLSEARCH
              const AASequence& peptide_sequence = fixed_and_variable_modified_peptide.getSuffix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_y_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("b"))
          {
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("b", "");
              ion_nr_string.substitute("+", "");
              auto ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_ONUXLSEARCH
              const AASequence& peptide_sequence = aas.getPrefix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_b_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("a"))
          {
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("a", "");
              ion_nr_string.substitute("+", "");
              auto ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_ONUXLSEARCH
              const AASequence& peptide_sequence = aas.getPrefix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_a_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("[M+")) // precursor ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1; // for visualion charge is not really important so we set it to 0
            fa.annotation = ion_name;
            peak_is_annotated.insert(aligned.second);
            annotated_precursor_ions.push_back(fa);
          }
          else if (ion_name.hasPrefix("i")) // immonium ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1;
            fa.annotation = ion_name;
            peak_is_annotated.insert(aligned.second);
            annotated_immonium_ions.push_back(fa);
          }
        }

        // generate fragment annotation strings for unshifted ions
        vector<PeptideHit::PeakAnnotation> fas;
        if (!unshifted_b_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("b", unshifted_b_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!unshifted_y_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("y", unshifted_y_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!unshifted_a_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("a", unshifted_a_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!annotated_immonium_ions.empty())
        {
          fas.insert(fas.end(), annotated_immonium_ions.begin(), annotated_immonium_ions.end());          
        }

        vector<double> sites_sum_score(aas.size(), 0);

        /////////////////
        // Align partial-loss-spectrum to the experimental measured one
        alignment.clear();

        spectrum_aligner.getSpectrumAlignment(alignment, partial_loss_spectrum, exp_spectrum);

        const PeakSpectrum::StringDataArray& partial_loss_annotations = partial_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& partial_loss_charges = partial_loss_spectrum.getIntegerDataArrays()[0];

        if (alignment.empty())
        {
          a.fragment_annotations = fas;
          continue;
        }

        for (vector<std::pair<Size, Size> >::const_iterator pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
        {
          // only annotate experimental peaks with shift - i.e. do not annotated complete loss peaks again
          if (peak_is_annotated.find(pair_it->second) != peak_is_annotated.end()) { continue; }

          // information on the experimental fragment in the alignment
          const Size & fragment_index = pair_it->second;
          const Peak1D & fragment = exp_spectrum[fragment_index];
          const double & fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double & fragment_mz = fragment.getMZ();
          const int & fragment_charge = exp_spectrum.getIntegerDataArrays().back()[fragment_index];

          String ion_name = partial_loss_annotations[pair_it->first];
          const int charge = partial_loss_charges[pair_it->first];

          vector<String> f;

          ion_name.split(' ', f);  // e.g. "y3 C3O" or just "y2"
          String fragment_shift_name;
          if (f.size() == 2) { fragment_shift_name = f[1]; }

          String fragment_ion_name = f[0]; // e.g. y3

          #ifdef DEBUG_ONUXLSEARCH
            LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
          #endif

          // define which ion names are annotated
          if (fragment_ion_name.hasPrefix("y"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("y", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_y_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (fragment_ion_name.hasPrefix("b"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("b", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_b_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (fragment_ion_name.hasPrefix("a"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("a", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              ONuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_a_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix(ONuXLFragmentIonGenerator::ANNOTATIONS_MARKER_ION_PREFIX))
          {
            if (fragment_charge <= 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              annotated_marker_ions.push_back(fa);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << 1 << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("i"))
          {
            if (fragment_charge <= 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              shifted_immonium_ions.push_back(fa);
            }
            #ifdef DEBUG_ONUXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << 1 << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("[M+"))
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1; // for visualion charge is not really important so we set it to 1, TODO: read out charge from ion name and set here
            fa.annotation = ion_name;
            annotated_precursor_ions.push_back(fa);
          }
        }

        // track shifts in n- and c-term ladders (in AAs coordinates)
        // n_shifts and c_shifts will contain the summed intensities over all observed shifts at that position
        // the distinction allows to easily detect prefix and suffix ladders in the next step
        vector<double> n_shifts(sites_sum_score.size(), 0);
        vector<double> c_shifts(sites_sum_score.size(), 0);

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_b_ions.find(i + 1) == shifted_b_ions.end()) { continue; }
          for (auto& k : shifted_b_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_a_ions.find(i + 1) == shifted_a_ions.end()) { continue; }
          for (auto& k : shifted_a_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_shifts.size(); ++i)
        {
          const Size ion_index = c_shifts.size() - i;
          if (shifted_y_ions.find(ion_index) == shifted_y_ions.end()) { continue; }
          for (auto& k : shifted_y_ions[ion_index]) { c_shifts[i] += k.intensity; }
        }

        vector<double> n_noshifts(sites_sum_score.size(), 0);
        vector<double> c_noshifts(sites_sum_score.size(), 0);
        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_b_ions.find(i + 1) == unshifted_b_ions.end()) { continue; }
          for (auto& k : unshifted_b_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_a_ions.find(i + 1) == unshifted_a_ions.end()) { continue; }
          for (auto& k : unshifted_a_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_noshifts.size(); ++i)
        {
          const Size ion_index = c_noshifts.size() - i;
          if (unshifted_y_ions.find(ion_index) == unshifted_y_ions.end()) { continue; }
          for (auto& k : unshifted_y_ions[ion_index]) { c_noshifts[i] += k.intensity; }
        }

#ifdef DEBUG_ONUXLSEARCH
        cout << "n:";
        for (auto& k : n_shifts) cout << k << " ";
        cout << endl;
        cout << "c:";
        for (auto& k : c_shifts) cout << k << " ";
        cout << endl;
        cout << "n0:";
        for (auto& k : n_noshifts) cout << k << " ";
        cout << endl;
        cout << "c0:";
        for (auto& k : c_noshifts) cout << k << " ";
        cout << endl;
#endif

        // Rules implemented:
        // 1. if cross-link on AA, then the prefix or suffix ending at this AA must be shifted
        // 2. if the previous AA in the prefix / suffix had a stronger shifted signal, then the current on is not the correct one
        // 3. if the current AA is cross-linked, then the previous AA is not cross-linked and we should observe an unshifted prefix / suffix ion
        // 4. if the current AA is cross-linked, we should observe a shifted prefix / suffix ion for the next AA, too
        // 5. Sum up all intensities of shifted prefix / suffix ions
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          sites_sum_score[i] = 0.0;
          if (n_shifts[i] == 0 && c_shifts[i] == 0) { continue; } // no shifts? no cross-link at this AA

          if (n_shifts[i] > 0)
          {
            if (i >= 1 && n_shifts[i - 1] > n_shifts[i]) continue; // Strong signal from shifted AA before the current one? Then skip it.
            if (i >= 1 && n_noshifts[i - 1] == 0) continue; // continue if unshifted AA is missing before (left of) the shifted one.
            if (i < n_shifts.size()-1 && n_shifts[i + 1] == 0) continue; // Need a shifted ladder after (maybe too conservative?)
            // sum up all intensities from this position and all longer prefixes that also carry the NA
            for (Size j = i; j != sites_sum_score.size(); ++j) { sites_sum_score[i] += n_shifts[j]; }
          }

          if (c_shifts[i] > 0)
          {
            if (i < c_shifts.size()-1 && c_shifts[i + 1] > c_shifts[i]) continue; // AA after already shifted. So ignore this one.
            if (i < c_noshifts.size()-1 && c_noshifts[i + 1] == 0) continue; // continue if unshifted AA is missing before (right of) the shifted one.
            if (i >=1 && c_shifts[i - 1] == 0) continue; // Need a shifted ladder after (maybe too conservative?)
            sites_sum_score[i] += c_shifts[i];
            // sum up all intensities from this position and all longer suffixes that also carry the NA
            for (int j = i; j >= 0; --j) { sites_sum_score[i] += c_shifts[j]; }
          }
        }
#ifdef DEBUG_ONUXLSEARCH
        cout << "site sum score (shifted a/b/y-ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        #ifdef DEBUG_ONUXLSEARCH
          LOG_DEBUG << "Localisation based on immonium ions: ";
        #endif
        String aas_unmodified = aas.toUnmodifiedString();
        for (Size i = 0; i != aas_unmodified.size(); ++i)
        {
          String origin = String(aas_unmodified[i]);

          for (auto& a : shifted_immonium_ions)
          {
            // compare origin (the AA) of immonium ion to current AA
            if (a.annotation[1] == aas_unmodified[i])
            {
              sites_sum_score[i] += a.intensity;
            }
          }
        }
#ifdef DEBUG_ONUXLSEARCH
        cout << "site sum score (shifted a/b/y-ions & immonium ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        String best_localization = unmodified_sequence;
        double best_localization_score = 0;
        String localization_scores;
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          if (sites_sum_score[i] > best_localization_score) { best_localization_score = sites_sum_score[i]; }
        }

        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          #ifdef DEBUG_ONUXLSEARCH
            LOG_DEBUG << String::number(100.0 * sites_sum_score[i], 2);
          #endif

          if (i != 0) localization_scores += ' ';
          if (sites_sum_score[i] > 0 )
          {
            localization_scores += String::number(100.0 * sites_sum_score[i], 2);
          }
          else
          {
            localization_scores += "0";
          }

          if (best_localization_score > 0.0 && sites_sum_score[i] >= best_localization_score - 1e-6)
          {
            best_localization[i] = tolower(best_localization[i]);
          }
        }
        #ifdef DEBUG_ONUXLSEARCH
          LOG_DEBUG << endl;
        #endif

        // create annotation strings for shifted fragment ions
        ONuXLFragmentAnnotationHelper::addShiftedPeakFragmentAnnotation_(shifted_b_ions,
                                          shifted_y_ions,
                                          shifted_a_ions,
                                          shifted_immonium_ions,
                                          annotated_marker_ions,
                                          annotated_precursor_ions,
                                          fas);

        // store score of best localization(s)
        a.localization_scores = localization_scores;
        a.best_localization = best_localization;
        a.best_localization_score = best_localization_score;
        a.fragment_annotations = fas;

        #ifdef DEBUG_ONUXLSEARCH
          LOG_DEBUG << "Ion centric annotation: " << endl;
          LOG_DEBUG << "unshifted b ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", unshifted_b_ions) << endl;
          LOG_DEBUG << "unshifted y ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", unshifted_y_ions) << endl;
          LOG_DEBUG << "unshifted a ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", unshifted_a_ions) << endl;
          LOG_DEBUG << "shifted b ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", shifted_b_ions) << endl;
          LOG_DEBUG << "shifted y ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", shifted_y_ions) << endl;
          LOG_DEBUG << "shifted a ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", shifted_a_ions) << endl;
          LOG_DEBUG << "shifted immonium ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::shiftedIonsToString(shifted_immonium_ions) << endl;
          LOG_DEBUG << "shifted marker ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::shiftedIonsToString(annotated_marker_ions) << endl;
          LOG_DEBUG << "shifted precursor ions: " << endl;
          LOG_DEBUG << ONuXLFragmentAnnotationHelper::shiftedIonsToString(annotated_precursor_ions) << endl;
          LOG_DEBUG << "Localization scores: ";
          LOG_DEBUG << localization_scores << endl;
          LOG_DEBUG << "Localisation based on ion series and immonium ions of all observed fragments: ";
          LOG_DEBUG << best_localization << endl;
        #endif
      }
    }
  }

  /// Filter by top scoring hits, reconstruct original peptide from memory efficient structure, and add additional meta information.
  void ONuXLSearchEngine::postProcessHits_(const PeakMap& exp, 
    vector<vector<AnnotatedHit_> >& annotated_hits, 
    vector<ProteinIdentification>& protein_ids, 
    vector<PeptideIdentification>& peptide_ids, 
    Size top_hits, 
    const ONuXLModificationMassesResult& mm, 
    const vector<ResidueModification>& fixed_modifications, 
    const vector<ResidueModification>& variable_modifications, 
    Size max_variable_mods_per_peptide,
    const String& database_filename)
  {
      // remove all but top n scoring (Note: this is currently necessary as postScoreHits_ might reintroduce nucleotide specific hits for fast scoring)
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
        PeptideIdentification pi;
        pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        pi.setScoreType("hyperscore");
        pi.setHigherScoreBetter(true);
        pi.setRT(exp[scan_index].getRT());
        pi.setMZ(exp[scan_index].getPrecursors()[0].getMZ());
        pi.setMetaValue("precursor_intensity", exp[scan_index].getPrecursors()[0].getIntensity());
        Size charge = exp[scan_index].getPrecursors()[0].getCharge();

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        size_t rank(0);
        for (auto const & ah : annotated_hits[scan_index])
        {
          PeptideHit ph;
          ph.setCharge(charge);

          // get unmodified string
          const String & s = ah.sequence.getString();

          OPENMS_POSTCONDITION(!s.empty(), "Error: empty sequence in annotated hits.");
          AASequence aas = AASequence::fromString(s);

          // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
          vector<AASequence> all_modified_peptides;
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);

          // reannotate much more memory heavy AASequence object
          AASequence fixed_and_variable_modified_peptide = all_modified_peptides[ah.peptide_mod_index]; 
          ph.setScore(ah.score);
          ph.setMetaValue(String("RNPxl:total_loss_score"), ph.getScore()); // important for Percolator feature set

          // determine RNA modification from index in map
          std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, ah.rna_mod_index);
          ph.setMetaValue(String("RNPxl:immonium_score"), ah.immonium_score);
          ph.setMetaValue(String("RNPxl:precursor_score"), ah.precursor_score);
          ph.setMetaValue(String("RNPxl:a_ion_score"), ah.a_ion_score);
          ph.setMetaValue(String("RNPxl:marker_ions_score"), ah.marker_ions_score);
          ph.setMetaValue(String("RNPxl:partial_loss_score"), ah.partial_loss_score);

          // total loss and partial loss (pl) related subscores (matched ion current, avg. fragment error, morpheus score)
          ph.setMetaValue(String("RNPxl:MIC"), ah.MIC);
          ph.setMetaValue(String("RNPxl:err"), ah.err);
          ph.setMetaValue(String("RNPxl:Morph"), ah.Morph);
          ph.setMetaValue(String("RNPxl:pl_MIC"), ah.pl_MIC);
          ph.setMetaValue(String("RNPxl:pl_err"), ah.pl_err);
          ph.setMetaValue(String("RNPxl:pl_Morph"), ah.pl_Morph);
          ph.setMetaValue(String("RNPxl:total_MIC"), ah.total_MIC);  // fraction of matched ion current from total + partial losses

          ph.setMetaValue(String("RNPxl:RNA"), *mod_combinations_it->second.begin()); // return first nucleotide formula matching the index of the empirical formula
          ph.setMetaValue(String("RNPxl:NT"), String(ah.cross_linked_nucleotide));  // the cross-linked nucleotide
          ph.setMetaValue(String("RNPxl:RNA_MASS_z0"), EmpiricalFormula(mod_combinations_it->first).getMonoWeight()); // RNA uncharged mass via empirical formula

          ph.setMetaValue(String("RNPxl:best_localization_score"), ah.best_localization_score);
          ph.setMetaValue(String("RNPxl:localization_scores"), ah.localization_scores);
          ph.setMetaValue(String("RNPxl:best_localization"), ah.best_localization);

          ph.setPeakAnnotations(ah.fragment_annotations);
          ph.setMetaValue("isotope_error", static_cast<int>(ah.isotope_error));
          ph.setMetaValue("rank", rank);
          // set the amino acid sequence (for complete loss spectra this is just the variable and modified peptide. For partial loss spectra it additionally contains the loss induced modification)
          ph.setSequence(fixed_and_variable_modified_peptide);
          phs.push_back(ph);
          ++rank;
        }

        pi.setHits(phs);
        pi.assignRanks();

#ifdef _OPENMP
#pragma omp critical (peptide_ids_access)
#endif
        {
          peptide_ids.push_back(pi);
        }
      }
    }

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("RNPxlSearch");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_filename;
    search_parameters.charges = String((int)param_.getValue("precursor:min_charge")) + ":" + String((int)param_.getValue("precursor:max_charge"));
    search_parameters.fixed_modifications = param_.getValue("modifications:fixed").toStringList();
    search_parameters.variable_modifications = param_.getValue("modifications:variable").toStringList();
    search_parameters.missed_cleavages = (int)param_.getValue("peptide:missed_cleavages");
    search_parameters.fragment_mass_tolerance = (double)param_.getValue("fragment:mass_tolerance");
    search_parameters.precursor_mass_tolerance = (double)param_.getValue("precursor:mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = (String)param_.getValue("precursor:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.fragment_mass_tolerance_ppm = (String)param_.getValue("fragment:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme((String)param_.getValue("peptide:enzyme"));

    /* default features added in PercolatorAdapter:
     * SpecId, ScanNr, ExpMass, CalcMass, mass, 
     * peplen, charge#min..#max, enzN, enzC, enzInt, dm, absdm
     */     
    StringList feature_set;
    feature_set
       << "isotope_error"
       << "RNPxl:total_loss_score"
       << "RNPxl:immonium_score"
       << "RNPxl:precursor_score"
       << "RNPxl:a_ion_score"
       << "RNPxl:marker_ions_score"
       << "RNPxl:partial_loss_score"
       << "RNPxl:MIC"
       << "RNPxl:err"
       << "RNPxl:Morph"
       << "RNPxl:pl_MIC"
       << "RNPxl:pl_err"
       << "RNPxl:pl_Morph"
       << "RNPxl:total_MIC"
       << "RNPxl:RNA_MASS_z0";

    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));

    protein_ids[0].setSearchParameters(search_parameters);
  }

  void ONuXLSearchEngine::mapPrecursorMassesToScans_(const Int min_precursor_charge,
                                 const Int max_precursor_charge,
                                 const IntList &precursor_isotopes,
                                 const double small_peptide_mass_filter_threshold,
                                 const Size peptide_min_size,
                                 const PeakMap & spectra,
                                 multimap<double, pair<Size, int>> & multimap_mass_2_scan_index) const
  {
    Size fractional_mass_filtered(0), small_peptide_mass_filtered(0);

    for (MSExperiment::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();

        if (precursor_charge < min_precursor_charge
         || precursor_charge > max_precursor_charge)
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();

        // map (corrected) precursor mass to spectra
        for (int i : precursor_isotopes)
        {
          double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

          // corrected for monoisotopic misassignments of the precursor annotation
          if (i != 0) { precursor_mass -= i * Constants::C13C12_MASSDIFF_U; }

          if (param_.getValue("RNPxl:filter_fractional_mass") == "true")
          {
            if (precursor_mass < 1750.0 && precursor_mass - floor(precursor_mass) < 0.2)
            {
              fractional_mass_filtered++;
              continue;
            }
          }

          if (precursor_mass < small_peptide_mass_filter_threshold)
          {
            small_peptide_mass_filtered++;
            continue;
          }

          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, make_pair(scan_index, i)));
        }
      }
    }
  }

  void ONuXLSearchEngine::initializeSpectrumGenerators(TheoreticalSpectrumGenerator &total_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &partial_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &a_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &immonium_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &precursor_ion_sub_score_spectrum_generator) const
  {
    Param param(total_loss_spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    total_loss_spectrum_generator.setParameters(param);

    param = partial_loss_spectrum_generator.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    partial_loss_spectrum_generator.setParameters(param);

    // generator for sub scores for a-ions, immonium and precursor peaksparam = a_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    a_ion_sub_score_spectrum_generator.setParameters(param);
    param = immonium_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "true");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    immonium_ion_sub_score_spectrum_generator.setParameters(param);
    param = precursor_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    precursor_ion_sub_score_spectrum_generator.setParameters(param);
  }
