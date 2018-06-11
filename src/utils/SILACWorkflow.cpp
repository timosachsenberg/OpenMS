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
// $Maintainer: Nantia Leonidou $
// $Authors: Nantia Leonidou $
// --------------------------------------------------------------------------

#include <map>
#include <list>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_SILACWorkflow
    @brief Template for a new Tool
    This tool can be used for scientific stuff and more scientific applications.


    <B>The command line parameters of this tool are NANTIA SDF: </B>
    @verbinclude UTILS_SILACWorkflow.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SILACWorkflow.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class SILACWorkflow :
  public TOPPBase

{
public:

  SILACWorkflow() :
  TOPPBase("SILACWorkflow", "SILACWorkflow Tool", false)
  {
  }

protected:
  // register the tool parameters, it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    //Input mzML file list
    registerInputFileList_("in", "<files>", StringList(), "Input MzMLFile");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    //Input idxML file list (light)
    registerInputFileList_("in_ids_light", "<files>", StringList(), "Input idXML light");
    setValidFormats_("in_ids_light", ListUtils::create<String>("idXML"));

    //Input idXML file list (heavy)
    registerInputFileList_("in_ids_heavy", "<files>", StringList(), "Input idXML heavy");
    setValidFormats_("in_ids_heavy", ListUtils::create<String>("idXML"));

    //Input database fasta
    registerInputFile_("database", "<file>", "", "input database");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    // Optional experimental design file
    registerInputFile_("design", "<file>", "", "design file", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    //Output
    registerOutputFile_("out", "<file>", "", "Output consensusXML file");
    setValidFormats_("out",ListUtils::create<String>(".consensusXML"));

/*
    registerStringOption_("decoy_string", "<text>", "DECOY_", "String to indicate decoy protein", false);
    setValidStrings_("decoy_string", ListUtils::create<String>("DECOY_,dec_"));

    registerStringOption_("decoy_string_position", "<choice>", "prefix", "Should the string be prepended or appended?", false);
    setValidStrings_("decoy_string_position", ListUtils::create<String>("prefix,suffix"));

    registerStringOption_("enzyme_name", "<choice>", "Trypsin/P", "Enzyme which determines valid cleavage sites", false);
    StringList enzymes;
    ProteaseDB::getInstance()->getAllNames(enzymes);
    setValidStrings_("enzyme_name", enzymes);
*/
    //get all parameters from different algorithms at once
    Param pp_defaults = PeakPickerHiRes().getDefaults();
    Param ffm_defaults = FeatureFinderMultiplexAlgorithm().getDefaults();
    Param pep_defaults = Math::PosteriorErrorProbabilityModel().getParameters();
    Param index_defaults = PeptideIndexing().getParameters();
    Param fdr_defaults = FalseDiscoveryRate().getParameters();
    Param idMapper_defaults = IDMapper().getParameters();
    Param pq_defaults = PeptideAndProteinQuant().getDefaults();
    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("Quantification:", ffm_defaults);
    combined.insert("PosteriorErrorProbability:", pep_defaults);
    combined.insert("PeptideIndexing:", index_defaults);
    combined.insert("FalseDiscoveryRate:", fdr_defaults);
    combined.insert("IDMapper:", idMapper_defaults);
    combined.insert("Protein Quantification:", pq_defaults);
    registerFullParam_(combined);
    }

/*
***************************************************************
P e p t i d e  I d e n t i f i c a t i o n
***************************************************************
*/
protected:
  PeptideIndexing::ExitCodes indexPepAndProtIds_(
    const vector<FASTAFile::FASTAEntry>& fasta_db,
    vector<ProteinIdentification>& protein_ids,
    vector<PeptideIdentification>& peptide_ids
  )
  /**
   * [indexPepAndProtIds_ give indexes to all protein and peptide ids]
   * @param  fasta_db    [database]
   * @param  protein_ids [vector of protein IDS]
   * @param  peptide_ids [vector of peptide IDS]
   * @return
   */
   {
     Param pp_param = getParam_().copy("PeptideIndexing:", true);
     writeDebug_("Parameters passed to PeptideIndexing algorithm", pp_param, 3);
     PeptideIndexing indexer;
     indexer.setLogType(log_type_);
     indexer.setParameters(pp_param);

     vector<FASTAFile::FASTAEntry> fasta_db_tmp(fasta_db);
     PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db_tmp, protein_ids, peptide_ids);

     return indexer_exit;

/* PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string", getStringOption_("decoy_string"));
    param_pi.setValue("decoy_string_position", getStringOption_("decoy_string_position"));
    param_pi.setValue("enzyme:specificity", "none");
    param_pi.setValue("missing_decoy_action", "warn");
    param_pi.setValue("enzyme:name", getStringOption_("enzyme"));
    indexer.setParameters(param_pi);
*/

  };

protected:
  void calculatePEP_(
    vector<ProteinIdentification>& protein_identifications,
    vector<PeptideIdentification>& peptide_identifications
  )
  /**
   * [calculatePEP_ PEP calculation]
   * @param protein_identifications [vector of protein IDS]
   * @param peptide_identifications [vector of peptide IDS]
   */
  {
    //generate PEP_model
    Param pep_param = getParam_().copy("PosteriorErrorProbability:", true);
    writeDebug_("Parameters passed to PEP algorithm", pep_param, 3);
    Math::PosteriorErrorProbabilityModel PEP_model;
    PEP_model.setParameters(pep_param);

    //-------------------------------------------------------------
    // Posterior Error Probability calculation
    //-------------------------------------------------------------
    map<String, vector<vector<double> > > all_scores = Math::PosteriorErrorProbabilityModel::extractAndTransformScores(
      protein_identifications,
      peptide_identifications,
      false,
      true,
      true,
      0.05);

      for (auto & score : all_scores)  // for all search engine scores (should only be 1)
      {
        vector<String> engine_info;
        score.first.split(',', engine_info);
        String engine = engine_info[0];
        Int charge = (engine_info.size() == 2) ? engine_info[1].toInt() : -1;

        // fit to score vector
        bool return_value = PEP_model.fit(score.second[0]);
        if (!return_value)
        {
          writeLog_("Unable to fit data. Algorithm did not run through for the following search engine: " + engine);
        }

        bool unable_to_fit_data(true), data_might_not_be_well_fit(true);
        Math::PosteriorErrorProbabilityModel::updateScores(
          PEP_model,
          engine,
          charge,
          true, // prob_correct
          false, //split_charge
          protein_identifications,
          peptide_identifications,
          unable_to_fit_data,
          data_might_not_be_well_fit);

          if (unable_to_fit_data)
          {
            writeLog_(String("Unable to fit data for search engine: ") + engine);
          }
          else if (data_might_not_be_well_fit)
          {
            writeLog_(String("Data might not be well fitted for search engine: ") + engine);
          }
        }
      };

protected:
  //pair of protein and peptide identifications (vectors), as stored in an idXML file
  using ProtsPepsPair = std::pair< vector<ProteinIdentification>, vector<PeptideIdentification> >;

  //vector, of pair of IDs vectors
  using ProtsPepsPairs = vector<ProtsPepsPair>;

  ProtsPepsPairs prepareIDFiles_(
    const vector<FASTAFile::FASTAEntry>& fasta_db,
    const StringList & light,
    const StringList & heavy
  )
  /**
   * [prepareIDFiles_ prepare ID files]
   * @param  fasta_db [database]
   * @param  light
   * @param  heavy
   * @return          [vector, of pair, of IDSs vectors]
   */
  {

    //if stringList sizes not equal, throw an exception
    if (light.size() != heavy.size())
    {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not equal-sized lists");
    }

    IdXMLFile idXML_file;
    ProtsPepsPairs ret; //return statement

    //load light/heavy files and calculates PEP
    for (unsigned int i = 0; i < light.size(); i++)
    {
      /////////////////////////////////////////////////
      /////             light                    /////
      ///////////////////////////////////////////////
      vector<ProteinIdentification> light_protein_identifications;
      vector<PeptideIdentification> light_peptide_identifications;
      // read the identification file for light (contains both protein as well as peptide identifications)
      LOG_INFO << "Loading ID file (light): " << light[i] << endl;
      idXML_file.load(light[i], light_protein_identifications,  light_peptide_identifications);
      PeptideIndexing::ExitCodes indexer_exit = indexPepAndProtIds_(fasta_db, light_protein_identifications, light_peptide_identifications);
      if (indexer_exit != PeptideIndexing::EXECUTION_OK)
      {
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unable to execute");
      }
      // call calculatePEP_
      calculatePEP_(light_protein_identifications, light_peptide_identifications);

      /////////////////////////////////////////////////
      /////             heavy                    /////
      ///////////////////////////////////////////////
      vector<ProteinIdentification> heavy_protein_identifications;
      vector<PeptideIdentification> heavy_peptide_identifications;
      // read the identification file for heavy (contains both protein as well as peptide identifications)
      LOG_INFO << "Loading ID file (heavy): " << heavy[i] << endl;
      idXML_file.load(heavy[i], heavy_protein_identifications,  heavy_peptide_identifications);
      indexer_exit = indexPepAndProtIds_(fasta_db, heavy_protein_identifications, heavy_peptide_identifications);

      if (indexer_exit != PeptideIndexing::EXECUTION_OK)
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unable to execute");
      }
      // call calculatePEP_
      calculatePEP_(heavy_protein_identifications, heavy_peptide_identifications);

      //merge light and heavy (for PROTEIN IDs) and stores the merged files in a vector
      vector<ProteinIdentification>  mergedProtIDs;
      mergedProtIDs = mergeProteinIDs_(light_protein_identifications,heavy_protein_identifications);

      //merge light and heavy (for PEPTIDE IDs) and stores the merged files in a vector
      vector<PeptideIdentification>  mergedPeptIDs;
      mergedPeptIDs = mergePeptideIDs_(light_peptide_identifications,heavy_peptide_identifications);

      //create a pair of the merged files
      std::pair <vector<ProteinIdentification>, vector<PeptideIdentification>> mergedProtPepIDs;
      mergedProtPepIDs.first = mergedProtIDs;
      mergedProtPepIDs.second = mergedPeptIDs;

      //store the pairs in the vector, of pair of vectors, ret
      ret.push_back(mergedProtPepIDs);
    }
    return ret;
  };


protected:
  vector<PeptideIdentification> mergePeptideIDs_(
    const vector<PeptideIdentification>& light,
    const vector<PeptideIdentification>& heavy
  )
  /**
   * [mergePeptideIDs_ merge two vectors with peptide identifications]
   * @param  light
   * @param  heavy
   * @return       [vector of PeptideIdentification]
   */
  {
    std::map<String, PeptideIdentification> spectrum_to_id;

    /////////////////////////////////////////////////
    /////             light                    /////
    ///////////////////////////////////////////////
    for (PeptideIdentification const & l : light)
    {
      // insert all the spectrum references with the value l for the light ids in a map
      spectrum_to_id[l.getMetaValue("spectrum_reference")] = l;
    }

    /////////////////////////////////////////////////
    /////             heavy                    /////
    ///////////////////////////////////////////////
    for (PeptideIdentification const & h : heavy)
    {
      const String& h_ref = h.getMetaValue("spectrum_reference"); // get spectrum references

      // if the same id is not found, store heavy id result in map
      if (spectrum_to_id.find(h_ref) == spectrum_to_id.end())
      {
        spectrum_to_id[h_ref] = h;
      }
      else // already a light id in map? replace if score of heavy is better
      {
        const vector<PeptideHit> & l_hits = spectrum_to_id[h_ref].getHits();
        const vector<PeptideHit> & h_hits = h.getHits();

        bool higher_better = h.isHigherScoreBetter();

        // replace if first hit of h has a better score
        if ((l_hits[0].getScore() < h_hits[0].getScore() && higher_better)
        || (l_hits[0].getScore() > h_hits[0].getScore() && !higher_better))
        {
          spectrum_to_id[h_ref] = h;
        }
      }

    }
    // copy values from map into vector
    std::vector<PeptideIdentification> ret;
    for (auto spectrum_id_pair : spectrum_to_id) //range-based for loop for all items in map
    {
      ret.push_back(spectrum_id_pair.second);
    }

    return ret;
  };


protected:
  vector<ProteinIdentification> mergeProteinIDs_(
    const vector<ProteinIdentification>& light,
    const vector<ProteinIdentification>& heavy
  )
  /**
   * [mergeProteinIDs merge two vectors (light and heavy) of protein identifications]
   * @param  light
   * @param  heavy
   * @return       [vector of ProteinIdentification]
   */
  {
    vector<ProteinIdentification> ret;
    //Concatenate two vectors (light and heavy)
    ret.insert(ret.end(), light.begin(), light.end());
    ret.insert(ret.end(), heavy.begin(), heavy.end());

    return ret;
  };


protected:
   void calculateFDR__(vector<PeptideIdentification>& peptide_ids)
   /**
    * [calculateFDR__ estimates false discovery rate of peptide]
    * @param peptide_ids
    */
   {
     FalseDiscoveryRate fdr;
     fdr.apply(peptide_ids);
     IDFilter::filterHitsByScore(peptide_ids, 0.01);
   }


protected:
  void peptFDR_(ProtsPepsPairs& files)
  /**
   * [peptFDR_ from vector of pairs of vectors take peptide vector and compute FDR]
   * @param files
   */
  {
    for (ProtsPepsPairs::iterator j = files.begin(); j != files.end(); j++) //for each pair
    {
      // pass vector of PeptideIdentification to FDR calculation
      calculateFDR__(j->second); // j->second takes the second element of each pair
    }
  }


  /// the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in = getStringList_("in");  // read the mzML filenames
    String database(getStringOption_("database"));  // read the database filename

    StringList in_ids_light = getStringList_("in_ids_light"); //read idXML light file
    StringList in_ids_heavy = getStringList_("in_ids_heavy"); //read idXML heavy file

    String out = getStringOption_("out"); //output in .consensusXML format

    // load (optional) experimental design file
    String design_file = getStringOption_("design");
    ExperimentalDesign design;
    if (!design_file.empty())
    { 
      design = ExperimentalDesignFile::load(design_file, false);
    }
    else
    {  // default to unfractionated design
      ExperimentalDesign::MSFileSection msfs;
      Size count{1};
      for (String & s : in)
      {
        ExperimentalDesign::MSFileSectionEntry e;
        e.fraction = 1;
        e.fraction_group = count;
        e.label = 1;
        e.path = s;
        e.sample = count;
        msfs.push_back(e);
      }      
      design.setMSFileSection(msfs);
    }

    //-------------------------------------------------------------
    // reading input (read protein database)
    //-------------------------------------------------------------
    StopWatch a;
    a.start();
    vector<FASTAFile::FASTAEntry> fasta_db;
    //Loads the identifications of a fasta file
    FASTAFile fasta_reader;
    fasta_reader.load(database, fasta_db);

    ProtsPepsPairs id_files = prepareIDFiles_(fasta_db, in_ids_light, in_ids_heavy);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    peptFDR_(id_files);
    a.stop();
    LOG_INFO << "Peptide Identification tooks: " << a.getClockTime() << " seconds\n" "CPU time: " << a.getCPUTime() << " seconds\n";
    a.reset();

    //-------------------------------------------------------------
    // writing output (after FDR calculation)
    //-------------------------------------------------------------
    IdXMLFile idXML_file;
    ConsensusMap map;

    a.start();
    for (unsigned i = 0; i < in.size(); i++) // for all file names
    {
      const vector<ProteinIdentification> &prot_ids = id_files[i].first; //take the first element from each pair
      const vector<PeptideIdentification> &pept_ids = id_files[i].second; //take the second element from each pair
      // remove extension .mzML
      String out_filename = File::removeExtension(in[i]);
      // add the extension .idXML
      out_filename = out_filename + "_fdr.idXML";
      LOG_INFO << "Writing to file: " << out_filename << endl;
      // test if a file exists, if yes throw exception.
      //if (File::exists(out_filename) == true)
      //{
        try
          {
            (File::exists(out_filename) == true);
          }
        catch (Exception::UnableToCreateFile)
          {
            writeLog_("Error: Unable to create output file.");
            return CANNOT_WRITE_OUTPUT_FILE;
          }
      //}
      idXML_file.store(out_filename, prot_ids, pept_ids);
    }
    a.stop();
    LOG_INFO << "Extension remove tooks: " << a.getClockTime() << " seconds\n" "CPU time: " << a.getCPUTime() << " seconds\n";
    a.reset();

/*
***************************************************************
Q u a n t i f i c a t i o n  &  M a p p i n g
***************************************************************
*/
    a.start();

    // get parameters for all algorithms
    Param paramq = getParam_().copy("Quantification:", true); //set the parameters of the algorithm
    writeDebug_("Parameters passed to FeatureFinderMultiplex algorithm", paramq, 3);
    FeatureFinderMultiplexAlgorithm ffm;
    ffm.setParameters(paramq);
    ffm.setLogType(this->log_type_);

    Param parami = getParam_().copy("IDMapper:", true); //set the parameters of the algorithm
    writeDebug_("Parameters passed to IDMapper algorithm", parami, 3);
    IDMapper id_mapper;
    id_mapper.setParameters(parami);

    MzMLFile mzML_input;

    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    mzML_input.getOptions().setMSLevels(levels);

    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);
    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pp_param);

    Param pq_param = getParam_().copy("Protein Quantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);

    ConsensusXMLFile cons_file;
    ConsensusMap merged_map;

    for (unsigned i = 0; i < in.size(); i++) //for all file names
    {
      MSExperiment exp;
      mzML_input.load(in[i], exp); //load mzML file from stringList

      //** PeakPickerHiRes ** //
      // only pick if not already picked
      PeakMap ms_centroided;
      pp.pickExperiment(exp, ms_centroided, true);
      exp.clear(true);// free memory of profile PeakMaps


      //** FeatureFinderMultiplex ** //
      ffm.run(ms_centroided, true);  // run feature detection algorithm
      ConsensusMap cons_map = ffm.getConsensusMap();

      //** IDMapper **//
      id_mapper.annotate(cons_map, id_files[i].second, id_files[i].first, true, true);
      // annotate output with data processing info
      addDataProcessing_(cons_map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));

      //** IDConflictResolver **//
      IDConflictResolverAlgorithm::resolve(cons_map);
      IDConflictResolverAlgorithm::resolveBetweenFeatures(cons_map);
      addDataProcessing_(cons_map, getProcessingInfo_(DataProcessing::FILTERING));

      // TODO: MultiplexResolver
      //** MultiplexResolver **//

      // ** FileFilter **//
      // remove unassigned peptide identifications
      cons_map.getUnassignedPeptideIdentifications().clear();
      
      //remove unannotated features
      auto removeUnannotatedFeatures = [&](ConsensusFeature feature) -> bool
      {
        return feature.getPeptideIdentifications().empty();
      };
      auto iterator = remove_if(cons_map.begin(), cons_map.end(), removeUnannotatedFeatures);
      cons_map.erase(iterator, cons_map.end());

      //  (we don need this for storing )
      out = File::removeExtension(in[i]);
      out = out + ".consensusXML"; // add extension .idXML
      LOG_INFO << "Writing to file: " << out << endl;
      cons_file.store(out, cons_map); //store results

      // merge consensus maps
      merged_map.appendColumns(cons_map);
    }

    cons_file.store(out, merged_map);

    a.stop();
    LOG_INFO << "Quantification and Mapping took: " << a.getClockTime() << " seconds\n" "CPU time: " << a.getCPUTime() << " seconds\n";
    a.reset();



    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: Implement inference

    vector<PeptideIdentification> infered_peptides;
    ProteinIdentification infered_protein_groups;

    // currenty no proper inference implemented - just extract from consensusMap
    infered_peptides = merged_map.getUnassignedPeptideIdentifications();
    for (ConsensusMap::Iterator cons_it = merged_map.begin();
          cons_it != merged_map.end(); ++cons_it)
    {
      infered_peptides.insert(infered_peptides.end(),
                      cons_it->getPeptideIdentifications().begin(),
                      cons_it->getPeptideIdentifications().end());
      cons_it->getPeptideIdentifications().clear();
    }
    // only keep unique peptides (for now)
    IDFilter::keepUniquePeptidesPerProtein(infered_peptides);

    // merge protein identifications (for now)
    for (auto & p : merged_map.getProteinIdentifications())
    {
      for (auto & h : p.getHits()) { infered_protein_groups.insertHit(h); }
    }

    //-------------------------------------------------------------
    // Peptide quantification
    //-------------------------------------------------------------
    PeptideAndProteinQuant quantifier;
    quantifier.setParameters(pq_param);
    quantifier.readQuantData(merged_map, design);
    quantifier.quantifyPeptides(infered_peptides);

    //-------------------------------------------------------------
    // Protein quantification
    //-------------------------------------------------------------
    // TODO: ProteinQuantifier on (merged?) consensusXML (with 1% FDR?) + inference ids (unfiltered?)? 
    if (infered_protein_groups.getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(
       __FILE__, 
       __LINE__, 
       OPENMS_PRETTY_FUNCTION, 
       "No information on indistinguishable protein groups found.");
    }

    quantifier.quantifyProteins(infered_protein_groups);

    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------

    // Annotate quants to protein(groups) for easier export in mzTab
    auto const & protein_quants = quantifier.getProteinResults();
    PeptideAndProteinQuant::annotateQuantificationsToProteins(protein_quants, infered_protein_groups);
    vector<ProteinIdentification>& proteins = merged_map.getProteinIdentifications();
    proteins.insert(proteins.begin(), infered_protein_groups); // insert inference information as first protein identification

    // Fill MzTab with meta data and quants annotated in identification data structure
    const bool report_unmapped(true);
    const bool report_unidentified_features(false);
    MzTab m = MzTab::exportConsensusMapToMzTab(merged_map, String("null"), report_unidentified_features, report_unmapped);
    MzTabFile().store(out, m);

    return  EXECUTION_OK;
  }

};

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  SILACWorkflow tool;
  return tool.main(argc, argv);
}
/// @endcond
