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

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_SILACWorkflow

    @brief Template for a new Tool
   This tool can be used for scientific stuff.

    And more scientific applications.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_SILACWOrkflow.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SILACWOrkflow.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNewTool :
  public TOPPBase

{
public:
  TOPPNewTool() :
    TOPPBase("SILACWorklfow", "Template for Tool creation", false)
    {}

protected:
  // register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    //Input mzML file list
    registerInputFileList_("in", "<files>", StringList(), "Input MzMLFile");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_ids_light", "<files>", StringList(), "Input idXML light");
    setValidFormats_("in_ids_light", ListUtils::create<String>("idXML"));

    registerInputFileList_("in_ids_heavy", "<files>", StringList(), "Input idXML heavy");
    setValidFormats_("in_ids_heavy", ListUtils::create<String>("idXML"));

    registerInputFile_("database", "<file>", "", "input file ");
    setValidFormats_("database", ListUtils::create<String>("fasta"));
    //registerInputFile_("accession", "<file>", "","Input IdXML file, containing the identified peptides.", true);
    //setValidFormats_("accession", ListUtils::create<String>("idXML"));

    //registerStringOption_("method", "<choice>", "whitelist", "Switch between white/blacklisting", false);

    //setValidStrings_("method", ListUtils::create<String>("whitelist, blacklist"));

    //registerOutputFile_("out", "<file>", "", "Output FASTA file where the reduced database will be written to.");
    //setValidFormats_("out", ListUtils::create<String>("fasta"));
  }

//***************************************************************
// P e p t i d e  I d e n t i f i c a t i o n
//***************************************************************

 //gives indexes to all protein and peptide ids
public:
  // TODO: fasta datenbank als const & übergeben
  PeptideIndexing::ExitCodes indexPepAndProtIds(
    vector<FASTAFile::FASTAEntry>& fasta_db,
    vector<ProteinIdentification>& protein_ids,
    vector<PeptideIdentification>& peptide_ids
    )
  {
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:specificity", "none");
    param_pi.setValue("missing_decoy_action", "warn");
    param_pi.setValue("log", getStringOption_("log"));
    //set as default enzyme TrypsinP
    param_pi.setValue("enzyme:name", "Trypsin/P");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run( fasta_db,protein_ids, peptide_ids);

    return indexer_exit;
  };


//iterates over idXML files
public:
  void calculatePEP(
    vector<ProteinIdentification>& protein_identifications,
    vector<PeptideIdentification>& peptide_identifications)
  {
    //generate PEP_model
    Param pep_param = getParam_().copy("Posterior Error Probability:", true);
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

//prepares ID Files
public:
  // pair of protein and peptide identifications, as stored in an idXML file
  using ProtsPepsPair = std::pair<vector<ProteinIdentification>, vector<PeptideIdentification>>;

  using ProtsPepsPairs = vector<ProtsPepsPair>;

  ProtsPepsPairs prepareIDFiles(
    vector<FASTAFile::FASTAEntry>& fasta_db,
    const StringList & light,
    const StringList & heavy
    )
  {
     // TODO: checken ob light und heavy Stringlist gleich viele einträge haben
     //if stringList sizes not equal, throw an exception
     if (light.size() != heavy.size())
     {
        throw string("Not equal-sized lists"); //(?????)
     }
     IdXMLFile idXML_file;

     ProtsPepsPairs ret;
     for (unsigned int i = 0; i < light.size(); i++)
     {
      vector<ProteinIdentification> light_protein_identifications;
      vector<PeptideIdentification> light_peptide_identifications;
      // read the identification file for light (contains both protein as well as peptide identifications)
      idXML_file.load(light[i], light_protein_identifications,  light_peptide_identifications);
      PeptideIndexing::ExitCodes indexer_exit = indexPepAndProtIds(fasta_db, light_protein_identifications, light_peptide_identifications);
      if (indexer_exit != PeptideIndexing::EXECUTION_OK)
      {
        // TODO: write error message and return
      }
      // calls calculatePEP
      calculatePEP(light_protein_identifications, light_peptide_identifications);

      //the same for heavy
      vector<ProteinIdentification> heavy_protein_identifications;
      vector<PeptideIdentification> heavy_peptide_identifications;
      // read the identification file for heavy (contains both protein as well as peptide identifications)
      idXML_file.load(heavy[i], heavy_protein_identifications,  heavy_peptide_identifications);
      indexer_exit = indexPepAndProtIds( fasta_db,heavy_protein_identifications, heavy_peptide_identifications);

      if (indexer_exit != PeptideIndexing::EXECUTION_OK)
      {
        // TODO: write error message and return
      }

      // calls calculatePEP
      calculatePEP(heavy_protein_identifications, heavy_peptide_identifications);

      vector<ProteinIdentification>  mergedProtIDs;
      mergedProtIDs = mergeProteinIDs(light_protein_identifications,heavy_protein_identifications);
      //mergedProtIDs.push_back(mergeProteinIDs(light_protein_identifications,heavy_protein_identifications));
      //oder mergedProteinIDs.push_back(mergeProtIDs);
      vector<PeptideIdentification>  mergedPeptIDs;
      mergedPeptIDs = mergePeptideIDs(light_peptide_identifications,heavy_peptide_identifications);
      //mergedPeptIDs.push_back(mergePeptideIDs(light_peptide_identifications,heavy_peptide_identifications));
      ////oder mergedPeptideIDs.push_back(mergePeptIDs);
      std::pair <vector<ProteinIdentification>, vector<PeptideIdentification>> mergedProtPepIDs;
      mergedProtPepIDs.first = mergedProtIDs;
      mergedProtPepIDs.second = mergedPeptIDs;
      ret.push_back(mergedProtPepIDs);
     }
    return ret;
  };

  //merges two vectors with peptide identifications
public:
  vector<PeptideIdentification> mergePeptideIDs(
    const vector<PeptideIdentification>& light,
    const vector<PeptideIdentification>& heavy
  )
  {
    std::map<String, PeptideIdentification> spectrum_to_id;

    // inserts all the spectrum references with the value l for the light ids in a map
    for (PeptideIdentification const & l : light)
    {
      spectrum_to_id[l.getMetaValue("spectrum_reference")] = l;
    }

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


//merges two vectors with protein identifications
public:
  vector<ProteinIdentification> mergeProteinIDs(
    const vector<ProteinIdentification>& light,
    const vector<ProteinIdentification>& heavy)
  {
    vector<ProteinIdentification> ret;
    //Concatenates two vectors (light and heavy)
    ret.insert(ret.end(), light.begin(), light.end());
    ret.insert(ret.end(), heavy.begin(), heavy.end());

    return ret;
  };


//estimates false discovery rate of peptide
 private:
   void calculateFDR_(vector<PeptideIdentification>& peptide_ids)
   {
    FalseDiscoveryRate fdr;
    Param p = fdr.getDefaults();
    p.setValue("FDR:PSM", 0.01);
    fdr.setParameters(p);
    fdr.apply(peptide_ids);
   }


//*************************************************************
// Q u a n t i f i c a t i o n
//*************************************************************





// the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in = getStringList_("in");;  // read the mzML filenames
    String database(getStringOption_("database"));  // read the database filename
    StringList in_ids_heavy = getStringList_("in_ids_heavy"); //read idXML heavy file
    StringList in_ids_light = getStringList_("in_ids_light"); //read idXML light file

    //-------------------------------------------------------------
    // reading input (read protein database)
    //-------------------------------------------------------------
    vector<FASTAFile::FASTAEntry> fasta_db;

   //Loads the identifications of a fasta file
    FASTAFile fasta_reader;
    fasta_reader.load(database, fasta_db);

    ProtsPepsPairs id_files = prepareIDFiles(fasta_db, in_ids_light, in_ids_heavy);

    // For FIDO Adapter: merge all
/*
    //-------------------------------------------------------------
    // calculations
   //-------------------------------------------------------------


     // iteration in idXML- heavy file
     for(i = 0; i < in_ids_heavy.size(); i++)
     {
          // read the identification file (contains both protein as well as peptide identifications)
          vector<protein_identification> protein_identifications;
          vector<PeptideIdentification> peptide_identifications;

          //Loads the identifications of an idXML file
          IdXMLFile idXML_file;
          idXML_file.load(in_ids_heavy[i], protein_identifications, peptide_identifications);

          //PeptideIndexing.run
          PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

          if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
              (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
          {
            if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
            {
              return INPUT_FILE_EMPTY;
            }
            else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
            {
              return UNEXPECTED_RESULT;
            }
            else
            {
              return UNKNOWN_ERROR;
            }
          }
}*/

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    //this function compares set of proteins to fasta_accession
    /** vector<FASTAFile::FASTAEntry> db_new;

    for (Size i = 0; i != db.size() ; ++i)
    {
       const String& fasta_accession = db[i].identifier;
       const bool found = id_accessions.find(fasta_accession) != id_accessions.end();
       if ( (found && whitelist) || (!found && !whitelist) )
       {
          db_new.push_back(db[i]);
       }
    }
*/
    return  EXECUTION_OK;
  }

};

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPNewTool tool; // TODO: change name to SILACWorklfow
  return tool.main(argc, argv);
}
/// @endcond
