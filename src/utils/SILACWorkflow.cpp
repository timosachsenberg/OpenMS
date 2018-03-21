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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <vector>

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
// this function will be used to register the tool parameters
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


public:
  PeptideIndexing::ExitCodes indexPepAndProtIds(vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids)
  {
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:specificity", "none");
    param_pi.setValue("missing_decoy_action", "warn");
    param_pi.setValue("log", getStringOption_("log"));
    indexer.setParameters(param_pi);

    vector<FASTAFile::FASTAEntry> fasta_db;

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    return indexer_exit;
  }



//iterates over idXML files
public:
  void indexFiles(const StringList & inp)
  {
     for(unsigned int i = 0; i < inp.size(); i++)
     {
      // read the identification file (contains both protein as well as peptide identifications)
      vector<ProteinIdentification> protein_identifications;
      vector<PeptideIdentification> peptide_identifications;

      //Loads the identifications of an idXML file named idXML_file
      IdXMLFile idXML_file;
      idXML_file.load(inp[i], protein_identifications,  peptide_identifications);

      PeptideIndexing::ExitCodes indexer_exit = indexPepAndProtIds(protein_identifications, peptide_identifications);

      if (indexer_exit != PeptideIndexing::EXECUTION_OK)
      {
        // TODO:
      }
    }
}

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
    // reading input
    //-------------------------------------------------------------
    // read the protein database
    vector<FASTAFile::FASTAEntry> fasta_db;

   //Loads the identifications of a fasta file
    FASTAFile fasta_reader;
    fasta_reader.load(database, fasta_db);

    indexFiles(in_ids_heavy);
    indexFiles(in_ids_light);

    //-------------------------------------------------------------
    // calculations
/**    //-------------------------------------------------------------
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("decoy_string_position", "prefix");
      param_pi.setValue("enzyme:specificity", "none");
      param_pi.setValue("missing_decoy_action", "warn");
      param_pi.setValue("log", getStringOption_("log"));
      indexer.setParameters(param_pi);



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
  TOPPNewTool tool;
  return tool.main(argc, argv);
}
/// @endcond
