// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <iostream>
#include <iomanip>
#include <ostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MassCalculator MassCalculator

    @brief Calculates masses and mass-to-charge ratios of peptide sequences.

    Given a peptide sequence and a charge state, the charged mass (including H+ adducts) and the mass-to-charge ratio are computed.
  The peptide sequence can include modifications (for information on valid notation see the @ref OpenMS::AASequence "AASequence" class documentation).
  Neutral masses can be computed by using "0" as charge state.

    Input can be given directly as values of the parameters: @p in_seq for peptide sequences and @p charge for charge states.
  Alternatively, it can be read from a file (see parameter @p in) with the following format: A peptide sequence at the beginning of each line, optionally followed by any number of charge states.
  Whitespace, commas or semicolons can de used to delimit the different items. Parts of the input that cannot be understood will be skipped.
  If charge states are given in the input file as well as via the @p charge parameter, results are returned for the union of both sets of charge states.

    Output can be written to a file or to the screen (see parameter @p out). Results for different charge states are always ordered from lowest to highest charge.
A number of different output formats are available via the parameter @p format:
    - @p list writes a human-readable list of the form "ABCDEF: z=1 m=566.192 m/z=566.192, z=2 m=567.199 m/z=283.599";
    - @p table produces a CSV-like table (using parameter @p separator to delimit fields) with the columns "peptide", "charge", "mass", and "mass-to-charge", and with one row per peptide and charge state;
    - @p mass_only writes only mass values (one line per peptide, values for different charge states separated by spaces);
    - @p mz_only writes only mass-to-charge ratios (one line per peptide, values for different charge states separated by spaces).


    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MassCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MassCalculator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMassCalculator :
  public TOPPBase
{
public:

  TOPPMassCalculator() :
    TOPPBase("MassCalculator", "Calculates masses and mass-to-charge ratios of peptide sequences", false), use_avg_mass_(false), output_(0), format_(), res_type_(Residue::Full)
  {
    for (Size i = 0; i < Residue::SizeOfResidueType; i++)
    {
      Residue::ResidueType res_type = Residue::ResidueType(i);
      res_type_names_[Residue::getResidueTypeName(res_type)] = res_type;
    }
  }

protected:

  bool use_avg_mass_;
  ostream* output_;  // pointer to output stream (stdout or file)
  String format_, separator_;
  Residue::ResidueType res_type_;
  map<String, Residue::ResidueType> res_type_names_;

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file with peptide sequences and optionally charge numbers (mutually exclusive to 'in_seq')", false);
    setValidFormats_("in",ListUtils::create<String>("txt,fasta"));

    registerStringList_("in_seq", "<peptide_sequences>", StringList(), "List of peptide sequences (mutually exclusive to 'in')", false, false);

    registerOutputFile_("out", "<file>", "", "Output file; if empty, output is written to the screen", false);
    setValidFormats_("out",ListUtils::create<String>("txt"));

    registerIntList_("charge", "<numbers>", ListUtils::create<Int>("0"), "List of charge states; required if 'in_seq' is given", false);
    registerStringOption_("format", "<choice>", "list", "Output format ('list': human-readable list, 'table': CSV-like table, 'mass_only': mass values only, 'mz_only': m/z values only)\n", false);
    setValidStrings_("format", ListUtils::create<String>("list,table,mass_only,mz_only"));
    registerFlag_("average_mass", "Compute average (instead of monoisotopic) peptide masses");
    registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
    setValidStrings_("fragment_type", ListUtils::create<String>("full,internal,N-terminal,C-terminal,a-ion,b-ion,c-ion,x-ion,y-ion,z-ion"));
    registerStringOption_("separator", "<sep>", "", "Field separator for 'table' output format; by default, the 'tab' character is used", false);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 1, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("enzyme", all_enzymes);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:max_size", "<num>", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 0, "Number of missed cleavages.", false, false);
    registerStringOption_("peptide:filter_duplicates", "<toggle>", "true", "Peptides with identical amino acid sequence are counted only once ('true') or multiple ('false') times.", false, false);
    setValidStrings_("peptide:filter_duplicates", ListUtils::create<String>("true,false"));

  }

  double computeMass_(const AASequence& seq, Int charge) const
  {
    if (use_avg_mass_) return seq.getAverageWeight(res_type_, charge);
    else return seq.getMonoWeight(res_type_, charge);
  }

  void writeTable_(const AASequence& seq, const set<Int>& charges)
  {
    SVOutStream sv_out(*output_, separator_);
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      sv_out << seq.toString() << *it << mass;
      sv_out.writeValueOrNan(mass / *it);
      sv_out << endl;
    }
  }

  void writeList_(const AASequence& seq, const set<Int>& charges)
  {
    *output_ << seq.toString() << ": ";
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      if (it != charges.begin()) *output_ << ", ";
      *output_ << "z=" << *it << " m=" << mass << " m/z=";
      if (*it != 0) *output_ << (mass / *it);
      else *output_ << "inf";
    }
    *output_ << endl;
  }

  void writeMassOnly_(const AASequence& seq, const set<Int>& charges,
                      bool mz = false)
  {
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      if (it != charges.begin()) *output_ << " ";
      if (!mz) *output_ << mass;
      else if (*it == 0) *output_ << "inf";
      else *output_ << mass / *it;
    }
    *output_ << endl;
  }

  void writeLine_(const AASequence& seq, const set<Int>& charges)
  {
    if (format_ == "list") writeList_(seq, charges);
    else if (format_ == "table") writeTable_(seq, charges);
    else if (format_ == "mass_only") writeMassOnly_(seq, charges);
    else writeMassOnly_(seq, charges, true); // "mz_only"
  }

  String getItem_(String& line, const String& skip = " \t,;")
  {
    Size pos = line.find_first_of(skip);
    String prefix = line.substr(0, pos);
    pos = line.find_first_not_of(skip, pos);
    if (pos == String::npos) line = "";
    else line = line.substr(pos);
    return prefix;
  }

  void readFile_(const String& filename, const set<Int>& charges)
  {
    ifstream input(filename.c_str());
    String line;
    Size line_count(0);
    while (getline(input, line))
    {
      ++line_count;
      String item = getItem_(line);
      if ((item[0] == '"') && (item[item.size() - 1] == '"'))
      {
        item.unquote();
      }

      AASequence seq;
      try
      {
        seq = AASequence::fromString(item);
      } 
      catch (Exception::ParseError& /*e*/)
      {
        LOG_WARN << "Warning: '" << item << "' is not a valid peptide sequence - skipping\n";
        continue;
      }

      set<Int> local_charges(charges);
      Size conversion_failed_count(0);
      while (!line.empty())
      {
        item = getItem_(line);
        try
        {
          local_charges.insert(item.toInt());
        }
        catch (Exception::ConversionError& /*e*/)
        {
          ++conversion_failed_count;
        }
      }
      if (conversion_failed_count)
      {
        LOG_WARN << "Warning: Invalid charge state specified in line:" << line_count << ".\n";
      }
      if (local_charges.empty())
      {
        LOG_WARN << "Warning: No charge state specified - skipping (line:" << line_count << ")\n";
        continue;
      }
      writeLine_(seq, local_charges);
    }
    input.close();
  }

  ExitCodes main_(int, const char**)
  {
    String in = getStringOption_("in");
    StringList in_seq = getStringList_("in_seq");
    String out = getStringOption_("out");
    IntList charge_list = getIntList_("charge");
    set<Int> charges(charge_list.begin(), charge_list.end());
    use_avg_mass_ = getFlag_("average_mass");
    res_type_ = res_type_names_[getStringOption_("fragment_type")];

    ofstream outfile;
    if (out.empty())
    {
      output_ = &cout;
    }
    else
    {
      outfile.open(out.c_str());
      output_ = &outfile;
    }

    format_ = getStringOption_("format");
    if (format_ == "table")
    {
      separator_ = getStringOption_("separator");
      if (separator_.empty()) separator_ = "\t";
      // write header:
      SVOutStream sv_out(*output_, separator_);
      sv_out << "peptide" << "charge" << "mass" << "mass-to-charge" << endl;
    }

    if ((in.size() > 0) && (in_seq.size() > 0))
    {
      LOG_ERROR << "Specifying an input file and input sequences at the same time is not allowed!";
      return ILLEGAL_PARAMETERS;
    }

    if (in.hasSuffix(".fasta"))
    {
      FASTAFile fastaFile;
      vector<FASTAFile::FASTAEntry> fasta_db;
      fastaFile.load(in, fasta_db);

      const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
      EnzymaticDigestion digestor;
      digestor.setEnzyme(getStringOption_("enzyme"));
      digestor.setMissedCleavages(missed_cleavages);

      StringList fixedModNames = getStringList_("modifications:fixed");
      set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

      if (fixed_unique.size() != fixedModNames.size())
      {
        cout << "duplicate fixed modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      StringList varModNames = getStringList_("modifications:variable");
      set<String> var_unique(varModNames.begin(), varModNames.end());
      if (var_unique.size() != varModNames.size())
      {
        cout << "duplicate variable modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      vector<ResidueModification> fixedMods = getModifications_(fixedModNames);
      vector<ResidueModification> varMods = getModifications_(varModNames);
      Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

      // lookup for processed peptides. must be defined outside of omp section and synchronized
      set<StringView> processed_petides;

      // set minimum / maximum size of peptide after digestion
      Size min_peptide_length = getIntOption_("peptide:min_size");
      Size max_peptide_length = getIntOption_("peptide:max_size");

      map<double, Size> masses;

      bool skip_duplicates = getStringOption_("peptide:filter_duplicates") == "true" ? true : false;

      for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
      {
        vector<StringView> current_digest;
        digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length, max_peptide_length);
        for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
        {
          if (cit->getString().has('X')) continue;

          // if we want to skip duplicates we need to take track of already processes sequences.
          if (skip_duplicates)
          {
            bool already_processed = false;
            if (processed_petides.find(*cit) != processed_petides.end())
            {
              // peptide (and all modified variants) already processed so skip it
              already_processed = true;
            }
 
            if (already_processed)
            {
              continue;
            }
        
            processed_petides.insert(*cit);
          }
      
          vector<AASequence> all_modified_peptides;

          AASequence aas = AASequence::fromString(cit->getString());
          ModifiedPeptideGenerator::applyFixedModifications(fixedMods.begin(), fixedMods.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(varMods.begin(), varMods.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
          
          for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
          {
            const AASequence& candidate = all_modified_peptides[mod_pep_idx];
            double current_peptide_mass = candidate.getMonoWeight();

            // round to 6 decimal digits
            current_peptide_mass = round( current_peptide_mass * 1e6 ) / 1e6; 
            masses[current_peptide_mass]++;  
          }
        }
      }

      // write header:
      outfile << "mass\tcount" << fixed << setprecision(6) << endl;
      for (map<double, Size>::const_iterator mit = masses.begin(); mit != masses.end(); ++mit)
      {
        outfile << mit->first << "\t" << mit->second << endl;
      }    
    } else if (in.size() > 0)
    {
      readFile_(in, charges);
    }
    else
    {
      if (charges.empty())
      {
        LOG_ERROR << "Error: No charge state specified";
        return ILLEGAL_PARAMETERS;
      }
      for (StringList::iterator it = in_seq.begin(); it != in_seq.end(); ++it)
      {
        AASequence seq;
        try
        {
         seq = AASequence::fromString(*it);
        }
        catch (Exception::ParseError& /*e*/)
        {
          LOG_WARN << "Warning: '" << *it << "' is not a valid peptide sequence - skipping\n";
          continue;
        }

        writeLine_(seq, charges);
      }
    }

    if (!out.empty()) outfile.close();

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMassCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
