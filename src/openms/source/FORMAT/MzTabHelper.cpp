// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Timo Sachsenberg, Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzTabHelper.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <vector>
#include <set>

using namespace std;

namespace OpenMS
{

  MzTabHelper::MzTabHelper()
  {

  }

  MzTabHelper::~MzTabHelper()
  {

  }
  
  void MzTabHelper::addPepEvidenceToRows(const vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows)
  {
    if (!peptide_evidences.empty())
    {
      for (Size i = 0; i != peptide_evidences.size(); ++i)
      {
        // get AABefore and AAAfter as well as start and end for all pep evidences

        // pre/post
        // from spec: Amino acid preceding the peptide (coming from the PSM) in the protein
        // sequence. If unknown “null” MUST be used, if the peptide is N-terminal “-“
        // MUST be used.
        if (peptide_evidences[i].getAABefore() == PeptideEvidence::UNKNOWN_AA)
        {
          row.pre = MzTabString("null");
        }
        else if (peptide_evidences[i].getAABefore() == PeptideEvidence::N_TERMINAL_AA)
        {
          row.pre = MzTabString("-");
        }
        else
        {
          row.pre = MzTabString(String(peptide_evidences[i].getAABefore()));
        }

        if (peptide_evidences[i].getAAAfter() == PeptideEvidence::UNKNOWN_AA)
        {
          row.post = MzTabString("null");
        }
        else if (peptide_evidences[i].getAAAfter() == PeptideEvidence::C_TERMINAL_AA)
        {
          row.post = MzTabString("-");
        }
        else
        {
          row.post = MzTabString(String(peptide_evidences[i].getAAAfter()));
        }

        // start/end
        if (peptide_evidences[i].getStart() == PeptideEvidence::UNKNOWN_POSITION)
        {
          row.start = MzTabString("null");
        }
        else
        {
          row.start = MzTabString(String(peptide_evidences[i].getStart() + 1)); // counting in mzTab starts at 1
        }

        if (peptide_evidences[i].getEnd() == PeptideEvidence::UNKNOWN_POSITION)
        {
          row.end = MzTabString("null");
        }
        else
        {
          row.end = MzTabString(String(peptide_evidences[i].getEnd() + 1)); // counting in mzTab starts at 1
        }

        row.accession = MzTabString(peptide_evidences[i].getProteinAccession());

        rows.push_back(row);
      }
    }
    else
    { // report without pep evidence information
      row.pre = MzTabString("null");
      row.post = MzTabString("null");
      row.start = MzTabString("null");
      row.end = MzTabString("null");
      rows.push_back(row);
    }
  }
  
  void MzTabHelper::addMetaInfoToOptionalColumns(const set<String>& keys, vector<MzTabOptionalColumnEntry>& opt, const String id, const MetaInfoInterface meta)
  {
    for (set<String>::const_iterator sit = keys.begin(); sit != keys.end(); ++sit)
    {
      const String& key = *sit;
      MzTabOptionalColumnEntry opt_entry;
      opt_entry.first = String("opt_") + id + String("_") + String(key).substitute(' ','_');
      if (meta.metaValueExists(key))
      {
        opt_entry.second = MzTabString(meta.getMetaValue(key).toString().substitute(' ','_'));
      } // otherwise it is default ("null")
      opt.push_back(opt_entry);
    }
  }

  map<Size, MzTabModificationMetaData> MzTabHelper::generateMzTabStringFromModifications(const vector<String>& mods)
  {
    map<Size, MzTabModificationMetaData> mods_mztab;
    for (vector<String>::const_iterator sit = mods.begin(); sit != mods.end(); ++sit)
    {
      Size index = (sit - mods.begin()) + 1;
      MzTabModificationMetaData mod;
      MzTabParameter mp;
      mp.setCVLabel("UNIMOD");
      ModificationsDB* mod_db = ModificationsDB::getInstance();
      // MzTab standard is to just report Unimod accession.
      ResidueModification m = mod_db->getModification(*sit);
      String unimod_accession = m.getUniModAccession();
      mp.setAccession(unimod_accession.toUpper());
      mp.setName(m.getId());
      mod.modification = mp;

      if (m.getTermSpecificity() == ResidueModification::C_TERM)
      {
        mod.position = MzTabString("Any C-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::N_TERM)
      {
        mod.position = MzTabString("Any N-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::ANYWHERE)
      {
        mod.position = MzTabString("Anywhere");
      }
      else if (m.getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        mod.position = MzTabString("Protein C-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        mod.position = MzTabString("Protein N-term");
      }

      mod.site = MzTabString(m.getOrigin());
      mods_mztab[index] = mod;
    }
    return mods_mztab;
  }

  map<Size, MzTabModificationMetaData> MzTabHelper::generateMzTabStringFromVariableModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002454, No variable modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  map<Size, MzTabModificationMetaData> MzTabHelper::generateMzTabStringFromFixedModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002453, No fixed modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  // Generate MzTab style list of PTMs from AASequence object. 
  // All passed fixed modifications are not reported (as suggested by the standard for the PRT and PEP section).
  // In contrast, all modifications are reported in the PSM section (see standard document for details).
  MzTabModificationList MzTabHelper::extractModificationListFromAASequence(const AASequence& aas, const vector<String>& fixed_mods = vector<String>())
  {
    MzTabModificationList mod_list;
    vector<MzTabModification> mods;
  
    if (aas.isModified())
    {
      if (aas.hasNTerminalModification())
      {
        MzTabModification mod;
        const ResidueModification& res_mod = *(aas.getNTerminalModification());
        if (find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
        {
          MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
          vector<pair<Size, MzTabParameter> > pos;
          pos.push_back(make_pair(0, MzTabParameter()));
          mod.setModificationIdentifier(unimod_accession);
          mod.setPositionsAndParameters(pos);
          mods.push_back(mod);
        }
      }
  
      for (Size ai = 0; ai != aas.size(); ++ai)
      {
        if (aas[ai].isModified())
        {
          MzTabModification mod;
          const ResidueModification& res_mod = *(aas[ai].getModification());
          if (find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
          {
            // MzTab standard is to just report Unimod accession.
            MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
            vector<pair<Size, MzTabParameter> > pos;
            pos.push_back(make_pair(ai + 1, MzTabParameter()));
            mod.setPositionsAndParameters(pos);
            mod.setModificationIdentifier(unimod_accession);
            mods.push_back(mod);
          }
        }
      }
  
  if (aas.hasCTerminalModification())
  {
  MzTabModification mod;
  const ResidueModification& res_mod = *(aas.getCTerminalModification());
  if (find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
  {
  MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
  vector<pair<Size, MzTabParameter> > pos;
  pos.push_back(make_pair(aas.size() + 1, MzTabParameter()));
  mod.setPositionsAndParameters(pos);
  mod.setModificationIdentifier(unimod_accession);
  mods.push_back(mod);
  }
  }
  }
  mod_list.set(mods);
  return mod_list;
  }
  

  MzTab MzTabHelper::exportFeatureMapToMzTab(const FeatureMap& feature_map, const String& filename)
  {
    if (!filename.empty())
    {
      LOG_INFO << "exporting feature map: \"" << filename << "\" to mzTab: " << endl;
    }
    MzTab mztab;
    MzTabMetaData meta_data;

    vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        
    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;
    if (!prot_ids.empty())
    {
      ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
      var_mods = sp.variable_modifications;
      fixed_mods = sp.fixed_modifications;
      db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
    }

    meta_data.variable_mod = generateMzTabStringFromVariableModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromFixedModifications(fixed_mods);

    // mandatory meta values
    meta_data.mz_tab_type = MzTabString("Quantification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("Export from featureXML");

    MzTabMSRunMetaData ms_run;
    ms_run.location = feature_map.getPrimaryMSRunPath().empty() ? MzTabString("null") : MzTabString(feature_map.getPrimaryMSRunPath()[0]);
    meta_data.ms_run[1] = ms_run;
    meta_data.uri[1] = MzTabString(filename);
    meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO: we currently only support psm search engine scores annotated to the identification run
    meta_data.peptide_search_engine_score[1] = MzTabParameter();

    mztab.setMetaData(meta_data);

    // pre-analyze data for occuring meta values at feature and peptide hit level
    // these are used to build optional columns containing the meta values in internal data structures

    set<String> feature_user_value_keys;
    set<String> peptide_hit_user_value_keys;
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      const Feature& f = feature_map[i];
      vector<String> keys;
      f.getKeys(keys); //TODO: why not just return it?
      feature_user_value_keys.insert(keys.begin(), keys.end());

      const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
      for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          vector<String> ph_keys;
          hit->getKeys(ph_keys);
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

    MzTabPeptideSectionRows rows;
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      MzTabPeptideSectionRow row;
      const Feature& f = feature_map[i];
      row.mass_to_charge = MzTabDouble(f.getMZ());
      MzTabDoubleList rt_list;
      vector<MzTabDouble> rts;
      rts.push_back(MzTabDouble(f.getRT()));
      rt_list.set(rts);
      row.retention_time = rt_list;

      // set rt window if a bounding box has been set
      vector<MzTabDouble> window;
      if (f.getConvexHull().getBoundingBox() != DBoundingBox<2>())
      {
        window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().minX()));
        window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().maxX()));
      }

      MzTabDoubleList rt_window;
      rt_window.set(window);
      row.retention_time_window = rt_window;
      row.charge = MzTabInteger(f.getCharge());
      row.peptide_abundance_stdev_study_variable[1];
      row.peptide_abundance_std_error_study_variable[1];
      row.peptide_abundance_study_variable[1] = MzTabDouble(f.getIntensity());
      row.best_search_engine_score[1] = MzTabDouble();
      row.search_engine_score_ms_run[1][1] = MzTabDouble();

      // create opt_ column for peptide sequence containing modification
      MzTabOptionalColumnEntry opt_global_modified_sequence;
      opt_global_modified_sequence.first = String("opt_global_modified_sequence");
      row.opt_.push_back(opt_global_modified_sequence);

      // create and fill opt_ columns for feature (peptide) user values
      addMetaInfoToOptionalColumns(feature_user_value_keys, row.opt_, String("global"), f);

      vector<PeptideIdentification> pep_ids = f.getPeptideIdentifications();
      if (pep_ids.empty())
      {
        rows.push_back(row);
        continue;
      }

      // TODO: here we assume that all have the same score type etc.
      vector<PeptideHit> all_hits;
      for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
      }

      if (all_hits.empty())
      {
        rows.push_back(row);
        continue;
      }

      // create new peptide id object to assist in sorting
      PeptideIdentification new_pep_id = pep_ids[0];
      new_pep_id.setHits(all_hits);
      new_pep_id.assignRanks();

      const PeptideHit& best_ph = new_pep_id.getHits()[0];
      const AASequence& aas = best_ph.getSequence();
      row.sequence = MzTabString(aas.toUnmodifiedString());

      row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

      const set<String>& accessions = best_ph.extractProteinAccessions();
      const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

      row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
      // select accession of first peptide_evidence as representative ("leading") accession
      row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());
      row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());
      row.search_engine_score_ms_run[1][1] = MzTabDouble(best_ph.getScore());

      // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
      for (Size i = 0; i != row.opt_.size(); ++i)
      {
        MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

        if (opt_entry.first == String("opt_global_modified_sequence"))
        {
          opt_entry.second = MzTabString(aas.toString());
        }
      }

      // create and fill opt_ columns for psm (PeptideHit) user values
      addMetaInfoToOptionalColumns(peptide_hit_user_value_keys, row.opt_, String("global"), best_ph);

      rows.push_back(row);
    }
    mztab.setPeptideSectionRows(rows);
    return mztab;
  }

  MzTab MzTabHelper::exportIdentificationsToMzTab(const vector<ProteinIdentification>& prot_ids, const vector<PeptideIdentification>& peptide_ids, const String& filename)
  {
    if (!filename.empty())
    {
      LOG_INFO << "exporting identifications: \"" << filename << "\" to mzTab: " << endl;
    }
    vector<PeptideIdentification> pep_ids = peptide_ids;
    MzTab mztab;
    MzTabMetaData meta_data;
    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;
    String search_engine;
    String search_engine_version;

    if (!prot_ids.empty())
    {
      search_engine = prot_ids[0].getSearchEngine();
      search_engine_version = prot_ids[0].getSearchEngineVersion();
    }

    if (!prot_ids.empty())
    {
      MzTabParameter protein_score_type;
      protein_score_type.fromCellString("[,,custom score,]"); // TODO at least it should be noted if higher score is better. Better document type of score
      meta_data.protein_search_engine_score[1] = protein_score_type; // TODO add meta value to ProteinIdentification
      ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
      var_mods = sp.variable_modifications;
      fixed_mods = sp.fixed_modifications;
      db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);

      //sp.digestion_enzyme
      //sp.missed_cleavages
      // generate protein section
      MzTabProteinSectionRows protein_rows;

      Size current_run_index(1);
      for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin();
       it != prot_ids.end(); ++it, ++current_run_index)
      {
        const vector<ProteinIdentification::ProteinGroup>& protein_groups = it->getProteinGroups();
        const vector<ProteinIdentification::ProteinGroup>& indist_groups = it->getIndistinguishableProteins();
        const vector<ProteinHit>& protein_hits = it->getHits();

        MzTabMSRunMetaData ms_run;
        ms_run.location = it->getPrimaryMSRunPath().empty() ? MzTabString("null") : MzTabString(it->getPrimaryMSRunPath()[0]);
        // TODO: add processing information that this file has been exported from "filename"
        meta_data.ms_run[current_run_index] = ms_run;

        // pre-analyze data for occuring meta values at protein hit level
        // these are used to build optional columns containing the meta values in internal data structures
        set<String> protein_hit_user_value_keys = 
          MetaInfoInterfaceUtils::findCommonMetaKeys<vector<ProteinHit>, set<String> >(protein_hits.begin(), protein_hits.end(), 100.0);

        // we do not want descriptions twice
        protein_hit_user_value_keys.erase("Description");

        for (Size i = 0; i != protein_hits.size(); ++i)
        {
          const ProteinHit& hit = protein_hits[i];
          MzTabProteinSectionRow protein_row;

          protein_row.accession = MzTabString(hit.getAccession());
          protein_row.description = MzTabString(hit.getDescription()); 
       // protein_row.taxid = hit.getTaxonomyID(); // TODO add as meta value to protein hitNEWT taxonomy for the species.
       // MzTabString species = hit.getSpecies(); // Human readable name of the species
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.
          protein_row.best_search_engine_score[1] = MzTabDouble(hit.getScore());
       // MzTabParameterList search_engine; // Search engine(s) identifying the protein.
       // map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
       // map<Size, map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
       // MzTabInteger reliability;
       // map<Size, MzTabInteger> num_psms_ms_run;
       // map<Size, MzTabInteger> num_peptides_distinct_ms_run;
       // map<Size, MzTabInteger> num_peptides_unique_ms_run;
       // MzTabModificationList modifications; // Modifications identified in the protein.
       // MzTabString uri; // Location of the protein’s source entry.
       // MzTabStringList go_terms; // List of GO terms for the protein.
          double coverage = hit.getCoverage();
          protein_row.protein_coverage = coverage >= 0 ? MzTabDouble(coverage) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
       // vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”

          // create and fill opt_ columns for protein hit user values
          addMetaInfoToOptionalColumns(protein_hit_user_value_keys, protein_row.opt_, String("global"), hit);

          // optional column for protein groups
          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("single_protein");
          protein_row.opt_.push_back(opt_column_entry);
           
          protein_rows.push_back(protein_row);
        }

        // Protein groups are currently simply PRT rows with extra opt columns
        for (Size i = 0; i != protein_groups.size(); ++i)
        {
          const ProteinIdentification::ProteinGroup& group = protein_groups[i];
          MzTabProteinSectionRow protein_row;
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.
            
          MzTabStringList ambiguity_members;
          ambiguity_members.setSeparator(',');
          vector<MzTabString> entries;
          for (Size j = 0; j != group.accessions.size() ; ++j)
          {
            // set accession and description to first element of group
            if (j == 0)
            {
              protein_row.accession = MzTabString(group.accessions[j]);
              // protein_row.description  // TODO: how to set description? information not contained in group
            }
            entries.push_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);

          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("protein_group");
          protein_row.opt_.push_back(opt_column_entry);
          protein_rows.push_back(protein_row);
        }

        for (Size i = 0; i != indist_groups.size(); ++i)
        {
          const ProteinIdentification::ProteinGroup& group = indist_groups[i];
          MzTabProteinSectionRow protein_row;
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.
            
          MzTabStringList ambiguity_members;
          ambiguity_members.setSeparator(',');
          vector<MzTabString> entries;
          for (Size j = 0; j != group.accessions.size() ; ++j)
          {
            // set accession and description to first element of group
            if (j == 0)
            {
              protein_row.accession = MzTabString(group.accessions[j]);
            }
            entries.push_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("indistinguishable_group");
          protein_row.opt_.push_back(opt_column_entry);
          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);
          //vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”
          protein_rows.push_back(protein_row);
        }

      }
      mztab.setProteinSectionRows(protein_rows);
    }
    // end protein groups

    // start PSMs

    // mandatory meta values
    meta_data.mz_tab_type = MzTabString("Identification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("Export from idXML");

    meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);
    MzTabParameter psm_search_engine_score;
    psm_search_engine_score.fromCellString("[,," + search_engine + "," + search_engine_version + "]");
    meta_data.psm_search_engine_score[1] = psm_search_engine_score;

    mztab.setMetaData(meta_data);

    MzTabPSMSectionRows rows;
    Size psm_id(0);
    for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it, ++psm_id)
    {
      // skip empty peptide identification objects
      if (it->getHits().empty())
      {
        continue;
      }

      // sort by rank
      it->assignRanks();

      MzTabPSMSectionRow row;

      // only consider best peptide hit for export
      const PeptideHit& best_ph = it->getHits()[0];
      const AASequence& aas = best_ph.getSequence();
      row.sequence = MzTabString(aas.toUnmodifiedString());

      // extract all modifications in the current sequence for reporting. In contrast to peptide and protein section all modifications are reported.
      row.modifications = extractModificationListFromAASequence(aas);

      row.PSM_ID = MzTabInteger(psm_id);
      row.database = db;
      row.database_version = db_version;
      MzTabParameterList search_engines;
      search_engines.fromCellString("[,," + search_engine + "," + search_engine_version + "]");
      row.search_engine = search_engines;

      row.search_engine_score[1] = MzTabDouble(best_ph.getScore());
      vector<MzTabDouble> rts_vector;
      rts_vector.push_back(MzTabDouble(it->getRT()));
      MzTabDoubleList rts;
      rts.set(rts_vector);
      row.retention_time = rts;
      row.charge = MzTabInteger(best_ph.getCharge());
      row.exp_mass_to_charge = MzTabDouble(it->getMZ());
      row.calc_mass_to_charge = best_ph.getCharge() != 0 ? MzTabDouble(aas.getMonoWeight(Residue::Full, best_ph.getCharge()) / best_ph.getCharge()) : MzTabDouble();

      // add opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
      MzTabOptionalColumnEntry opt_entry;
      opt_entry.first = String("opt_global_modified_sequence");
      opt_entry.second = MzTabString(aas.toString());
      row.opt_.push_back(opt_entry);

      // currently write all keys
      // TODO: percentage procedure with MetaInfoInterfaceUtils
      vector<String> ph_keys;
      best_ph.getKeys(ph_keys);
      // TODO: no conversion but make funtion on collections
      set<String> ph_key_set(ph_keys.begin(), ph_keys.end());
      addMetaInfoToOptionalColumns(ph_key_set, row.opt_, String("global"), best_ph);

      // TODO Think about if the uniqueness can be determined by # of peptide evidences
      // b/c this would only differ when evidences come from different DBs
      const set<String>& accessions = best_ph.extractProteinAccessions();
      row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);

      // create row for every PeptideEvidence entry (mapping to a protein)
      const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

      // pass common row entries and create rows for all peptide evidences
      addPepEvidenceToRows(peptide_evidences, row, rows);
    }

    mztab.setPSMSectionRows(rows);

    return mztab;
  }

  MzTab MzTabHelper::exportConsensusMapToMzTab(const ConsensusMap& consensus_map, const String& filename)
  {
    if (!filename.empty())
    {
      LOG_INFO << "exporting consensus map: \"" << filename << "\" to mzTab: " << endl;
    }
    MzTab mztab;

    // TODO:refactor from IDSplitter
    vector<PeptideIdentification> pep_ids(consensus_map.getUnassignedPeptideIdentifications());

    for (ConsensusMap::ConstIterator cons_it = consensus_map.begin();
         cons_it != consensus_map.end(); ++cons_it)
    {
      pep_ids.insert(pep_ids.end(),
                      cons_it->getPeptideIdentifications().begin(),
                      cons_it->getPeptideIdentifications().end());
    }

    vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();
    mztab = exportIdentificationsToMzTab(prot_ids, pep_ids, "");

    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;
    if (!prot_ids.empty())
    {
      ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
      var_mods = sp.variable_modifications;
      fixed_mods = sp.fixed_modifications;
      db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
    }

    // determine number of channels
    Size n_study_variables = consensus_map.getFileDescriptions().size();

    MzTabMetaData meta_data;

    // mandatory meta values
    meta_data.mz_tab_type = MzTabString("Quantification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("Export from consensusXML");
    meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);
    meta_data.peptide_search_engine_score[1] = MzTabParameter();
    meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO insert search engine information
    MzTabMSRunMetaData ms_run;
    StringList ms_runs = consensus_map.getPrimaryMSRunPath();
    for (Size i = 0; i != ms_runs.size(); ++i)
    {
      ms_run.location = MzTabString(ms_runs[i]);
      meta_data.ms_run[i + 1] = ms_run;
    }

    mztab.setMetaData(meta_data);

    // pre-analyze data for occuring meta values at consensus feature and peptide hit level
    // these are used to build optional columns containing the meta values in internal data structures
    set<String> consensus_feature_user_value_keys;
    set<String> peptide_hit_user_value_keys;
    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      const ConsensusFeature& c = consensus_map[i];
      vector<String> keys;
      c.getKeys(keys);
      consensus_feature_user_value_keys.insert(keys.begin(), keys.end());

      const vector<PeptideIdentification>& pep_ids = c.getPeptideIdentifications();
      for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          vector<String> ph_keys;
          hit->getKeys(ph_keys);
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

    MzTabPeptideSectionRows rows;

    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      MzTabPeptideSectionRow row;
      const ConsensusFeature& c = consensus_map[i];

      // create opt_ column for peptide sequence containing modification
      MzTabOptionalColumnEntry opt_global_modified_sequence;
      opt_global_modified_sequence.first = String("opt_global_modified_sequence");
      row.opt_.push_back(opt_global_modified_sequence);

      // create opt_ columns for consensus feature (peptide) user values
      for (set<String>::const_iterator mit = consensus_feature_user_value_keys.begin(); mit != consensus_feature_user_value_keys.end(); ++mit)
      {
        MzTabOptionalColumnEntry opt_entry;
        const String& key = *mit;
        opt_entry.first = String("opt_global_") + key;
        if (c.metaValueExists(key))
        {
          opt_entry.second = MzTabString(c.getMetaValue(key).toString());
        } // otherwise it is default ("null")
        row.opt_.push_back(opt_entry);
      }

      // create opt_ columns for psm (PeptideHit) user values
      for (set<String>::const_iterator mit = peptide_hit_user_value_keys.begin(); mit != peptide_hit_user_value_keys.end(); ++mit)
      {
        MzTabOptionalColumnEntry opt_entry;
        const String& key = *mit;
        opt_entry.first = String("opt_global_") + key;
        // leave value empty as we have to fill it with the value from the best peptide hit
        row.opt_.push_back(opt_entry);
      }

      row.mass_to_charge = MzTabDouble(c.getMZ());
      MzTabDoubleList rt_list;
      vector<MzTabDouble> rts;
      rts.push_back(MzTabDouble(c.getRT()));
      rt_list.set(rts);
      row.retention_time = rt_list;
      MzTabDoubleList rt_window;
      row.retention_time_window = rt_window;
      row.charge = MzTabInteger(c.getCharge());
      row.best_search_engine_score[1] = MzTabDouble();

      // initialize columns
      for (Size study_variable = 1; study_variable <= n_study_variables; ++study_variable)
      {
        row.peptide_abundance_stdev_study_variable[study_variable] = MzTabDouble();
        row.peptide_abundance_std_error_study_variable[study_variable] = MzTabDouble();
        row.peptide_abundance_study_variable[study_variable] = MzTabDouble();
        row.search_engine_score_ms_run[1][study_variable] = MzTabDouble();
      }

      ConsensusFeature::HandleSetType fs = c.getFeatures();
      for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
      {
        Size study_variable = fit->getMapIndex() + 1;
        row.peptide_abundance_stdev_study_variable[study_variable];
        row.peptide_abundance_std_error_study_variable[study_variable];
        row.peptide_abundance_study_variable[study_variable] = MzTabDouble(fit->getIntensity());
      }

      vector<PeptideIdentification> pep_ids = c.getPeptideIdentifications();
      if (!pep_ids.empty())
      {
        if (pep_ids.size() != 1)
        {
          throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Consensus features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!");
        }

        pep_ids[0].assignRanks();
        const PeptideHit& best_ph = pep_ids[0].getHits()[0];
        const AASequence& aas = best_ph.getSequence();
        row.sequence = MzTabString(aas.toUnmodifiedString());

        row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

        const set<String>& accessions = best_ph.extractProteinAccessions();
        const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

        row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
        // select accession of first peptide_evidence as representative ("leading") accession
        row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());

        row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());

        // fill opt_ columns

        // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
        for (Size i = 0; i != row.opt_.size(); ++i)
        {
          MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

          if (opt_entry.first == String("opt_global_modified_sequence"))
          {
            opt_entry.second = MzTabString(aas.toString());
          }
        }

        // fill opt_ column of psm
        vector<String> ph_keys;
        best_ph.getKeys(ph_keys);
        for (Size k = 0; k != ph_keys.size(); ++k)
        {
          const String& key = ph_keys[k];

          // find matching entry in opt_ (TODO: speed this up)
          for (Size i = 0; i != row.opt_.size(); ++i)
          {
            MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

            if (opt_entry.first == String("opt_global_") + key)
            {
              opt_entry.second = MzTabString(best_ph.getMetaValue(key).toString());
            }
          }
        }
      }

      rows.push_back(row);
    }

    mztab.setPeptideSectionRows(rows);
    return mztab;
  }

  void MzTabHelper::setMSRuns(const StringList& ms_run_locations, MzTab& mztab)
  {
    MzTabMetaData md = mztab.getMetaData();
    Size i(1);
    for (StringList::const_iterator sit = ms_run_locations.begin(); sit != ms_run_locations.end(); ++sit)
    {
      md.ms_run[i].location = MzTabString("file://" + File::absolutePath(*sit));
      ++i;
    }
    mztab.setMetaData(md);
  }

  void MzTabHelper::exportQuants(MzTab& mztab, 
    const PeptideAndProteinQuant::PeptideQuant& peptide_quants, 
    const PeptideAndProteinQuant::ProteinQuant& protein_quants,
    const ConsensusMap& consensus_map,
    const vector<ProteinIdentification::ProteinGroup>& indist_groups
    )
  {
    // Export proteins, the raw peptide features and the PSMs
    vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();
    mztab = MzTabHelper::exportConsensusMapToMzTab(consensus_map, "exported from ProteinQuantifier");

    // build map to retrieve detailed protein level information from consensus map
    map<String, ProteinHit> acc2prot_id;
    for (vector<ProteinIdentification>::const_iterator p_it = prot_ids.begin(); p_it != prot_ids.end(); ++p_it)
    {
      const vector<ProteinHit>& hits = p_it->getHits();
      for (vector<ProteinHit>::const_iterator h_it = hits.begin(); h_it != hits.end(); ++h_it)
      {
        acc2prot_id[h_it->getAccession()] = *h_it;
      }
    }
  
    // build map from protein group leader to group 
    map<String, ProteinIdentification::ProteinGroup> leader2group;
    for (vector<ProteinIdentification::ProteinGroup>::const_iterator p_it = indist_groups.begin(); p_it != indist_groups.end(); ++p_it)
    {
      const vector<String>& accessions = p_it->accessions;
      leader2group[accessions[0]] = *p_it;
    }

    MzTabProteinSectionRows rows;

    // determine number of channels
    Size n_study_variables = consensus_map.getFileDescriptions().size();

    // loop over protein quants
    for (PeptideAndProteinQuant::ProteinQuant::const_iterator q_it = protein_quants.begin(); q_it != protein_quants.end(); ++q_it)
    {
      // protein or protein group leader
      const String & protein_accession = q_it->first;
      const PeptideAndProteinQuant::ProteinData & protein_data = q_it->second;
      const PeptideAndProteinQuant::SampleAbundances & sample_abundances = protein_data.total_abundances;

      MzTabProteinSectionRow row;
      row.accession = MzTabString(protein_accession);

      String protein_group_type("single_protein");

      // check if single protein or group of indistinguishable proteins
      if (leader2group.find(protein_accession) != leader2group.end())
      {
        const ProteinIdentification::ProteinGroup& protein_group = leader2group.at(protein_accession);
        if (protein_group.accessions.size() >= 2)
        {
          protein_group_type = "indistinguishable_group";

          // add ambiguity members (exclude leader protein)
          String ambiguity_members;
          for (Size i = 1; i != protein_group.accessions.size(); ++i)
          {
            if (i != 1) ambiguity_members += ",";
            ambiguity_members += protein_group.accessions[i];
          }
          row.ambiguity_members.setSeparator(',');
          row.ambiguity_members.fromCellString(ambiguity_members);
        }
      }
      
      MzTabOptionalColumnEntry opt_column_entry;
      opt_column_entry.first = "opt_global_protein_group_type";
      opt_column_entry.second = protein_group_type;
      row.opt_.push_back(opt_column_entry);

      const ProteinHit& ph = acc2prot_id.at(protein_accession);
      row.protein_coverage = ph.getCoverage() >= 0 ? MzTabDouble(ph.getCoverage()) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
      row.description = ph.getDescription();

      // note: we put the protein score here as it makes more sense than the PSM score
      row.best_search_engine_score[1] = ph.getScore();

      // initialize columns
      for (Size study_variable = 1; study_variable <= n_study_variables; ++study_variable)
      {
        row.protein_abundance_stdev_study_variable[study_variable] = MzTabDouble();
        row.protein_abundance_std_error_study_variable[study_variable] = MzTabDouble();
        row.protein_abundance_study_variable[study_variable] = MzTabDouble();
        row.search_engine_score_ms_run[1][study_variable] = MzTabDouble();
      }

      // loop over sample abundances
      Size study_variable(1);
      for (PeptideAndProteinQuant::SampleAbundances::const_iterator a_it = sample_abundances.begin(); a_it != sample_abundances.end(); ++a_it)
      {
        row.protein_abundance_stdev_study_variable[study_variable];
        row.protein_abundance_std_error_study_variable[study_variable];
        row.protein_abundance_study_variable[study_variable] = MzTabDouble(a_it->second);
        ++study_variable;
      }
      rows.push_back(row);
    }
    mztab.setProteinSectionRows(rows);
  }
}

