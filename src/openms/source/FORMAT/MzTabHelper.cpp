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

namespace OpenMS
{

  MzTabHelper::MzTabHelper()
  {

  }

  MzTabHelper::~MzTabHelper()
  {

  }
  
  void addPepEvidenceToRows(const std::vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows)
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
  
  void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String id, const MetaInfoInterface meta)
  {
    for (std::set<String>::const_iterator sit = keys.begin(); sit != keys.end(); ++sit)
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

  std::map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const std::vector<String>& mods)
  {
    std::map<Size, MzTabModificationMetaData> mods_mztab;
    for (std::vector<String>::const_iterator sit = mods.begin(); sit != mods.end(); ++sit)
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

  std::map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const std::vector<String>& mods)
  {
    if (mods.empty())
    {
      std::map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002454, No variable modifications searched, ]");
      mods_mztab.insert(std::make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  std::map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const std::vector<String>& mods)
  {
    if (mods.empty())
    {
      std::map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002453, No fixed modifications searched, ]");
      mods_mztab.insert(std::make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }


}
