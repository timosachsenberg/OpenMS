// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

#include <vector>
#include <map>

using std::vector;
using std::map;

namespace OpenMS
{
  // static
  void ModifiedPeptideGenerator::applyFixedModifications(
    const vector<Ribonucleotide>::const_iterator& fixed_mods_begin,
    const vector<Ribonucleotide>::const_iterator& fixed_mods_end,
    NASequence& seq)
  {
    // set terminal modifications for modifications without amino acid preference
    for (vector<Ribonucleotide>::const_iterator fixed_it = fixed_mods_begin; fixed_it != fixed_mods_end; ++fixed_it)
    {
      if (fixed_it->getTermSpecificity() == Ribonucleotide::N_TERM)
      {
        if (!seq.hasNTerminalModification())
        {
          seq.setNTerminalModification(fixed_it);
        }
      }
      else if (fixed_it->getTermSpecificity() == Ribonucleotide::C_TERM)
      {
        if (!seq.hasCTerminalModification())
        {
          seq.setCTerminalModification(fixed_it);
        }
      }
    }

    //iterate over each ribo
    size_t residue_index(0);
    for (auto const & r : seq)
    {
      // skip already modified residue
      if (r.isModified()) { ++residue_index; continue; }

      //set fixed modifications
      for (vector<Ribonucleotide>::const_iterator fixed_it = fixed_mods_begin; fixed_it != fixed_mods_end; ++fixed_it)
      {
        // check if amino acid match between modification and current residue
        if (residue_it->getOneLetterCode()[0] != fixed_it->getOrigin())
        {
          ++residue_index;
          continue;
        }

        // Term specificity is ANYWHERE on the seq, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
        const RibonucleotideModification::TermSpecificity& term_spec = fixed_it->getTermSpecificity();
        if (term_spec == RibonucleotideModification::ANYWHERE)
        {
          seq.setModification(residue_index, fixed_it->getFullName());
        }
        else if (term_spec == RibonucleotideModification::C_TERM && residue_index == (seq.size() - 1))
        {
          seq.setCTerminalModification(fixed_it->getFullName());
        }
        else if (term_spec == RibonucleotideModification::N_TERM && residue_index == 0)
        {
          seq.setNTerminalModification(fixed_it->getFullName());
        }
      }
      ++residue_index;
    }
  }

  // static
  void ModifiedPeptideGenerator::applyVariableModifications(
    const vector<Ribonucleotide>::const_iterator& var_mods_begin,
    const vector<Ribonucleotide>::const_iterator& var_mods_end,
    const NASequence& seq,
    size_t max_variable_mods_per_seq,
    vector<NASequence>& all_modified_seqs,
    bool keep_unmodified)
  {
    // no variable modifications specified or no variable mods allowed? no compatibility map needs to be build
    if (var_mods_begin == var_mods_end || max_variable_mods_per_seq == 0)
    {
      // if unmodified seqs should be kept return the original list of digested seqs
      if (keep_unmodified)
      {
        all_modified_seqs.push_back(seq);
      }
      return;
    }

    // if there is at most one variable modification allowed for a seq we don't need combinatoric placement and can reside to a faster implementation
    if (max_variable_mods_per_seq == 1)
    {
      applyAtMostOneVariableModification_(var_mods_begin,
                                          var_mods_end,
                                          seq,
                                          all_modified_seqs,
                                          keep_unmodified);
      return;
    }

    const int N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    const int C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    //keep a list of all possible modifications of this seq
    vector<NASequence> modified_seqs;

    // only add unmodified version if flag is set (default)
    if (keep_unmodified) { modified_seqs.push_back(seq); }

    //iterate over each residue and build compatibility mapping describing
    //which amino acid (seq index) is compatible with which modification
    map<int, vector<Ribonucleotide> > map_compatibility;

    // set terminal modifications for modifications without amino acid preference
    for (vector<Ribonucleotide>::const_iterator variable_it = var_mods_begin; variable_it != var_mods_end; ++variable_it)
    {
      if (variable_it->getTermSpecificity() == Ribonucleotide::N_TERM)
      {
        if (!seq.hasNTerminalModification())
        {
          map_compatibility[N_TERM_MODIFICATION_INDEX].push_back(*variable_it);
        }
      }
      else if (variable_it->getTermSpecificity() == Ribonucleotide::C_TERM)
      {
        if (!seq.hasCTerminalModification())
        {
          map_compatibility[C_TERM_MODIFICATION_INDEX].push_back(*variable_it);
        }
      }
    }

    size_t residue_index(0);
    for (auto const & r : seq)
    {
      // skip already modified residues
      if (r.isModified()) { ++residue_index; continue; }

      //determine compatibility of variable modifications
      for (vector<Ribonucleotide>::const_iterator variable_it = var_mods_begin; variable_it != var_mods_end; ++variable_it)
      {
        // check if amino acid match between modification and current residue
        if (residue_it->getOneLetterCode()[0] != variable_it->getOrigin()) { ++residue_index; continue; }

        // Term specificity is ANYWHERE on the seq, C_TERM or N_TERM
        // (currently no explicit support in OpenMS for protein C-term and
        // protein N-term)
        const Ribonucleotide::TermSpecificity& term_spec = variable_it->getTermSpecificity();
        if (term_spec == Ribonucleotide::ANYWHERE)
        {
          map_compatibility[static_cast<int>(residue_index)].push_back(*variable_it);
        }
        else if (term_spec == Ribonucleotide::C_TERM && residue_index == (seq.size() - 1))
        {
          map_compatibility[C_TERM_MODIFICATION_INDEX].push_back(*variable_it);
        }
        else if (term_spec == Ribonucleotide::N_TERM && residue_index == 0)
        {
          map_compatibility[N_TERM_MODIFICATION_INDEX].push_back(*variable_it);
        }
      }
      ++residue_index;
    }

    // Check if no compatible site that can be modified by variable
    // modification. If so just return seqs without variable modifications.
    const size_t & compatible_mod_sites = map_compatibility.size();
    if (compatible_mod_sites == 0)
    {
      if (keep_unmodified) { all_modified_seqs.push_back(seq); }
      return;
    }

    // generate powerset of max_variable_mods_per_seq sized subset of all compatible modification sites
    size_t max_placements = std::min(max_variable_mods_per_seq, compatible_mod_sites);
    for (size_t n_var_mods = 1; n_var_mods <= max_placements; ++n_var_mods)
    {
      // enumerate all modified seqs with n_var_mods variable modified residues
      size_t zeros = std::max((size_t)0, compatible_mod_sites - n_var_mods);
      vector<bool> subset_mask;

      for (size_t i = 0; i != compatible_mod_sites; ++i)
      {
        // create mask 000011 to select last (e.g. n_var_mods = 2) two compatible sites as subset from the set of all compatible sites
        if (i < zeros)
        {
          subset_mask.push_back(false);
        }
        else
        {
          subset_mask.push_back(true);
        }
      }

      // generate all subsets of compatible sites {000011, ... , 101000, 110000} with current number of allowed variable modifications per seq
      do
      {
        // create subset indices e.g.{4,12} from subset mask e.g. 1010000 corresponding to the positions in the seq sequence
        vector<int> subset_indices;
        map<int, vector<Ribonucleotide> >::const_iterator mit = map_compatibility.begin();
        for (size_t i = 0; i != compatible_mod_sites; ++i, ++mit)
        {
          if (subset_mask[i])
          {
            subset_indices.push_back(mit->first);
          }
        }

        // now enumerate all modifications
        recurseAndGenerateVariableModifiedSequences_(subset_indices, map_compatibility, 0, seq, modified_seqs);
      }
      while (next_permutation(subset_mask.begin(), subset_mask.end()));
    }
    // add modified version of the current seq to the list of all seqs
    all_modified_seqs.insert(all_modified_seqs.end(), modified_seqs.begin(), modified_seqs.end());
  }


  // static
  void ModifiedPeptideGenerator::recurseAndGenerateVariableModifiedPeptides_(
    const vector<int>& subset_indices,
    const map<int, vector<Ribonucleotide> >& map_compatibility,
    int depth,
    const NASequence& current_seq,
    vector<NASequence>& modified_seqs)
  {
    const int N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    const int C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    // cout << depth << " " << subset_indices.size() << " " << current_seq.toString() << endl;

    // end of recursion. Add the modified seq and return
    if (depth == (int)subset_indices.size())
    {
      modified_seqs.push_back(current_seq);
      return;
    }

    // get modifications compatible to residue at current seq position
    const int current_index = subset_indices[depth];

    map<int, vector<Ribonucleotide> >::const_iterator pos_mod_it = map_compatibility.find(current_index);
    const vector<Ribonucleotide>& mods = pos_mod_it->second; // we don't need to check for .end as entry is guaranteed to exist

    for (auto const & m : mods)
    {
      // copy seq and apply modification
      NASequence new_seq = current_seq;
      if (current_index == C_TERM_MODIFICATION_INDEX)
      {
        new_seq.setCTerminalModification(m.getFullName());
      }
      else if (current_index == N_TERM_MODIFICATION_INDEX)
      {
        new_seq.setNTerminalModification(m.getFullName());
      }
      else
      {
        new_seq.setModification(current_index, m.getFullName());
      }

      // recurse with modified seq
      recurseAndGenerateVariableModifiedSequences_(subset_indices, map_compatibility, depth + 1, new_seq, modified_seqs);
    }
  }

  // static
  void ModifiedPeptideGenerator::applyAtMostOneVariableModification_(
    const vector<Ribonucleotide>::const_iterator& var_mods_begin,
    const vector<Ribonucleotide>::const_iterator& var_mods_end,
    const NASequence& seq,
    vector<NASequence>& all_modified_seqs,
    bool keep_unmodified)
  {
    if (keep_unmodified)
    {
      all_modified_seqs.push_back(seq);
    }

    // we want the same behavior as for the slower function... we would need a reverse iterator here that NASequence doesn't provide
    for (NASequence::ConstIterator residue_it = seq.end() - 1; residue_it != seq.begin() - 1; --residue_it)
    {
      // skip already modified residues
      if (residue_it->isModified())
      {
        continue;
      }

      size_t residue_index = residue_it - seq.begin();

      //determine compatibility of variable modifications
      for (vector<Ribonucleotide>::const_iterator variable_it = var_mods_begin; variable_it != var_mods_end; ++variable_it)
      {
        // check if amino acid match between modification and current residue
        if (residue_it->getOneLetterCode()[0] != variable_it->getOrigin()) { continue; }

        bool is_compatible = false;

        // Term specificity is ANYWHERE on the seq, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
        const Ribonucleotide::TermSpecificity& term_spec = variable_it->getTermSpecificity();
        if (term_spec == Ribonucleotide::ANYWHERE)
        {
          is_compatible = true;
        }
        else if (term_spec == Ribonucleotide::C_TERM && residue_index == (seq.size() - 1))
        {
          is_compatible = true;
        }
        else if (term_spec == Ribonucleotide::N_TERM && residue_index == 0)
        {
          is_compatible = true;
        }

        // residue modification an be placed at current position? Then generate modified seq.
        if (is_compatible)
        {
          NASequence new_seq = seq;
          new_seq.setModification(residue_index, variable_it->getFullName());
          all_modified_seqs.push_back(new_seq);
        }
      }
    }
  }
}
