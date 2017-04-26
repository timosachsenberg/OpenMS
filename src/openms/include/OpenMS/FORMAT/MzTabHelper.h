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

#ifndef OPENMS_FORMAT_MZTABHELPER_H
#define OPENMS_FORMAT_MZTABHELPER_H

#include <OpenMS/FORMAT/MzTabHelper.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <map>
#include <vector>
#include <list>
#include <algorithm>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
/**
      @brief Common helper functions used by MzTab writing tools or converters.

      @ingroup FileIO
  */

  class OPENMS_DLLAPI MzTabHelper
  {
  public:
    /// Default constructor
    MzTabHelper();

    /// Destructor
    virtual ~MzTabHelper();
    
    /**
      @brief Gets peptide_evidences with data from internal structures adds their info to an MzTabPSMSectionRow (pre- or unfilled)

      @param peptide_evidences Vector of PeptideEvidence holding internal data.
      @param row Pre- or unfilled MzTabPSMSectionRow to be filled with the data.
      @param rows Vector of MzTabPSMSectionRow to add the differently updated rows to.
    */
    static void addPepEvidenceToRows(const std::vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows);
    
    /**
      @brief Inserts values from MetaInfoInterface objects matching a (precalculated or filtered) set of keys to optional columns of an MzTab row.
  
      @param keys Only values matching those keys will be extracted from the object inheriting from MetaInfoInterface.
      @param opt A vector of optional columns to add to.
      @param id The identifier for this optional value according to the mzTab standard (like global, MS_Run, assay, etc.)
      @param meta The object holding the MetaInformation (like PeptideHit, ProteinHit, etc.)
      @return void: Only updates the values of the columns in opt
    */
    static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String id, const MetaInfoInterface meta);
    
    /// @brief helper function to convert OpenMS mod names to MzTab mod strings    
    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const std::vector<String>& mods);
    
    /// @brief helper function to convert OpenMS mod names to MzTab mod strings    
    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const std::vector<String>& mods);

    /// @brief helper function to convert OpenMS mod names to MzTab mod strings    
    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const std::vector<String>& mods);
    
    /**
      @brief Generate MzTab style list of PTMs from AASequence object. 
      @param fixed_mods Fixed modifications are not reported.
      According to mzTab standard, fixed modifications are omitted in the PRT and PEP section.
      In contrast, all modifications are reported in the PSM section (see standard document for details).
    */
    static MzTabModificationList extractModificationListFromAASequence(const AASequence& aas, const std::vector<String>& fixed_mods);
    
    /**
       @brief Convert a feature map to MzTab
       @param feature_map Feature map to be converted
       @param filename Name of the featureXML to be annotated in the mzTab file as conversion source.
    */ 
    static MzTab exportFeatureMapToMzTab(const FeatureMap& feature_map, const String& filename);
    
    /**
       @brief Convert identifications to MzTab
       @param prot_ids ProteinIdentifications to be converted
       @param peptide_ids PeptideIdentifications to be converted
       @param filename Name of the identification file to be annotated in the mzTab file as conversion source.
    */ 
    static MzTab exportIdentificationsToMzTab(const std::vector<ProteinIdentification>& prot_ids, const std::vector<PeptideIdentification>& peptide_ids, const String& filename);
    
    /**
       @brief Convert a consensus map to MzTab
       @param consensus_map Consensus map to be converted
       @param filename Name of the consensusXML to be annotated in the mzTab file as conversion source.
    */ 
    static MzTab exportConsensusMapToMzTab(const ConsensusMap& consensus_map, const String& filename);
   
    /**
      @brief Reannotate MS runs to get proper link to the spectra
    */
    static void setMSRuns(const StringList& ms_run_locations, MzTab& mztab);
 
    /**
      @brief Add peptide and protein quantifications to an existing MzTab object
  
      @param mztab  MzTab object to be appended
    */
    static void exportQuants(MzTab& mztab, const PeptideAndProteinQuant::PeptideQuant& peptide_quants, const PeptideAndProteinQuant::ProteinQuant& protein_quants, const ConsensusMap& consensus_map);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZTABHELPER_H
