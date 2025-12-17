// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <memory>
#include <vector>

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  class IsobaricQuantifierStatistics;
  class ConsensusMap;
  class ConsensusFeature;

  /**
    @brief Performs isotope impurity correction on intensities extracted from isobaric labeling experiments.

    This class implements algorithms for correcting isotope impurities in quantitative proteomics data
    obtained from isobaric labeling experiments such as iTRAQ or TMT. Isotope impurities arise from
    the fact that the reagents used for labeling are not 100% pure and contain isotopic variants that
    can contribute to the signal in neighboring channels.

    The correction is performed using a non-negative least squares (NNLS) approach, which solves the
    linear system Ax = b, where:
    - A is the correction matrix derived from the isotope impurity information provided by the reagent manufacturer
    - b is the vector of observed intensities in each channel
    - x is the vector of corrected intensities

    The NNLS approach ensures that the corrected intensities remain non-negative, which is physically
    meaningful for mass spectrometry data.

    @see IsobaricQuantitationMethod
    @see IsobaricQuantifierStatistics
  */
  class OPENMS_DLLAPI IsobaricIsotopeCorrector
  {
public:
    /**
     @brief Apply isotope correction to the given input map and store the corrected values in the output map.

     @param consensus_map_in The map containing the values that should be corrected.
     @param consensus_map_out The map where the corrected values should be stored.
     @param quant_method IsobaricQuantitationMethod (e.g., iTRAQ 4 plex)

     @throws Exception::FailedAPICall If the least-squares fit fails.
     @throws Exception::InvalidParameter If the given correction matrix is invalid.
     */
    static IsobaricQuantifierStatistics correctIsotopicImpurities(const ConsensusMap& consensus_map_in,
                                                                  ConsensusMap& consensus_map_out,
                                                                  const IsobaricQuantitationMethod* quant_method);

    /**
     @brief Apply isotope correction to a vector of channel intensities.

     This method applies the isotope correction directly to a vector of intensities representing
     the different isobaric channels. The vector is modified in-place to contain the corrected values.

     @param intensities Vector of channel intensities to be corrected (modified in-place)
     @param quant_method IsobaricQuantitationMethod providing the correction matrix (e.g., iTRAQ 4 plex)

     @throws Exception::FailedAPICall If the least-squares fit fails.
     @throws Exception::InvalidParameter If the given correction matrix is invalid.

     @note The size of the intensities vector must match the number of channels in the quantitation method.
     */
    static void
    correctIsotopicImpurities(std::vector<double>& intensities,
                              const IsobaricQuantitationMethod* quant_method);

    // No instance methods as this is a purely static class

private:
    /**
     * @brief Fills the input vectors for the NNLS step given the ConsensusFeature.
     *
     * @param[out] b Vector to be filled with intensities
     * @param[out] m_b OpenMS matrix to be filled with intensities (alternative representation)
     * @param[in] cf ConsensusFeature containing the channel intensities
     * @param[in] cm ConsensusMap containing the feature
     */
    static void fillInputVector_(std::vector<double>& b,
                                 Matrix<double>& m_b,
                                 const ConsensusFeature& cf,
                                 const ConsensusMap& cm);

    /**
     * @brief Extract channel intensities from a ConsensusFeature.
     */
    static std::vector<double> getIntensities_(const IsobaricQuantitationMethod* quant_method,
                                               const ConsensusFeature& cf,
                                               const ConsensusMap& cm);

    /**
     @brief Solve the non-negative least squares problem using OpenMS matrices.
     */
    static void solveNNLS_(const Matrix<double>& correction_matrix,
                           const Matrix<double>& m_b, Matrix<double>& m_x);

    /**
     @brief Solve the non-negative least squares problem using Matrix and vectors.
     */
    static void solveNNLS_(Matrix<double>& correction_matrix,
                           std::vector<double>& b,
                           std::vector<double>& x);

    /**
     @brief Compute statistics for the correction process.
     */
    static void computeStats_(const std::vector<double>& m_x,
                              const std::vector<double>& x_naive,
                              const float cf_intensity,
                              const IsobaricQuantitationMethod* quant_method,
                              IsobaricQuantifierStatistics& stats);

    /**
     @brief Compute statistics for the correction process using OpenMS matrices.
     */
    static void computeStats_(const Matrix<double>& m_x,
                              const std::vector<double>& x_naive,
                              const float cf_intensity,
                              const IsobaricQuantitationMethod* quant_method,
                              IsobaricQuantifierStatistics& stats);

    /**
     @brief Update the output consensus map with corrected intensities using std::vector.
     */
    static float updateOutputMap_(const ConsensusMap& consensus_map_in,
                                 ConsensusMap& consensus_map_out,
                                 Size current_cf,
                                 const std::vector<double>& m_x);

    /**
     @brief Update the output consensus map with corrected intensities using OpenMS Matrix.
     */
    static float updateOutputMap_(const ConsensusMap& consensus_map_in,
                                 ConsensusMap& consensus_map_out,
                                 Size current_cf,
                                 const Matrix<double>& m_x);
  };
} // namespace
