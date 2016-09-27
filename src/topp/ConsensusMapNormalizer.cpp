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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>

#include <Eigen/Dense>

using namespace OpenMS;
using namespace std;
using namespace Eigen;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_ConsensusMapNormalizer ConsensusMapNormalizer

    @brief Normalization of intensities in a set of maps using robust regression.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ConsensusMapNormalizer \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
    </table>
</CENTER>

The tool normalizes the intensities of a set of maps (consensusXML file). The following normalization algorithms are available:

- Robust regression: Maps are normalized pair-wise relative to the map with the most features. Given two maps, peptide featues are classified as non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) or outliers. From the non-outliers an average intensity ratio is calculated and used for normalization.

- Median correction: The median of all maps is set to the median of the map with the most features.

- Quantile normalization: Performs an exact quantile normalization if the number of features is equal across all maps. Otherwise, an approximate quantile normalization using resampling is applied.

<B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_ConsensusMapNormalizer.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_ConsensusMapNormalizer.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPConsensusMapNormalizer :
  public TOPPBase
{

public:
  TOPPConsensusMapNormalizer() :
    TOPPBase("ConsensusMapNormalizer", "Normalizes maps of one consensusXML file")
  {
  }

protected:
/*
#' SVT Imputation
#'
#' Imputation using Singular Value Thresholding a la Cai, Candes, Shen.
#' @param x a data frame or matrix with size n1 x n2 where each row represents a different record
#' @param lambda the penalty on the singular values
#' @param stepsize optional.  If not provided, uses 1.2 * (n1 * n2) / (number of missing elements)
#' @param threshold convergence threshold
#' @param max.iters maximum number of iterations.  Note that each iteration will require computing
#'   an SVD
#' @param verbose if TRUE print status updates
#' @references A Singular Value Thresholding Algorithm for Matrix Completion. Cai, Candes, Shen. 
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   SVTImpute(x, 3)
#' @export
SVTImpute = function(x, lambda, stepsize, threshold = 1e-3, max.iters = 10, verbose=F) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix

  x.zeroimpute = x
  x.zeroimpute[missing.matrix] = 0
  x.zeroimpute.norm = norm(x.zeroimpute, "F")

  if (missing(stepsize)) stepsize = min(1.2 * length(x) / sum(missing.matrix), 1.9)
  k = ceiling(lambda / ( stepsize * norm(x.zeroimpute, "F")))

  Y = k*stepsize*x.zeroimpute

  for (iter in 1:max.iters) {
    if (verbose) print(paste("Begin iteration", iter))
    y.svd = svd(Y)
    lambda.indices = which(y.svd$d < lambda)
    if(length(lambda.indices) > 0) {
      d.augmented = c(y.svd$d[-lambda.indices] - lambda, 
                      rep(0, length(lambda.indices)))
    } else {
      d.augmented = y.svd$d - lambda
    }
    X.k = (y.svd$u %*% diag(d.augmented) %*% t(y.svd$v))
    X.k.temp = X.k; X.k.temp[missing.matrix] = 0
    if (norm(X.k.temp - x.zeroimpute, "F") / x.zeroimpute.norm < threshold) {
      if (verbose) print(paste("Converging on iteration", iter))
      break
    }
    Y[!missing.matrix] = Y[!missing.matrix] + stepsize*(x[!missing.matrix] - X.k[!missing.matrix])
    Y[missing.matrix] = 0
  }

  return( list(
    x=X.k,
    missing.matrix = missing.matrix
  ))
}
*/
  static void imputeSVTSingleLambda_(MatrixXd & X, double lambda, Size max_iter = 1000, double stepsize = 0, double threshold = 1e-5)
  {
    // replace missing values by zeros to generate matrix Z
    MatrixXd X_Z(X);  

    vector<pair<Size, Size> > missing_values;
    for (int i = 0; i < X.cols(); ++i)
    {
      for (int j = 0; j < X.rows(); ++j)
      {
        if (std::isnan(X(j, i)) || X(j, i) == 0)
        {
          X_Z(j, i) = 0.0;
          missing_values.push_back(make_pair(j, i)); 
        }
      }
    }

    // if no values are missing just return
    if (missing_values.empty()) return;

    // calculate Frobenius norm of zero-imputed X
    double nZ = X_Z.norm();

    // set stepsize as default to min(1.2 * (nrows*ncols) / numnber of missing values, 1.9)
    if (stepsize <= 0.0) stepsize = std::min(1.2 * (X.rows() * X.cols()) / missing_values.size(), 1.9); 

    // k = ceiling(lambda / (stepsize * nZ))
    double k = ceil(lambda / (stepsize * nZ));    

    // scale zero imputed matrix Z to obtain Y: Y = k*stepsize*Z
    MatrixXd Y(k * stepsize * X_Z);

    MatrixXd X_k(X);

    // iterate max_iter times  
    for (Size i = 0; i != max_iter; ++i)
    { 
      //cout << "Begin iteration: " << i << endl;

      // perform SVD on Y
      JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV); 
      VectorXd sigma = svd.singularValues();

      // subtract lambda from entries > lambda and append zeros otherwise (=truncate to zero)
      // then reshape singular value vector into diagonal matrix
      MatrixXd D((sigma.array() >= lambda).select(sigma.array() - lambda, 0.0).matrix().asDiagonal());

      X_k = svd.matrixU() * D * svd.matrixV().transpose();

      MatrixXd X_k_temp(X_k); 
      for (Size m = 0; m != missing_values.size(); ++m)
      {
        X_k_temp(missing_values[m].first, missing_values[m].second) = 0.0;
      }
   
      if ((X_k_temp - X_Z).norm() / nZ < threshold) // converged?
      {
        break;
      }

      Y += stepsize * (X - X_k);

      for (Size m = 0; m != missing_values.size(); ++m)
      {
        Y(missing_values[m].first, missing_values[m].second) = 0.0;
      }
    }
  
    // return X
    X = X_k;
  }

  static double imputeSVTChooseLambda_(MatrixXd & m, Size max_iter = 1000, double min = 0, double max = 1000, double step = 10)
  {
    // find best lambda
    // 1. randomly erase 1/3 of valid entries and calculate error for different lambdas
    vector<pair<int, int> > valid_values;
    for (int i = 0; i < m.cols(); ++i)
    {
      for (int j = 0; j < m.rows(); ++j)
      {
        double v = m(j, i);
        if (!std::isnan(v) && v != 0)
        {
          valid_values.push_back(make_pair(j, i)); 
        }
      }
    }

    // no imputation possible as no valid values given - return zeros
    if (valid_values.empty()) 
    {
      m.setZero();
    }

    std::random_shuffle(valid_values.begin(), valid_values.end());

    valid_values.resize(valid_values.size() / 3);
   
    MatrixXd sample(m);
    for (Size i = 0; i != valid_values.size(); ++i)
    {
      sample(valid_values[i].first, valid_values[i].second) = 0; // artifically zero out valid value
    }

    double best_rmse(std::numeric_limits<double>::max());
    double best_lambda(min);
    MatrixXd best_matrix(m);
    imputeSVTSingleLambda_(best_matrix, 0); // generate a valid matrix

    // 2. choose lambda with lowest error
    for (double l = min; l <= max; l += step)
    { 
      // impute using current lambda
      MatrixXd m_imputed(m);
      imputeSVTSingleLambda_(m_imputed, l);

      // calculate RMSE on zero-ed out values
      double se(0); // sum of squared errors
      for (Size i = 0; i != valid_values.size(); ++i)
      {
        double error = m_imputed(valid_values[i].first, valid_values[i].second) - m(valid_values[i].first, valid_values[i].second);
        se += error * error;
      }

      double rmse = sqrt(se / valid_values.size());

      if (rmse < best_rmse)
      {
        best_rmse = rmse;
        best_lambda = l;
        best_matrix = m_imputed;
      }
    }

    m = best_matrix;
    return best_rmse;
  }

  void imputeSVT_(ConsensusMap & cm)
  {
    // copy consensus map to eigen matrix
    Size nrows(cm.size());
    Size ncols(cm.getFileDescriptions().size());

    MatrixXd m = MatrixXd::Zero(nrows, ncols);
    
    for (Size j = 0; j != cm.size(); ++j)
    {
      const ConsensusFeature & c = cm[j];
      const ConsensusFeature::HandleSetType & fs = c.getFeatures();      
      for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
      {
        Size column_index = fit->getMapIndex();
        m(j, column_index) = log10(fit->getIntensity() + 1); // log10(x+1) transform
      } 
    }

    // determine missing values in each row (=consensus feature)
    Size missing(0);
    map<Size, set<Size> > row_to_missing_column;
    for (int i = 0; i < m.cols(); ++i)
    {
      for (int j = 0; j < m.rows(); ++j)
      {
        if (std::isnan(m(j, i)) || m(j, i) == 0)
        {
          row_to_missing_column[j].insert(i);
          ++missing;
        }
      }
    }

    cout << "Imputing missing values: " << missing << endl;

    double rmse = imputeSVTChooseLambda_(m); 

    cout << "Imputation finished. RMSE of log10 values: " << rmse << endl;

    // copy back to consensus map
    for (Size j = 0; j != cm.size(); ++j)
    {
      // current consensus feature had missing values
      if (row_to_missing_column.find(j) != row_to_missing_column.end())
      {
        const set<Size> & missing = row_to_missing_column[j];
        ConsensusFeature & c = cm[j];
        for (set<Size>::iterator it = missing.begin(); it != missing.end(); ++it)
        {
          Size missing_index = *it; // missing column
          Feature f;
          f.setRT(c.getRT());
          f.setMZ(c.getMZ());
          f.setIntensity(pow(10.0, m(j, missing_index)) - 1.0); // transform log10(x+1) back
          FeatureHandle fh(missing_index, f);
          c.insert(fh);
        }
      }
    }
       
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    addEmptyLine_();
    registerStringOption_("algorithm_type", "<type>", "robust_regression", "The normalization algorithm that is applied. 'robust_regression' scales each map by a fator computed from the ratios of non-differential background features (as determined by the ratio_threshold parameter), 'quantile' performs quantile normalization, 'median' scales all maps to the same median intensity, 'median_shift' shifts the median instead of scaling (WARNING: if you have regular, log-normal MS data, 'median_shift' is probably the wrong choice. Use only if you know what you're doing!)", false, false);
    setValidStrings_("algorithm_type", ListUtils::create<String>("robust_regression,median,median_shift,quantile"));
    registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "Only for 'robust_regression': the parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers.", false);
    setMinFloat_("ratio_threshold", 0.001);
    setMaxFloat_("ratio_threshold", 1.0);
    registerStringOption_("accession_filter", "<regexp>", "", "Use only features with accessions (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);
    registerStringOption_("description_filter", "<regexp>", "", "Use only features with description (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);
    registerStringOption_("imputation_method", "<type>", "none", "Imputation method used.", false, false);
    setValidStrings_("imputation_method", ListUtils::create<String>("none,SVT"));
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String algo_type = getStringOption_("algorithm_type");
    String acc_filter = getStringOption_("accession_filter");
    String desc_filter = getStringOption_("description_filter");
    double ratio_threshold = getDoubleOption_("ratio_threshold");
    String imputation_method = getStringOption_("imputation_method");

    ConsensusXMLFile infile;
    infile.setLogType(log_type_);
    ConsensusMap map;
    infile.load(in, map);

    if (imputation_method == "SVT")
    {
      imputeSVT_(map);
      //addDataProcessing_(map, getProcessingInfo_(DataProcessing::IMPUTATION));
    }


    //map normalization
    if (algo_type == "robust_regression")
    {
      map.sortBySize();
      vector<double> results = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, ratio_threshold, acc_filter, desc_filter);
      ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, results);
    }
    else if (algo_type == "median")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, acc_filter, desc_filter);
    }
    else if (algo_type == "median_shift")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SHIFT, acc_filter, desc_filter);
    }
    else if (algo_type == "quantile")
    {
      if (acc_filter != "" || desc_filter != "")
      {
        LOG_WARN << endl << "NOTE: Accession / description filtering is not supported in quantile normalization mode. Ignoring filters." << endl << endl;
      }
      ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(map);
    }
    else
    {
      cerr << "Unknown algorithm type  '" << algo_type.c_str() << "'." << endl;
      return ILLEGAL_PARAMETERS;
    }

    //annotate output with data processing info and save output file
    addDataProcessing_(map, getProcessingInfo_(DataProcessing::NORMALIZATION));
    infile.store(out, map);

    return EXECUTION_OK;
  }
};


int main(int argc, const char ** argv)
{
  TOPPConsensusMapNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
