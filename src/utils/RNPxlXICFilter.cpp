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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLSpectrumProcessing.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/XLMS/ONuXLSpectrumProcessing.h>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_RNPxlXICFilter RNPxlXICFilter
    @brief Filters MS2 spectra based on XIC intensities in control and treatment. Used in RNPxl experiments to reduce candidate spectra.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ RNPxlXICFilter \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHires </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_RNPxl </td>
    </tr>
  </table>
</CENTER>

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RNPxlXICFilter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RNPxlXICFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRNPxlXICFilter :
  public TOPPBase
{
public:
  TOPPRNPxlXICFilter() :
    TOPPBase("RNPxlXICFilter", "Remove MS2 spectra from treatment based on the fold change between control and treatment.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // input files

    registerInputFile_("control", "<file>", "", "input mzML file");
    setValidFormats_("control", ListUtils::create<String>("mzML"));
    registerInputFile_("treatment", "<file>", "", "input mzML file");
    setValidFormats_("treatment", ListUtils::create<String>("mzML"));

    registerDoubleOption_("fold_change", "", 2.0, "fold change between XICs", false, false);
    registerDoubleOption_("rt_tol", "", 20, "RT tolerance in [s] for finding max peak (whole RT range around RT middle)", false, false);
    registerDoubleOption_("mz_tol", "", 10, "m/z tolerance in [ppm] for finding a peak", false, false);

    // output files
    registerOutputFile_("out", "<file>", "", "output of the treatment file after XIC filtering.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));    
  }


  ExitCodes main_(int, const char**) override
  {
    // Parameter parsing
    const string control_mzml(getStringOption_("control"));
    const string treatment_mzml(getStringOption_("treatment"));

    const string out_mzml(getStringOption_("out"));
    const double mz_tolerance_ppm = getDoubleOption_("mz_tol");
    const double fold_change = getDoubleOption_("fold_change");
    const double rt_tolerance_s = getDoubleOption_("rt_tol");

    // load experiments
    PeakMap exp_control;
    MzMLFile mzml_file;
    mzml_file.load(control_mzml, exp_control);

    PeakMap exp_treatment;
    mzml_file.load(treatment_mzml, exp_treatment);

   OpenNuXL::ONuXLSpectrumProcessing::filterByFoldChange(
                          exp_control, 
                          exp_treatment,
                          rt_tolerance_s, 
                          mz_tolerance_ppm, 
                          fold_change);

    mzml_file.store(out_mzml, exp_treatment);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPRNPxlXICFilter tool;
  return tool.main(argc, argv);
}
