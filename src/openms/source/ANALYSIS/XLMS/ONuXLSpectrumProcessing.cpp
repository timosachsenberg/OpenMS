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

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLDeisotoper.h>
#include <OpenMS/ANALYSIS/XLMS/ONuXLSpectrumProcessing.h>

namespace OpenMS
{
namespace OpenNuXL
{
  void ONuXLSpectrumProcessing::preprocess(
    PeakMap& exp, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    bool single_charge_spectra, 
    bool annotate_charge)
  {
    // filter MS2 map and remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // sort by mz
      exp[exp_index].sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(exp[exp_index], 
                                         fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
                                         1, 3, 
                                         false, 
                                         2, 10, 
                                         single_charge_spectra, 
                                         annotate_charge);
    #ifdef DEBUG_ONUXLSEARCH
      cout << "after deisotoping..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
//      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
      cout << endl;
    #endif

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);

    #ifdef DEBUG_ONUXLSEARCH
      cout << "after mower..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
    
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

    #ifdef DEBUG_ONUXLSEARCH
      cout << "after nlargest..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
 
      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
  
    #ifdef DEBUG_ONUXLSEARCH
      cout << "after sort..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
    }
//    MzMLFile().store(String("RNPxlSearch_a_") + String((int)annotate_charge) + ".mzML", exp);
  }
}
}