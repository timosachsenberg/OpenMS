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
#include <OpenMS/CONCEPT/LogStream.h>


#include <vector>
#include <set>


using namespace OpenMS;
using namespace OpenNuXL;
using namespace std;


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

void ONuXLSpectrumProcessing::filterByFoldChange(
                        const PeakMap& exp_control, 
                        PeakMap& exp_treatment,
                        double rt_tolerance_s, 
                        double mz_tolerance_ppm, 
                        double fold_change)
{
  // extract precursor mz and rts
  vector<double> pc_mzs;
  vector<double> pc_ms2_rts;
  for (Size i = 0; i != exp_treatment.size(); ++i)
  {
    if (exp_treatment[i].getMSLevel() == 2)
    {
      if (!exp_treatment[i].getPrecursors().empty())
      {
        double pc_mz = exp_treatment[i].getPrecursors()[0].getMZ();
        double ms2_rt = exp_treatment[i].getRT(); // use rt of MS2
        pc_mzs.push_back(pc_mz);
        pc_ms2_rts.push_back(ms2_rt);
      }
    }
  }

  vector<double> control_XIC_larger_rts;
  vector<double> treatment_XIC_larger_rts;
  vector<double> indifferent_XICs_rts;

  // search for each EIC and add up
  for (Size i = 0; i < pc_mzs.size(); ++i)
  {
    //cerr << "start" << endl;
    double pc_ms2_rt = pc_ms2_rts[i];
    double pc_mz = pc_mzs[i];

    //std::cerr << "Rt" << cm[i].getRT() << "  mz: " << cm[i].getMZ() << " R " <<  cm[i].getMetaValue("rank") << "\n";

    double mz_da = mz_tolerance_ppm * pc_mzs[i] / 1e6; // mz tolerance in Dalton
    double rt_start = pc_ms2_rts[i] - rt_tolerance_s / 2.0;

    // get area iterator (is MS1 only!) for rt and mz window
    PeakMap::ConstAreaIterator it1 = 
      exp_control.areaBeginConst(pc_ms2_rt - rt_tolerance_s / 2, pc_ms2_rt + rt_tolerance_s / 2, pc_mz - mz_da, pc_mz  + mz_da);
    PeakMap::ConstAreaIterator it2 = 
      exp_treatment.areaBeginConst(pc_ms2_rt - rt_tolerance_s / 2, pc_ms2_rt + rt_tolerance_s / 2, pc_mz - mz_da, pc_mz  + mz_da);

    // determine maximum number of MS1 scans in retention time window
    set<double> rts1;
    set<double> rts2;
    for (; it1 != exp_control.areaEndConst(); ++it1)
    {
      rts1.insert(it1.getRT());
    }

    for (; it2 != exp_treatment.areaEndConst(); ++it2)
    {
      rts2.insert(it2.getRT());
    }

    Size length = std::max(rts1.size(), rts2.size()) / 2.0;

    //cout << length << endl;
    if (length == 0)
    {
      LOG_INFO << "WARNING: no MS1 scans in retention time window found in both maps (mz: " << pc_mzs[i] << " / rt: " << pc_ms2_rts[i] << ")" << endl;
      continue;
    }

    vector<double> XIC1(length, 0.0);
    vector<double> XIC2(length, 0.0);

    it1 = exp_control.areaBeginConst(pc_ms2_rt - rt_tolerance_s / 2.0, pc_ms2_rt + rt_tolerance_s / 2.0, pc_mz - mz_da, pc_mz + mz_da);
    it2 = exp_treatment.areaBeginConst(pc_ms2_rt - rt_tolerance_s / 2.0, pc_ms2_rt + rt_tolerance_s / 2.0, pc_mz - mz_da, pc_mz + mz_da);

    for (; it1 != exp_control.areaEndConst(); ++it1)
    {
      double relative_rt = (it1.getRT() - rt_start) / rt_tolerance_s;
      Size bin = relative_rt * (length - 1);
      XIC1[bin] += it1->getIntensity();
      if (bin >= length)
      {
        bin = length - 1;
      }

    }

    for (; it2 != exp_treatment.areaEndConst(); ++it2)
    {
      double relative_rt = (it2.getRT() - rt_start) / rt_tolerance_s;
      Size bin = relative_rt * (length - 1);
      if (bin >= length)
      {
        bin = length - 1;
      }
      XIC2[bin] += it2->getIntensity();
    }

    double total_itensity1 = std::accumulate(XIC1.begin(), XIC1.end(), 0.0);
    double total_itensity2 = std::accumulate(XIC2.begin(), XIC2.end(), 0.0);

    double ratio = total_itensity2 / (total_itensity1 + 1);

    if (ratio < 1.0 / fold_change)
    {
      control_XIC_larger_rts.push_back(pc_ms2_rt);
    }
    else if (ratio > fold_change)
    {
      treatment_XIC_larger_rts.push_back(pc_ms2_rt);
    }
    else
    {
      indifferent_XICs_rts.push_back(pc_ms2_rt);
      continue;
    }
  }

  LOG_INFO << "control larger: " << control_XIC_larger_rts.size() 
           << " treatment larger: " << treatment_XIC_larger_rts.size() 
           << " indifferent: " << indifferent_XICs_rts.size() << endl;

  PeakMap exp_out = exp_treatment;
  exp_out.clear(false); // don't clear meta-data

  for (Size i = 0; i != exp_treatment.size(); ++i)
  {
    Size ms_level = exp_treatment[i].getMSLevel();

    if (ms_level == 1)
    {
      exp_out.addSpectrum(exp_treatment[i]);
      continue;
    }
    else if (ms_level == 2)
    {
      // determine if pc is in list -> passed
      double rt = exp_treatment[i].getRT();
      for (Size j = 0; j != treatment_XIC_larger_rts.size(); ++j)
      {
        if (fabs(rt - treatment_XIC_larger_rts[j]) <= 0.001)
        {
          exp_out.addSpectrum(exp_treatment[i]);
          break;
        }
      }
    }
  }

  exp_treatment.swap(exp_out);

}
