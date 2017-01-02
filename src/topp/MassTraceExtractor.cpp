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
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_MassTraceExtractor MassTraceExtractor

        @brief MassTraceExtractor extracts mass traces from a @ref MSExperiment map and stores them into a @ref FeatureXMLFile.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MassTraceExtractor \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
        </table>
        </CENTER>


        This TOPP tool detects mass traces in centroided LC-MS maps and stores them as features in
        a @ref FeatureMap. These features may be either used directly as input for an metabolite ID approach or further
        be assembled to aggregate features according to a theoretical isotope pattern. For metabolomics experiments,
        the @ref TOPP_FeatureFinderMetabo tool offers both mass trace extraction and isotope pattern assembly.
        For proteomics data, please refer to the @ref TOPP_FeatureFinderCentroided tool.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_MassTraceExtractor.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMassTraceExtractor :
  public TOPPBase
{
public:
  TOPPMassTraceExtractor() :
    TOPPBase("MassTraceExtractor", "Detects mass traces in centroided LC-MS data.")
  {
  }

protected:
  struct MassTracesStatistics_
  {
    double sd_ppm_median;
    double sd_ppm_lMAD;
    double sd_ppm_rMAD;
    double mz_wiggle_median;
    double mz_wiggle_lMAD;
    double mz_wiggle_rMAD;
  };

  static void calculateQuantiles_(const vector<double>& v, const String& header, double& median, double& mad_l, double& mad_r)
  {
    vector<double> ints(v);
    std::sort(ints.begin(), ints.end());

    LOG_INFO << header << ints.size() << "\n";
    if (!ints.empty())
    {
      LOG_INFO << "             min: " << ints[0] << "\n"
               << "     1% quantile: " << ints[ints.size() * 1 / 100] << "\n"
               << "     5% quantile: " << ints[ints.size() * 5 / 100] <<"\n"
               << "    10% quantile: " << ints[ints.size() * 1 / 10] << "\n"
               << "    20% quantile: " << ints[ints.size() * 2 / 10] <<"\n"
               << "    30% quantile: " << ints[ints.size() * 3 / 10]<< "\n"
               << "    40% quantile: " << ints[ints.size() * 4 / 10]<< "\n"
               << "          median: " << ints[ints.size() * 5 / 10]<< "\n"
               << "    60% quantile: " << ints[ints.size() * 6 / 10]<< "\n"
               << "    70% quantile: " << ints[ints.size() * 7 / 10]<< "\n"
               << "    80% quantile: " << ints[ints.size() * 8 / 10]<< "\n"
               << "    90% quantile: " << ints[ints.size() * 9 / 10]<< "\n"
               << "    95% quantile: " << ints[ints.size() * 95 / 100]<< "\n"
               << "    99% quantile: " << ints[ints.size() * 99 / 100]<< "\n"
               << "             max: " << ints[ints.size() - 1] << "\n";        
      median = ints[ints.size() * 5 / 10];
      vector<double>::iterator int_mid = ints.begin() + ints.size() * 0.5;

      mad_l = Math::MAD(ints.begin(), int_mid, median);
      mad_r = Math::MAD(int_mid, ints.end(), median);

      LOG_INFO << "left / right MAD  : " << mad_l << "\t" << mad_r << std::endl;
    }
  }

  static void getPeakMapStatistics_(MSExperiment<Peak1D> & peak_map, double ppm)
  {
    // overall intensity distribution
    vector<double> ints;
    for (Size i = 0; i != peak_map.size(); ++i)
    {
      for (Size j = 0; j != peak_map[i].size(); ++j)
      {
        ints.push_back(peak_map[i][j].getIntensity());
      }
    }
    double median, t2, rMAD;
    calculateQuantiles_(ints, "MS run intensities: ", median, t2, rMAD);

    // collect intensities from potential isotopic pairs
    vector<double> iso_ints;

    // consider up to charge 10 isotopic patterns
    vector<double> mz_shifts;
    for (Size z = 1; z <= 10; ++z)
    {
      for (Size n = 1; n <= 2; ++n)  
      {
        Size gcd = Math::gcd(n, z);
        double shift = (double)(n/gcd) / (double)(z / gcd) * Constants::C13C12_MASSDIFF_U;
        if (std::find(mz_shifts.begin(), mz_shifts.end(), shift) == mz_shifts.end()) mz_shifts.push_back(shift);
      }
    }
    
    for (Size i = 2; i < peak_map.size() - 2; ++i)
    {
      vector<Size> indices;
      for (Size j = 0; j != peak_map[i].size(); ++j)
      {
        const double mz = peak_map[i][j].getMZ();
        const double intensity = peak_map[i][j].getIntensity();  

        // keep all peaks higher than upper quartile
        if (intensity > median + rMAD) 
        {
          iso_ints.push_back(intensity);
          indices.push_back(j);
          continue;
        }
 
        // keep all peaks that have some isotopic neighbors
        double tol = Math::ppmToMass(ppm, mz);

        for (vector<double>::const_iterator it = mz_shifts.begin(); it != mz_shifts.end(); ++it)
        {
          // isotopic peak to the left (or right)?
          if (peak_map[i].findNearest(mz - *it, tol) != -1)
          {
            // check if at least one above or below is present
            if ( peak_map[i + 1].findNearest(mz - *it, tol) != -1
              || peak_map[i - 1].findNearest(mz - *it, tol) != -1
              || peak_map[i + 2].findNearest(mz - *it, tol) != -1
              || peak_map[i - 2].findNearest(mz - *it, tol) != -1
               )
            {
             iso_ints.push_back(intensity);
             indices.push_back(j);
             break;
            }         
          }
          else if (peak_map[i].findNearest(mz + *it, tol) != -1)
          {
            if ( peak_map[i + 1].findNearest(mz + *it, tol) != -1
              || peak_map[i - 1].findNearest(mz + *it, tol) != -1
              || peak_map[i + 2].findNearest(mz + *it, tol) != -1
              || peak_map[i - 2].findNearest(mz + *it, tol) != -1
               )
              {
               iso_ints.push_back(intensity);
               indices.push_back(j);
               break;
              }           
          }
        } 
      }
      peak_map[i].select(indices); // filter all without isotopic peak
    }
    calculateQuantiles_(iso_ints, "MS run isotopic intensities: ", median, t2, rMAD);
/* TODO: IFDEF
    MzMLFile().store("filtered.mzML", peak_map);
*/
  }

  static MassTracesStatistics_ getMassTraceStatistics_(const FeatureMap & feature_map)
  {
    vector<double> stats_sd, stats_sd_ppm, stats_fwhm, stats_mz_wiggle, stats_int_wiggle, stats_fwhm_high_int_frac, stats_fwhm_ints;

    for (FeatureMap::ConstIterator fm_it = feature_map.begin(); fm_it != feature_map.end(); ++fm_it)
    {
      stats_sd.push_back((double)fm_it->getMetaValue("SD"));
      stats_sd_ppm.push_back((double)fm_it->getMetaValue("SD_ppm"));
      stats_fwhm.push_back((double)fm_it->getWidth());
      stats_mz_wiggle.push_back((double)fm_it->getMetaValue("mz_wiggle"));
      stats_int_wiggle.push_back((double)fm_it->getMetaValue("int_wiggle"));
      stats_fwhm_high_int_frac.push_back((double)fm_it->getMetaValue("fwhm_high_int_frac"));
      stats_fwhm_ints.push_back((double)fm_it->getMetaValue("fwhm_int_left"));
      stats_fwhm_ints.push_back((double)fm_it->getMetaValue("fwhm_int_right"));
    }

    MassTracesStatistics_ stats;

    double t1,t2,t3;
    calculateQuantiles_(stats_sd, "Mass trace sd m/z: ", t1, t2, t3);
    calculateQuantiles_(stats_sd_ppm, "Mass trace sd ppm: ", stats.sd_ppm_median, stats.sd_ppm_lMAD, stats.sd_ppm_rMAD);
    calculateQuantiles_(stats_fwhm, "Mass trace FWHM: ", t1, t2, t3);
    calculateQuantiles_(stats_mz_wiggle, "Mass trace m/z wiggle in fwhm: ", stats.mz_wiggle_median, stats.mz_wiggle_lMAD, stats.mz_wiggle_rMAD);
    calculateQuantiles_(stats_int_wiggle, "Mass trace int wiggle in fwhm: ", t1, t2, t3);
    calculateQuantiles_(stats_fwhm_high_int_frac, "High intensity peak fraction outside of FWHM: ", t1, t2, t3);
    calculateQuantiles_(stats_fwhm_ints, "Mass trace FWHM border intensities: ", t1, t2, t3);

    return stats;
  }

  static void filterMassTracesByStatistics_(FeatureMap & feature_map, const MassTracesStatistics_& stats)
  {
    //copy all properties
    FeatureMap map_sm = feature_map;
    FeatureMap map_removed = feature_map;
    map_sm.clear(false);
    map_removed.clear(false); //TODO debug
    Size filtered(0);
    for (FeatureMap::Iterator fm_it = feature_map.begin(); fm_it != feature_map.end(); ++fm_it)
    {
      if ((double)fm_it->getMetaValue("SD_ppm") <= stats.sd_ppm_median + 5.0 * stats.sd_ppm_rMAD &&
          (double)fm_it->getMetaValue("mz_wiggle") <= stats.mz_wiggle_median + 5.0 * stats.mz_wiggle_rMAD &&
          (double)fm_it->getMetaValue("int_wiggle") < 1.0 &&
          (double)fm_it->getMetaValue("fwhm_high_int_frac") < 0.5)
      {
        map_sm.push_back(*fm_it);
      }
      else
      {
        map_removed.push_back(*fm_it);
        ++filtered;
      }
    }
    LOG_INFO << "Detected and filtered: " << filtered << " outliers (> 5MAD)." << std::endl;
    LOG_INFO << "Mass traces: " << map_sm.size() << std::endl;
    feature_map = map_sm;
    feature_map.updateRanges();
    map_removed.updateRanges();
    FeatureXMLFile().store("removed.featureXML", map_removed);
  }

  void estimateTraceParameters(MSExperiment<> ms_peakmap)
  {    
    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    bool use_epd = epd_param.getValue("enabled").toBool();
    double start_fwhm = (double) com_param.getValue("chrom_fwhm");

    //-------------------------------------------------------------
    // configure and run MTD
    //-------------------------------------------------------------
    MassTraceDetection mt_ext;
    mtd_param.insert("", com_param);
    mtd_param.remove("chrom_fwhm");
    mt_ext.setParameters(mtd_param);

    double start_sd_ppm = (double)mtd_param.getValue("mass_error_ppm");

    // filter peakmap todo: rename
    getPeakMapStatistics_(ms_peakmap, start_sd_ppm);

    vector<MassTrace> m_traces;
    mt_ext.run(ms_peakmap, m_traces, 1000);

    // 1a. run: estimate ppm error
    double median_sd_ppm, median_fwhm;
    double t2,t3;

    if (use_epd)
    {
      ElutionPeakDetection ep_det;

      Param p_epd = ElutionPeakDetection().getDefaults();
      p_epd.setValue("max_fwhm", 1e10);
      p_epd.setValue("chrom_fwhm", start_fwhm);
      ep_det.setParameters(p_epd);

      std::vector<MassTrace> split_mtraces;
      ep_det.detectPeaks(m_traces, split_mtraces);
      swap(split_mtraces, m_traces);
    }

    LOG_INFO << "Using " << m_traces.size() << " mass traces for parameter estimation." << endl;  

    std::vector<double> stats_sd, stats_sd_ppm, stats_fwhm;

    for (Size i = 0; i < m_traces.size(); ++i)
    {
      MassTrace & m = m_traces[i];

      if (m.getSize() == 0) continue;

      m.updateMeanMZ();
      m.updateWeightedMZsd();

      // SD absolute and ppm on full trace
      double sd = m.getCentroidSD();
      double sd_ppm = sd / m.getCentroidMZ() * 1e6;
      stats_sd.push_back(sd);
      stats_sd_ppm.push_back(sd_ppm);

      // FWHM
      double fwhm = m.estimateFWHM(use_epd);
      stats_fwhm.push_back(fwhm);
    }

    calculateQuantiles_(stats_sd_ppm, "Mass trace sd ppm: ", median_sd_ppm, t2, t3);
    calculateQuantiles_(stats_fwhm, "Mass trace FWHM: ", median_fwhm, t2, t3);

    // 1b. improve ppm estimate by narrowing tolerance window
    if (median_sd_ppm < start_sd_ppm) 
    {
      // limit median to 1 ppb (to cope with simulated data)
      if (median_sd_ppm < 1e-3) median_sd_ppm = 1e-3;

      // filter peakmap todo: rename
      getPeakMapStatistics_(ms_peakmap, median_sd_ppm);

      stats_sd.clear();
      stats_sd_ppm.clear();
      stats_fwhm.clear();
      m_traces.clear();

      mtd_param.setValue("mass_error_ppm", median_sd_ppm);
      mt_ext.setParameters(mtd_param);
      mt_ext.run(ms_peakmap, m_traces, 1000);

      if (use_epd)
      {
        ElutionPeakDetection ep_det;

        Param p_epd = ElutionPeakDetection().getDefaults();
        p_epd.setValue("max_fwhm", 1e10);
        p_epd.setValue("chrom_fwhm", median_fwhm);
        ep_det.setParameters(p_epd);

        std::vector<MassTrace> split_mtraces;
        ep_det.detectPeaks(m_traces, split_mtraces);
        swap(split_mtraces, m_traces);
      }

      std::vector<double> stats_sd, stats_sd_ppm, stats_fwhm;

      for (Size i = 0; i < m_traces.size(); ++i)
      {
        MassTrace & m = m_traces[i];

        if (m.getSize() == 0) continue;

        m.updateMeanMZ();
        m.updateWeightedMZsd();

        // SD absolute and ppm on full trace
        double sd = m.getCentroidSD();
        double sd_ppm = sd / m.getCentroidMZ() * 1e6;
        stats_sd.push_back(sd);
        stats_sd_ppm.push_back(sd_ppm);

        // FWHM
        double fwhm = m.estimateFWHM(use_epd);
        stats_fwhm.push_back(fwhm);
      }

      // print some stats about standard deviation of mass traces
      calculateQuantiles_(stats_sd_ppm, "Mass trace sd ppm: ", median_sd_ppm, t2, t3);
      calculateQuantiles_(stats_fwhm, "Mass trace FWHM: ", median_fwhm, t2, t3);
    }

    LOG_INFO << "Paramter estimated (sd ppm / chrom FWHM): " << median_sd_ppm << " / " << median_fwhm << endl;
  }
 
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input centroided mzML file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output featureXML file with mass traces");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));
    registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content", false);
    setValidStrings_("out_type", ListUtils::create<String>("featureXML,consensusXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");

  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    Param combined;
    Param p_com;
    p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
    p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");

    combined.insert("common:", p_com);

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.remove("noise_threshold_int");
    p_mtd.remove("chrom_peak_snr");

    combined.insert("mtd:", p_mtd);

    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.remove("noise_threshold_int");
    p_epd.remove("chrom_peak_snr");
    p_epd.remove("chrom_fwhm");

    p_epd.setValue("enabled", "true", "Enables/disables the chromatographic peak detection of mass traces");
    p_epd.setValidStrings("enabled", ListUtils::create<String>("true,false"));
    combined.insert("epd:", p_epd);

    return combined;
  }

  ExitCodes main_(int, const char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    MSExperiment<Peak1D> ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    (mz_data_file.getOptions()).setMSLevels(ms_level);
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.size() == 0)
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    
    // make sure that the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    //-------------------------------------------------------------
    // get params for MTD and EPD algorithms
    //-------------------------------------------------------------
    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);


    //-------------------------------------------------------------
    // configure and run MTD
    //-------------------------------------------------------------
    MassTraceDetection mt_ext;
    mtd_param.insert("", com_param);
    mtd_param.remove("chrom_fwhm");
    mt_ext.setParameters(mtd_param);

    estimateTraceParameters(ms_peakmap);

    vector<MassTrace> m_traces;
    mt_ext.run(ms_peakmap, m_traces);

    vector<MassTrace> m_traces_final;

    bool use_epd = epd_param.getValue("enabled").toBool();

    if (!use_epd)
    {
      swap(m_traces_final, m_traces);
    }
    else
    {
      ElutionPeakDetection ep_det;

      epd_param.remove("enabled"); // artificially added above
      epd_param.insert("", com_param);

      ep_det.setParameters(epd_param);

      std::vector<MassTrace> split_mtraces;
      // note: this step will destroy any meta data annotation (e.g. FWHM_mz_avg)
      ep_det.detectPeaks(m_traces, split_mtraces);

      if (ep_det.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        ep_det.filterByPeakWidth(split_mtraces, m_traces_final);

        LOG_INFO << "Notice: " << split_mtraces.size() - m_traces_final.size()
                 << " of total " << split_mtraces.size() 
                 << " were dropped because of too low peak width." << std::endl;
      }
      else
      {
        swap(m_traces_final, split_mtraces);
      }

    }

    //-------------------------------------------------------------
    // writing consensus map output
    //-------------------------------------------------------------
    if (out_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap consensus_map;
      consensus_map.setPrimaryMSRunPath(ms_peakmap.getPrimaryMSRunPath());
      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        if (m_traces_final[i].getSize() == 0) continue;

        ConsensusFeature fcons;
        int k = 0;
        for (MassTrace::const_iterator it = m_traces_final[i].begin(); it != m_traces_final[i].end(); ++it)
        {
          FeatureHandle fhandle;
          fhandle.setRT(it->getRT());
          fhandle.setMZ(it->getMZ());
          fhandle.setIntensity(it->getIntensity());
          fhandle.setUniqueId(++k);
          fcons.insert(fhandle);
        }

        fcons.setMetaValue(3, m_traces_final[i].getLabel());
        fcons.setCharge(0);
        fcons.setWidth(m_traces_final[i].estimateFWHM(use_epd));
        fcons.setQuality(1 - (1.0 / m_traces_final[i].getSize()));

        fcons.setRT(m_traces_final[i].getCentroidRT());
        fcons.setMZ(m_traces_final[i].getCentroidMZ());
        fcons.setIntensity(m_traces_final[i].getIntensity(false));
        consensus_map.push_back(fcons);
      }
      consensus_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      addDataProcessing_(consensus_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      consensus_map.setUniqueId();
      ConsensusXMLFile().store(out, consensus_map);

    }
    else //(out_type == FileTypes::FEATUREXML)
    {

      //-----------------------------------------------------------
      // convert mass traces to features
      //-----------------------------------------------------------

      std::vector<double> stats_sd, stats_sd_ppm, stats_fwhm, stats_mz_wiggle, stats_int_wiggle, stats_fwhm_int;
      FeatureMap ms_feat_map;
      ms_feat_map.setPrimaryMSRunPath(ms_peakmap.getPrimaryMSRunPath());

      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        MassTrace & m = m_traces_final[i];

        if (m.getSize() == 0) continue;

        m.updateMeanMZ();
        m.updateWeightedMZsd();

        Feature f;
        f.setMetaValue(3, m.getLabel());
        f.setCharge(0);
        f.setMZ(m.getCentroidMZ());
        f.setIntensity(m.getIntensity(false));
        f.setRT(m.getCentroidRT());
        f.setWidth(m.estimateFWHM(use_epd));
        f.setOverallQuality(1 - (1.0 / m.getSize()));
        f.getConvexHulls().push_back(m.getConvexhull());
        stats_fwhm.push_back(f.getWidth());

        // SD absolute and ppm on full trace
        double sd = m.getCentroidSD();
        double sd_ppm = sd / f.getMZ() * 1e6;
        stats_sd.push_back(sd);
        stats_sd_ppm.push_back(sd_ppm);
        f.setMetaValue("SD", sd);        
        f.setMetaValue("SD_ppm", sd_ppm);

        // stdev. of the m/z and int differences in FWHM region (wiggle)
        std::pair<Size, Size> fwhm_idx = m.getFWHMborders();
        stats_fwhm_int.push_back(m[fwhm_idx.first].getIntensity());
        stats_fwhm_int.push_back(m[fwhm_idx.second].getIntensity());
        f.setMetaValue("fwhm_int_left", m[fwhm_idx.first].getIntensity());
        f.setMetaValue("fwhm_int_right", m[fwhm_idx.second].getIntensity());
        vector<double> mz_wiggle;
        double int_wiggle(0);
        double mz_before = m[fwhm_idx.first].getMZ();
        double int_before = m[fwhm_idx.first].getIntensity();
        double sign_before(0);
        for (Size k = fwhm_idx.first + 1; k <= fwhm_idx.second; ++k)
        {
          double diff_ppm = fabs(mz_before - m[k].getMZ()) / m.getCentroidMZ() * 1e6;
          mz_wiggle.push_back(diff_ppm);

          double diff_int = m[k].getIntensity() - int_before;
          if (diff_int > 0 && sign_before < 0)  // local minimum?
          {
            int_wiggle += 1; // count minimum
            sign_before = 1;
          }
          else if (diff_int < 0 && sign_before > 0) // local maximum
          {
            int_wiggle += 1; // count maximum
            sign_before = -1;
          }
          else // either same slope, last peak had equal intensity, or first peak in fwhm
          {
            if (diff_int > 0)
            {
              sign_before = 1;
            }
            else if (diff_int < 0)
            {
              sign_before = -1;
            }
            else
            {
              sign_before = 0;
            }
          }

          mz_before = m[k].getMZ();
          int_before = m[k].getIntensity();
        }
        double sd_mz_wiggle(0);
        if (mz_wiggle.size() > 1) sd_mz_wiggle = Math::sd(mz_wiggle.begin(), mz_wiggle.end());
        stats_mz_wiggle.push_back(sd_mz_wiggle);

        if (fwhm_idx.second - fwhm_idx.first >= 2) 
        {
          int_wiggle /= (fwhm_idx.second - fwhm_idx.first - 1.0);
        }
        stats_int_wiggle.push_back(int_wiggle);

        f.setMetaValue("mz_wiggle", sd_mz_wiggle);
        f.setMetaValue("int_wiggle", int_wiggle);

        // determine fraction of high intensity peaks outside of FWHM region 
        double half_max_intensity = m.getMaxIntensity(false) / 2.0;
        Size high_int_outside(0);
        Size count(0);
        for (Size k = 0; k < fwhm_idx.first; ++k)
        {
          if (m[k].getIntensity() > half_max_intensity) {++high_int_outside;}
          count++;
        } 

        for (Size k = fwhm_idx.second + 1; k < m.getSize(); ++k)
        {
          if (m[k].getIntensity() > half_max_intensity) {++high_int_outside;}
          count++;
        } 

        double high_frac_outside = count > 0 ? (double) high_int_outside / (double)count : high_int_outside;

        f.setMetaValue("fwhm_high_int_frac", high_frac_outside);

        if (m.fwhm_mz_avg > 0) f.setMetaValue("FWHM_mz_avg", m.fwhm_mz_avg);
        
        ms_feat_map.push_back(f);
      }

      LOG_INFO << "Mass traces: " << m_traces_final.size() << "\n";

      // print some stats about standard deviation of mass traces
      MassTracesStatistics_ stats = getMassTraceStatistics_(ms_feat_map);

      // filter outliers
      filterMassTracesByStatistics_(ms_feat_map, stats);

      // print some stats about standard deviation of mass traces
      stats = getMassTraceStatistics_(ms_feat_map);

      ms_feat_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info TODO
      addDataProcessing_(ms_feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      //ms_feat_map.setUniqueId();

      FeatureXMLFile().store(out, ms_feat_map);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMassTraceExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
