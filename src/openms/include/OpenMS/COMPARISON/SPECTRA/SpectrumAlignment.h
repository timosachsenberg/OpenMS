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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <iostream>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#define ALIGNMENT_DEBUG
#undef  ALIGNMENT_DEBUG

namespace OpenMS
{

  /**
      @brief Aligns the peaks of two sorted spectra
      Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given.
                Scoring function is the m/z distance between peaks. Intensity does not play a role!

      Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
      Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window
      @note: a peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2.
      @note: intensity is ignored TODO: improve time complexity, currently O(|s1|*log(|s2|))

      @htmlinclude OpenMS_SpectrumAlignment.parameters

      @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI SpectrumAlignment :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectrumAlignment();

    /// copy constructor
    SpectrumAlignment(const SpectrumAlignment & source);

    /// destructor
    ~SpectrumAlignment() override;

    /// assignment operator
    SpectrumAlignment & operator=(const SpectrumAlignment & source);
    // @}

    template <typename SpectrumType1, typename SpectrumType2>
    void getSpectrumAlignmentFastCharge(
      std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance, 
      bool fragment_mass_tolerance_unit_ppm, 
      const SpectrumType1& theo_spectrum, 
      const SpectrumType2& exp_spectrum,
      const typename SpectrumType1::IntegerDataArray& theo_charges,
      const typename SpectrumType2::IntegerDataArray& exp_charges,
      const double min_ratio = 0)
    {
      OPENMS_PRECONDITION(exp_spectrum.isSorted(), "Spectrum needs to be sorted.");
      OPENMS_PRECONDITION(theo_spectrum.isSorted(), "Spectrum needs to be sorted.");
      OPENMS_PRECONDITION((alignment.empty() == true), "Alignment result vector needs to be empty.");

      const Size n_t(theo_spectrum.size());
      const Size n_e(exp_spectrum.size());
      const bool has_charge = !(exp_charges.empty() || theo_charges.empty());
                      
      if (n_t == 0 || n_e == 0) { return; }

      Size t(0), e(0);
      alignment.reserve(theo_spectrum.size());

      while (t < n_t && e < n_e)
      {
        const double theo_mz = theo_spectrum[t].getMZ();
        const double exp_mz = exp_spectrum[e].getMZ();

        double tz(0), ez(0); 
        if (has_charge)  
        {
          tz = theo_charges[t];
          ez = exp_charges[e];
        }
        const bool tz_matches_ez = (ez == tz || !ez || !tz);

        const double ratio = min_ratio > 0 ? theo_spectrum[t].getIntensity() / exp_spectrum[e].getIntensity() : 1.0;
        const bool ti_matches_ei = (ratio >= min_ratio || ratio >= 1.0/min_ratio);

        double d = exp_mz - theo_mz;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        if (fabs(d) <= max_dist_dalton) // match in tolerance window? 
        {
          // get first peak with matching charge in tolerance window
          if (!tz_matches_ez || !ti_matches_ei)
          {
            Size e_candidate(e);
            while (e_candidate < n_e-1)
            {
              ++e_candidate;
              double new_ez = has_charge ? exp_charges[e_candidate] : 0;
              const bool charge_matches = (new_ez == tz || !new_ez || !tz);
 
              double new_ratio = min_ratio > 0 ? theo_spectrum[t].getIntensity() / exp_spectrum[e_candidate].getIntensity() : 1.0;
              const bool ratio_matches = (new_ratio >= min_ratio || new_ratio >= 1.0/min_ratio);;

              double new_d = exp_spectrum[e].getMZ() - theo_mz;
              if (charge_matches && ratio_matches && new_d <= max_dist_dalton)
              { // found a match
                break;
              }
              else if (new_d > max_dist_dalton)
              { // no match found              
                e_candidate = e;                
                break;
              }
            }
            if (e == e_candidate) 
            { // no match found continue with next theo. peak
              ++t;
              continue;
            }
            else
            { // match found
              e = e_candidate;
            }
          }

          // Invariant: e now points to the first peak in tolerance window, that matches in charge

          // last peak? there can't be a better one in this tolerance window
          if (e >= n_e - 1) { alignment.emplace_back(std::make_pair(t, e)); return; }

          Size closest_exp_peak(e); 

          // Invariant: closest_exp_peak always point to best match 

          double new_d, new_ez(0);
          double best_d = exp_spectrum[closest_exp_peak].getMZ() - theo_mz;            

//        std::cerr << "first peak in window:" <<  exp_spectrum[closest_exp_peak].getMZ() << "\t theo: " << theo_mz << "\t best d:" << best_d << "\n";

          do // check for better match in tolerance window
          {
            // advance to next exp. peak
            ++e;

            // determine distance of next peak
            new_d = exp_spectrum[e].getMZ() - theo_mz;
            const bool in_tolerance_window = (fabs(new_d) < max_dist_dalton);
//          std::cerr << "  new peak:" << exp_spectrum[e].getMZ() << "\t theo: " << theo_mz << "\t new d: " << new_d << "\t in window:" << in_tolerance_window << "\n";

            if (!in_tolerance_window) { break; }

            // Invariant: e is in tolerance window

            // check if charge of next peak matches
            if (has_charge) { new_ez = exp_charges[e]; }
            const bool charge_matches = (new_ez == tz || !new_ez || !tz);
            // check if ratio of next peak matches
            const double new_ratio = min_ratio > 0 ? theo_spectrum[t].getIntensity() / exp_spectrum[e].getIntensity() : 1.0;
            const bool ratio_matches = (new_ratio >= min_ratio || new_ratio >= 1.0/min_ratio);
            if (!charge_matches || !ratio_matches) { continue; }

            // Invariant: charge and ratio matches

            const bool better_distance = (fabs(new_d) <= fabs(best_d));

            // better distance (and matching charge)? better match found
            if (better_distance) 
            { // found a better match
              closest_exp_peak = e;
              best_d = new_d;
            }
            else 
            { // distance got worse -> no additional matches!
              break; 
            }
          } 
          while (e < n_e - 1);
         
//        std::cerr << "added peak:" << exp_spectrum[closest_exp_peak].getMZ() << "\t theo: " << theo_mz << "\n";

          // search in tolerance window for an experimental peak closer to theoretical one 
          alignment.emplace_back(std::make_pair(t, closest_exp_peak));
     
          e = closest_exp_peak + 1;  // advance experimental peak to 1-after the best match
          ++t; // advance theoretical peak
        }
        else if (d < 0) // exp. peak is left of theo. peak (outside of tolerance window)
        {
          ++e; 
        }
        else if (d > 0) // theo. peak is left of exp. peak (outside of tolerance window)
        {
          ++t;
        }
      }
    }

    template <typename SpectrumType1, typename SpectrumType2>
    void getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const SpectrumType1 & s1, const SpectrumType2 & s2) const
    {
      if (!s1.isSorted() || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
      }

      // clear result
      alignment.clear();
      double tolerance = (double)param_.getValue("tolerance");

      if (!param_.getValue("is_relative_tolerance").toBool() )
      {
        std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
        std::map<Size, std::map<Size, double> > matrix;

        // init the matrix with "gap costs" tolerance
        matrix[0][0] = 0;
        for (Size i = 1; i <= s1.size(); ++i)
        {
          matrix[i][0] = i * tolerance;
          traceback[i][0]  = std::make_pair(i - 1, 0);
        }
        for (Size j = 1; j <= s2.size(); ++j)
        {
          matrix[0][j] = j * tolerance;
          traceback[0][j] = std::make_pair(0, j - 1);
        }

        // fill in the matrix
        Size left_ptr(1);
        Size last_i(0), last_j(0);

        //Size off_band_counter(0);
        for (Size i = 1; i <= s1.size(); ++i)
        {
          double pos1(s1[i - 1].getMZ());

          for (Size j = left_ptr; j <= s2.size(); ++j)
          {
            bool off_band(false);
            // find min of the three possible directions
            double pos2(s2[j - 1].getMZ());
            double diff_align = fabs(pos1 - pos2);

            // running off the right border of the band?
            if (pos2 > pos1 && diff_align > tolerance)
            {
              if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
              {
                off_band = true;
              }
            }

            // can we tighten the left border of the band?
            if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
            {
              ++left_ptr;
            }

            double score_align = diff_align;

            if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j - 1) != matrix[i - 1].end())
            {
              score_align += matrix[i - 1][j - 1];
            }
            else
            {
              score_align += (i - 1 + j - 1) * tolerance;
            }

            double score_up = tolerance;
            if (matrix.find(i) != matrix.end() && matrix[i].find(j - 1) != matrix[i].end())
            {
              score_up += matrix[i][j - 1];
            }
            else
            {
              score_up += (i + j - 1) * tolerance;
            }

            double score_left = tolerance;
            if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j) != matrix[i - 1].end())
            {
              score_left += matrix[i - 1][j];
            }
            else
            {
              score_left += (i - 1 + j) * tolerance;
            }

    #ifdef ALIGNMENT_DEBUG
          cerr << i << " " << j << " " << left_ptr << " " << pos1 << " " << pos2 << " " << score_align << " " << score_left << " " << score_up << endl;
    #endif

          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
             matrix[i][j] = score_align;
             traceback[i][j] = std::make_pair(i - 1, j - 1);
             last_i = i;
             last_j = j;
          }
          else
          {
            if (score_up <= score_left)
            {
              matrix[i][j] = score_up;
              traceback[i][j] = std::make_pair(i, j - 1);
            }
            else
            {
              matrix[i][j] = score_left;
              traceback[i][j] = std::make_pair(i - 1, j);
            }
          }

          if (off_band)
          {
            break;
          }
        }
      }

          //last_i = s1.size() + 1;
          //last_j = s2.size() + 1;

          //cerr << last_i << " " << last_j << endl;

    #ifdef ALIGNMENT_DEBUG
    #if 0
          cerr << "TheMatrix: " << endl << " \t  \t";
          for (Size j = 0; j != s2.size(); ++j)
          {
        cerr << s2[j].getPosition()[0] << " \t";
          }
          cerr << endl;
          for (Size i = 0; i <= s1.size(); ++i)
          {
        if (i != 0)
        {
          cerr << s1[i - 1].getPosition()[0] << " \t";
        }
        else
        {
          cerr << " \t";
        }
        for (Size j = 0; j <= s2.size(); ++j)
        {
          if (matrix.has(i) && matrix[i].has(j))
          {
            if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
            {
              cerr << "\\";
            }
            else
            {
              if (traceback[i][j].first == i - 1 && traceback[i][j].second == j)
              {
            cerr << "|";
              }
              else
              {
            cerr << "-";
              }
            }

            cerr << matrix[i][j] << "  \t";
          }
          else
          {
            cerr << "-1  \t";
          }
        }
        cerr << endl;
          }
    #endif
    #endif

      // do traceback
      Size i = last_i;
      Size j = last_j;

      while (i >= 1 && j >= 1)
      {
        if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
        {
          alignment.push_back(std::make_pair(i - 1, j - 1));
        }
        Size new_i = traceback[i][j].first;
        Size new_j = traceback[i][j].second;

        i = new_i;
        j = new_j;
      }

      std::reverse(alignment.begin(), alignment.end());

    #ifdef ALIGNMENT_DEBUG
    #if 0
          // print alignment
          cerr << "Alignment (size=" << alignment.size() << "): " << endl;

          Size i_s1(0), i_s2(0);
          for (vector<pair<Size, Size> >::const_reverse_iterator it = alignment.rbegin(); it != alignment.rend(); ++it, ++i_s1, ++i_s2)
          {
        while (i_s1 < it->first - 1)
        {
          cerr << i_s1 << " " << s1[i_s1].getPosition()[0] << " " << s1[i_s1].getIntensity() << endl;
          i_s1++;
        }
        while (i_s2 < it->second - 1)
        {
          cerr << " \t " <<  i_s2 << " " << s2[i_s2].getPosition()[0] << " " << s2[i_s2].getIntensity() << endl;
          i_s2++;
        }
        cerr << "(" << s1[it->first - 1].getPosition()[0] << " <-> " << s2[it->second - 1].getPosition()[0] << ") ("
             << it->first << "|" << it->second << ") (" << s1[it->first - 1].getIntensity() << "|" << s2[it->second - 1].getIntensity() << ")" << endl;
          }
    #endif
    #endif
    }
    else  // relative alignment (ppm tolerance)
    {        
      for (Size i = 0; i != s1.size(); ++i)
      {
        const double& theo_mz = s1[i].getMZ();
        double max_dist_dalton = theo_mz * tolerance * 1e-6;
 
        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size j = s2.findNearest(theo_mz);
        double exp_mz = s2[j].getMZ();

        // found peak match
        if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
        {
          alignment.push_back(std::make_pair(i, j));
        }
      }
    }
  }
 };
}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
