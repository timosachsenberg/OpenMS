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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/MorphologicalIsotopeFilter.h>

using namespace std;

namespace OpenMS
{

  MorphologicalIsotopeFilter::MorphologicalIsotopeFilter() :
    DefaultParamHandler("MorphologicalIsotopeFilter")
  {
    init_();
  }

  MorphologicalIsotopeFilter::MorphologicalIsotopeFilter(double tolerance) :
    DefaultParamHandler("MorphologicalIsotopeFilter")
  {
    init_();
    // after initialising with the default value, use the provided n
    param_.setValue("tolerance", tolerance);
    updateMembers_();
  }

  void MorphologicalIsotopeFilter::init_()
  {
    defaults_.setValue("tolerance", 10.0, "The m/z tolerance to isotopic peaks (ppm).");
    defaultsToParam_();
  }

  MorphologicalIsotopeFilter::~MorphologicalIsotopeFilter()
  {
  }

  MorphologicalIsotopeFilter::MorphologicalIsotopeFilter(const MorphologicalIsotopeFilter & source) :
    DefaultParamHandler(source)
  {
    updateMembers_();
  }

  MorphologicalIsotopeFilter & MorphologicalIsotopeFilter::operator=(const MorphologicalIsotopeFilter & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
      updateMembers_();
    }
    return *this;
  }

  void MorphologicalIsotopeFilter::filterPeakMap(PeakMap & peak_map)
  {
    // gather statistics on intensity distribution
    vector<double> ints;
    for (Size i = 0; i != peak_map.size(); ++i)
    {
      for (Size j = 0; j != peak_map[i].size(); ++j)
      {
        ints.push_back(peak_map[i][j].getIntensity());
      }
    }

    // calculate median and rMAD of intensities
    double median, rMAD;
    std::sort(ints.begin(), ints.end());
    median = ints[ints.size() * 5 / 10];
    vector<double>::iterator int_mid = ints.begin() + ints.size() * 0.5;
    rMAD = Math::MAD(int_mid, ints.end(), median);

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
        double tol = Math::ppmToMass(tol_, mz);

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
  }

  void MorphologicalIsotopeFilter::updateMembers_()
  {
    tol_ = (UInt)param_.getValue("tolerance");
  }

}
