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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_RIBONUCLEOTIDEDB_H
#define OPENMS_CHEMISTRY_RIBONUCLEOTIDEDB_H

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <boost/unordered_map.hpp>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Database of ribonucleotides (modified and unmodified)

      The information in this class comes from the Modomics database (http://modomics.genesilico.pl/modifications/) and is read from a tab-separated text file in @p data/CHEMISTRY/Modomics.tsv.
  */
  class OPENMS_DLLAPI RibonucleotideDB
  {
  public:
    /// const iterator type definition
    typedef std::vector<Ribonucleotide>::const_iterator ConstIterator;

    /// replacement for constructor (singleton pattern)
    inline static RibonucleotideDB* getInstance()
    {
      static RibonucleotideDB* db_ = 0;
      if (db_ == 0)
      {
        db_ = new RibonucleotideDB;
      }
      return db_;
    }

    /// destructor
    virtual ~RibonucleotideDB();

    /// copy constructor not available
    RibonucleotideDB(const RibonucleotideDB& other) = delete;

    /// assignment operator not available
    RibonucleotideDB& operator=(const RibonucleotideDB& other) = delete;

    /// Const iterator to beginning of database
    inline ConstIterator begin() const
    {
      return ribonucleotides_.begin();
    }

    /// Const iterator to end of database
    inline ConstIterator end() const
    {
      return ribonucleotides_.end();
    }

    /**
       @brief Get a ribonucleotide by its code (short name)

       @throw Exception::ElementNotFound if nothing was found
    */
    const Ribonucleotide& getRibonucleotide(const String& code);

    /**
       @brief Get the ribonucleotide with the longest code that matches a prefix of @p seq

       @throw Exception::ElementNotFound if nothing was found
    */
    const Ribonucleotide& getRibonucleotidePrefix(const String& seq);

  protected:
    /// default constructor
    RibonucleotideDB();

    /// read (modified) nucleotides from input file
    void readFromFile_(const String& path);

    /// create a (modified) nucleotide from an input row
    Ribonucleotide parseRow_(const String& row, Size line_count);

    /// list of known (modified) nucleotides
    std::vector<Ribonucleotide> ribonucleotides_;

    /// mapping of codes (short names) to indexes into @p ribonucleotides_
    boost::unordered_map<String, Size> code_map_;

    Size max_code_length_;
  };
}

#endif