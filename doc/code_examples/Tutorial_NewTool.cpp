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
// $Maintainer: FIRST LAST $
// $Authors: FIRST LAST
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MyPeptideSearchEngine MyPeptideSearchEngine

    @brief Perfroms simple peptide identification.

    This application is for demonstratin purposes only.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MyPeptideSearchEngine.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MyPeptideSearchEngine.html
*/

// We do not want the class to show up in the documentation (only the tool):
/// @cond TOPPCLASSES

class TOPPMyPeptideSearchEngine :
  public TOPPBase
{
public:

  // Note that name provided as first argument should match the cpp file.
  // The description should be short and descriptive as it will be displayed on the command line
  // as well as workflow systems.
  // The last boolean parameter marks it as unstable (=util)
  TOPPMyPeptideSearchEngine() :
    TOPPBase("MyPeptideSearchEngine", "Perform a simple peptide identification task.", false)
  {
  }

protected:

  ExitCodes main_(int, const char **)
  {
    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPMyPeptideSearchEngine tool;
  return tool.main(argc, argv);
}

/// @endcond
