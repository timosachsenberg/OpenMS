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
// $Authors: Matthias Fuchs $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/FeatureSQLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <string>
#include <cassert>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(FeatureMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureMap* pl_ptr = nullptr;
FeatureMap* nullPointer = nullptr;


START_SECTION((FeatureMap()))
	pl_ptr = new FeatureMap();
  TEST_NOT_EQUAL(pl_ptr, nullPointer)

	TEST_EQUAL(pl_ptr->getMin(), FeatureMap::PositionType::maxPositive())
	TEST_EQUAL(pl_ptr->getMax(), FeatureMap::PositionType::minNegative())
	TEST_REAL_SIMILAR(pl_ptr->getMinInt(), numeric_limits<double>::max())
	TEST_REAL_SIMILAR(pl_ptr->getMaxInt(), -numeric_limits<double>::max())
END_SECTION


/*
START_SECTION((FeatureMap(const FeatureMap &source)))
  Feature feature1;
  feature1.getPosition()[0] = 2.0;
  feature1.getPosition()[1] = 3.0;
  feature1.setIntensity(1.0f);

  Feature feature2;
  feature2.getPosition()[0] = 0.0;
  feature2.getPosition()[1] = 2.5;
  feature2.setIntensity(0.5f);

  Feature feature3;
  feature3.getPosition()[0] = 10.5;
  feature3.getPosition()[1] = 0.0;
  feature3.setIntensity(0.01f);
  
  FeatureMap map1;
  map1.setMetaValue("meta",String("value"));
  map1.push_back(feature1);
  map1.push_back(feature2);
  map1.push_back(feature3);
  map1.updateRanges();
  map1.setIdentifier("lsid");;
  map1.getDataProcessing().resize(1);


  // Calculate and output the ranges
  map1.updateRanges();
  cout << "Int: " << map1.getMinInt() << " - " << map1.getMaxInt() << endl;
  cout << "RT:  " << map1.getMin()[0] << " - " << map1.getMax()[0] << endl;
  cout << "m/z: " << map1.getMin()[1] << " - " << map1.getMax()[1] << endl;
  
END_SECTION
*/




START_SECTION((void write(std::string& out_fm, FeatureMap& feature_map)))
{
  FeatureXMLFile dfmap_file;
  
  Feature feature1;
  feature1.getPosition()[0] = 2.0;
  feature1.getPosition()[1] = 3.0;
  feature1.setIntensity(1.0f);

  Feature feature2;
  feature2.getPosition()[0] = 0.0;
  feature2.getPosition()[1] = 2.5;
  feature2.setIntensity(0.5f);

  Feature feature3;
  feature3.getPosition()[0] = 10.5;
  feature3.getPosition()[1] = 0.0;
  feature3.setIntensity(0.01f);
  
  FeatureMap map1;
  map1.setMetaValue("meta",String("value"));
  map1.push_back(feature1);
  map1.push_back(feature2);
  map1.updateRanges();
  map1.setIdentifier("lsid");
  map1.getDataProcessing().resize(1);

  FeatureMap e2;

  //############################################################################################################
  // default test case ####################################
  //dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases.featureXML"), e2);
  //dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases_userparams.featureXML"), e2);
  //dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_output.featureXML"), e2);

  // metabolomics 1.
  //dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("metabolomics_short.featureXML"), e2);
  //dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("metabolomics.featureXML"), e2);
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("proteomics.featureXML"), e2);

  
  //############################################################################################################

  // save as sql database
  FeatureSQLFile fsf;

  // default test case ####################################
  //fsf.write("MetaboIdent_1_allcases_outputSQL", e2);
  //fsf.write("MetaboIdent_1_allcases_userparams_outputSQL", e2);
  //fsf.write("MetaboIdent_1_output_outputSQL", e2);
  
  // metabolomics 1.
  //fsf.write("metabolomics_outputSQL", e2);
  fsf.write("proteomics_outputSQL", e2);
  //############################################################################################################

}
END_SECTION



START_SECTION((FeatureMap read(std::string& in_featureSQL)))
{
  FeatureMap output;
  FeatureSQLFile fsf;

  //############################################################################################################
  // default test case ####################################
  //output = fsf.read(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases_outputSQL"));
  //output = fsf.read(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases_userparams_outputSQL"));
  //output = fsf.read(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_output_outputSQL"));
  
  // metabolomics 1.
  //output = fsf.read(OPENMS_GET_TEST_DATA_PATH("metabolomics_short_outputSQL"));  
  //output = fsf.read(OPENMS_GET_TEST_DATA_PATH("metabolomics_outputSQL"));  
  output = fsf.read(OPENMS_GET_TEST_DATA_PATH("proteomics_outputSQL"));  
  //############################################################################################################

  FeatureXMLFile f;

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, output);
  WHITELIST("?xml-stylesheet")

  //############################################################################################################
  // test against input featureXML file
  // default test case ####################################
  //TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases.featureXML"), tmp_filename)  
  //TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_allcases_userparams.featureXML"), tmp_filename)  
  //TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("MetaboIdent_1_output.featureXML"), tmp_filename)  
  // metabolomics 1.
  //TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("metabolomics_short.featureXML"), tmp_filename)  
  //TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("metabolomics.featureXML"), tmp_filename)  
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("proteomics.featureXML"), tmp_filename)  
  //############################################################################################################


}
END_SECTION



/////////////////////////////////////////////////////////////
END_TEST










