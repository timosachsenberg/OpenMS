// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <fstream>
///////////////////////////

#define MULTI_LINE_STRING(...) #__VA_ARGS__ 

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IndexedMzMLDecoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexedMzMLDecoder* ptr = nullptr;
IndexedMzMLDecoder* nullPointer = nullptr;
START_SECTION((IndexedMzMLDecoder()))
  ptr = new IndexedMzMLDecoder;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IndexedMzMLDecoder()))
  delete ptr;
END_SECTION

START_SECTION((int parseOffsets(std::string filename, std::streampos indexoffset, OffsetVector & spectra_offsets, OffsetVector& chromatograms_offsets)))
  // see also IndexedMzMLFile_test.cpp
  std::streampos res = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(res, -1)

  IndexedMzMLDecoder::OffsetVector spectra_offsets;
  IndexedMzMLDecoder::OffsetVector chromatograms_offsets;

  int off = IndexedMzMLDecoder().parseOffsets(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), res, spectra_offsets, chromatograms_offsets);

  TEST_EQUAL(off, 0)
  TEST_EQUAL(spectra_offsets.size(), 2)
  TEST_EQUAL(chromatograms_offsets.size(), 1)
END_SECTION

START_SECTION([EXTRA] parseOffsets with no whitespace inside <index>)
  // the first <offset> is the first child of <index> (no text node before it)
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  {
    std::ofstream os(tmp_filename.c_str());
    os << "<indexList count=\"2\">"
          "<index name=\"spectrum\"><offset idRef=\"s1\">100</offset><offset idRef=\"s2\">200</offset></index>"
          "<index name=\"chromatogram\"><offset idRef=\"TIC\">300</offset></index>"
          "</indexList>\n<indexListOffset>0</indexListOffset>\n<fileChecksum>0</fileChecksum>\n</indexedmzML>\n";
  }
  IndexedMzMLDecoder::OffsetVector spectra_offsets;
  IndexedMzMLDecoder::OffsetVector chromatograms_offsets;
  int off = IndexedMzMLDecoder().parseOffsets(tmp_filename, 0, spectra_offsets, chromatograms_offsets);
  TEST_EQUAL(off, 0)
  TEST_EQUAL(spectra_offsets.size(), 2)
  ABORT_IF(spectra_offsets.size() != 2)
  TEST_EQUAL(spectra_offsets[0].first, "s1")
  TEST_EQUAL(spectra_offsets[0].second, 100)
  TEST_EQUAL(spectra_offsets[1].first, "s2")
  TEST_EQUAL(spectra_offsets[1].second, 200)
  TEST_EQUAL(chromatograms_offsets.size(), 1)
  ABORT_IF(chromatograms_offsets.size() != 1)
  TEST_EQUAL(chromatograms_offsets[0].first, "TIC")
  TEST_EQUAL(chromatograms_offsets[0].second, 300)
END_SECTION

    
START_SECTION((std::streampos findIndexListOffset(std::string filename, int buffersize = 1023)))
  // see also IndexedMzMLFile_test.cpp
  //
  std::streampos res = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(res, -1)

  // A std::streamoff value of -1 is also used to represent error conditions by some of the I/O library functions. 
  std::streampos nonindex = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(nonindex, -1)

END_SECTION

    

START_SECTION(([EXTRA] Tiny files clamp footer reads and parse directly from one buffer))
  std::string filename;
  NEW_TMP_FILE(filename)
  const std::string footer = "<indexList count='1'><index name='spectrum'><offset idRef='tiny'>37</offset>"
    "</index></indexList><indexListOffset>0</indexListOffset></indexedmzML>";
  { std::ofstream output(filename, std::ios::binary); output << footer; }
  IndexedMzMLDecoder decoder;
  TEST_EQUAL(decoder.findIndexListOffset(filename), 0)
  TEST_EQUAL(decoder.findIndexListOffset(filename, 1), -1)
  TEST_EQUAL(decoder.findIndexListOffset(filename, 0), -1)
  TEST_EQUAL(decoder.findIndexListOffset(filename, -17), -1)
  IndexedMzMLDecoder::OffsetVector spectra, chromatograms;
  TEST_EQUAL(decoder.parseOffsets(filename, 0, spectra, chromatograms), 0)
  TEST_EQUAL(spectra.size(), 1)
  if (!spectra.empty())
  {
    TEST_STRING_EQUAL(spectra[0].first, "tiny")
    TEST_EQUAL(spectra[0].second, 37)
  }
  TEST_TRUE(chromatograms.empty())
  TEST_EQUAL(decoder.parseOffsets(filename, -1, spectra, chromatograms), -1)
  TEST_EQUAL(decoder.parseOffsets(filename, std::streampos(footer.size()), spectra, chromatograms), -1)
  TEST_EQUAL(decoder.parseOffsets(filename, std::streampos(footer.size() + 1), spectra, chromatograms), -1)
END_SECTION

START_SECTION(([EXTRA] Empty missing and truncated footer tokens return no offset))
  std::string filename;
  NEW_TMP_FILE(filename)
  IndexedMzMLDecoder decoder;
  for (const std::string& text : {std::string(), std::string("<root/>"), std::string("<indexListOff"),
       std::string("<indexListOffset></indexListOffset>"), std::string("<indexListOffset>-1</indexListOffset>")})
  {
    { std::ofstream output(filename, std::ios::binary); output << text; }
    TEST_EQUAL(decoder.findIndexListOffset(filename), -1)
  }
  // A file truncated to empty must be rejected before entering a backend DOM parser.
  { std::ofstream output(filename, std::ios::binary); }
  IndexedMzMLDecoder::OffsetVector spectra, chromatograms;
  TEST_EQUAL(decoder.parseOffsets(filename, 0, spectra, chromatograms), -1)
END_SECTION

START_SECTION(([EXTRA] Namespaced footer pattern and overflowing offsets preserve conversion semantics))
  std::string filename;
  NEW_TMP_FILE(filename)
  IndexedMzMLDecoder decoder;
  { std::ofstream output(filename, std::ios::binary); output << "<ns:indexListOffset > 123</ns:indexListOffset>"; }
  TEST_EQUAL(decoder.findIndexListOffset(filename), 123)
  { std::ofstream output(filename, std::ios::binary); output << "<indexListOffset>" << std::string(100, '9') << "</indexListOffset>"; }
  TEST_EXCEPTION(Exception::ConversionError, decoder.findIndexListOffset(filename))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

