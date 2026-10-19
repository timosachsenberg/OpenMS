// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>

#include <libxml/parser.h>
#include <memory>
#include <limits>

namespace OpenMS
{

  namespace IndexedMzMLUtils
  {
    std::streampos stringToStreampos(const std::string& s)
    {
      // Try to cast the string to a type that can hold the integer value
      // stored in the std::streampos type (which can vary from 32 to 128 bits
      // depending on the system).
      //
      // long long has a minimal size of 64 bits and can address a range up to
      // 16 Exbibit (or 2 Exbibyte), we can hopefully expect our files to be
      // smaller than an Exabyte and should be safe.
      std::streampos res;
      try
      {

        res = boost::lexical_cast< unsigned long long >(s);
        // TESTING CAST: res = (int) res; // use this only when emulating 32 bit systems

      }
      catch (boost::bad_lexical_cast&)
      {
        std::cerr << "Trying to convert corrupted / unreadable value to std::streampos : " << s << std::endl;
        std::cerr << "This can also happen if the value exceeds 63 bits, please check your input." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to a 64 bit integer.");
      }

      // Check if the value can fit into std::streampos
      if ( std::abs( boost::lexical_cast< long double >(s) - res) > 0.1)
      {
        std::cerr << "Your system may not support addressing a file of this size,"
          << " only addresses that fit into a " << sizeof(std::streamsize)*8 <<
          " bit integer are supported on your system." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to an integer on your system.");
      }

      return res;
    }
  }

  int IndexedMzMLDecoder::parseOffsets(const std::string& filename, std::streampos indexoffset, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    //-------------------------------------------------------------
    // Open file, jump to end and read last indexoffset bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());
    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // get length of file:
    f.seekg(0, f.end);
    std::streampos length = f.tellg();

    if (indexoffset < 0 || indexoffset > length)
    {
      std::cerr << "IndexedMzMLDecoder::parseOffsets Error: Offset was " <<
        indexoffset << " (not between 0 and " << length << ")." << std::endl;
      return -1;
    }

    //-------------------------------------------------------------
    // Read full end of file to parse offsets for spectra and chroms
    //-------------------------------------------------------------
    // read data as a block into a buffer

    // allocate enough memory in buffer (+1 for string termination)
    std::streampos readl = length - indexoffset;
    char* buffer = new(std::nothrow) char[readl + std::streampos(1)];

    // catch case where not enough memory is available
    if (buffer == nullptr)
    {
      // Warning: Index takes up more than 10 % of the whole file, please check your input file.\n";
      std::cerr << "IndexedMzMLDecoder::parseOffsets Could not allocate enough memory to read in index of indexedMzML" << std::endl;
      std::cerr << "IndexedMzMLDecoder::parseOffsets calculated index offset " << indexoffset << " and file length " << length <<
        ", consequently tried to read into memory " << readl << " bytes." << std::endl;
      return -1;
    }

    // read into memory
    f.seekg(-readl, f.end);
    f.read(buffer, readl);
    buffer[readl] = '\0';

    //-------------------------------------------------------------
    // Add a sane start element and then give it to a DOM parser
    //-------------------------------------------------------------
    // http://stackoverflow.com/questions/4691039/making-xerces-parse-a-string-insted-of-a-file
    std::string tmp_fixed_xml = "<indexedmzML>" +  std::string(buffer) + "\n";
    int res = domParseIndexedEnd_(tmp_fixed_xml, spectra_offsets, chromatograms_offsets);

    delete[] buffer;

    return res;
  }

  std::streampos IndexedMzMLDecoder::findIndexListOffset(const std::string& filename, int buffersize)
  {
    // return value
    std::streampos indexoffset = -1;

    //-------------------------------------------------------------
    // Open file, jump to end and read last n bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());

    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // Read the last few bytes and hope our offset is there to be found
    std::unique_ptr<char[]> buffer(new char[buffersize + 1]);
    f.seekg(-buffersize, f.end);
    f.read(buffer.get(), buffersize);
    buffer.get()[buffersize] = '\0';

#ifdef DEBUG_READER
    std::cout << " reading file " << filename  << " with size " << buffersize << '\n';
    std::cout << buffer << '\n';
#endif

    //-------------------------------------------------------------
    // Since we could be anywhere in the XML structure, use regex to find
    // indexListOffset and read its content.
    //-------------------------------------------------------------
    boost::regex listoffset_rx(R"(<[^>/]*indexListOffset\s*>\s*(\d*))");
    boost::cmatch matches;
    boost::regex_search(buffer.get(), matches, listoffset_rx);
    std::string thismatch(matches[1].first, matches[1].second);
    if (!thismatch.empty())
    {
      try
      {
        indexoffset = OpenMS::IndexedMzMLUtils::stringToStreampos(thismatch);
      }
      catch (Exception::ConversionError& /*e*/)
      {
        std::cerr << "Corrupted / unreadable value in <indexListOffset> : " << thismatch << std::endl;
        // free resources and re-throw
        throw;  // re-throw conversion error
      }
    }
    else
    {
      std::cerr << "IndexedMzMLDecoder::findIndexListOffset Error: Could not find element indexListOffset in the last "
                << buffersize << " bytes. Maybe this is not a indexedMzML."
                << buffer.get() << std::endl;
    }

    return indexoffset;
  }

  int IndexedMzMLDecoder::domParseIndexedEnd_(const std::string& in, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {

    if (in.size() > static_cast<Size>(std::numeric_limits<int>::max())) { return -1; }
    std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
      xmlReadMemory(in.data(), static_cast<int>(in.size()), "index.xml", nullptr,
                    XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT | XML_PARSE_NOERROR | XML_PARSE_NOWARNING), xmlFreeDoc);
    xmlNode* root = document ? xmlDocGetRootElement(document.get()) : nullptr;
    if (!root) { return -1; }
    xmlNode* list = nullptr;
    for (xmlNode* node = root->children; node; node = node->next)
    {
      if (node->type == XML_ELEMENT_NODE && xmlStrEqual(node->name, BAD_CAST "indexList"))
      {
        if (list) { return -1; }
        list = node;
      }
    }
    if (!list) { return -1; }
    auto property = [](xmlNode* node, const char* key)
    {
      xmlChar* value = xmlGetProp(node, BAD_CAST key);
      std::string result = value ? reinterpret_cast<const char*>(value) : "";
      xmlFree(value);
      return result;
    };
    for (xmlNode* index = list->children; index; index = index->next)
    {
      if (index->type != XML_ELEMENT_NODE) { continue; }
      OffsetVector result;
      for (xmlNode* offset = index->children; offset; offset = offset->next)
      {
        if (offset->type != XML_ELEMENT_NODE) { continue; }
        xmlChar* content = xmlNodeGetContent(offset);
        std::string text = content ? reinterpret_cast<const char*>(content) : "";
        xmlFree(content);
        result.emplace_back(property(offset, "idRef"), IndexedMzMLUtils::stringToStreampos(text));
      }
      const auto name = property(index, "name");
      if (name == "spectrum") { spectra_offsets = std::move(result); }
      else if (name == "chromatogram") { chromatograms_offsets = std::move(result); }
      else { return -1; }
    }
    return 0;
  }

}
