// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

namespace
{
  struct Reader : OpenMS::Internal::XMLFile
  {
    using XMLFile::parseBuffer_;
  };
  struct Handler : OpenMS::Internal::XMLHandler
  {
    Handler() : XMLHandler("memory", "1.0") {}
    std::u16string text;
    void onStartElement(const char16_t*, const OpenMS::Internal::XMLAttributes&) override {}
    void onEndElement(const char16_t*) override {}
    void onCharacters(const char16_t* value, OpenMS::Size length) override { text.append(value,length); }
  };
}

int main(int argc, char** argv)
{
  if (argc != 2) return 2;
  std::ifstream input(argv[1],std::ios::binary);
  if (!input) return 2;
  const std::string xml((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
  Reader reader;
  Handler handler;
  try
  {
    reader.parseBuffer_(xml,&handler);
    std::cout << "{\"accepted\":true,\"utf16_code_units\":[";
    bool first=true;
    for (char16_t value:handler.text)
    {
      if (!first) std::cout << ',';
      first=false;
      std::cout << static_cast<unsigned int>(value);
    }
    std::cout << "]}\n";
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    std::cout << "{\"accepted\":false}\n";
  }
}
