// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <exception>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage: xml_schema_validation_probe XML_FILE XSD_FILE\n";
    return 2;
  }
  try
  {
    const bool valid = OpenMS::XMLValidator().isValid(argv[1], argv[2]);
    std::cout << "{\"valid\":" << (valid ? "true" : "false")
              << ",\"file\":" << std::quoted(argv[1])
              << ",\"schema\":" << std::quoted(argv[2]) << "}\n";
    return valid ? 0 : 1;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
