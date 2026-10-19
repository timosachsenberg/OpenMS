// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// Standalone baseline used by check_schemas.py; 0 = accepted, 1 = rejected.
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <memory>

struct Errors : xercesc::DefaultHandler
{
  bool failed = false;
  void error(const xercesc::SAXParseException&) override { failed = true; }
  void fatalError(const xercesc::SAXParseException&) override { failed = true; }
  void warning(const xercesc::SAXParseException&) override { failed = true; }
};
int main(int argc, char** argv)
{
  if (argc < 2 || argc > 3) { return 2; }
  xercesc::XMLPlatformUtils::Initialize();
  int result = 1;
  {
    Errors errors;
    std::unique_ptr<xercesc::SAX2XMLReader> parser(xercesc::XMLReaderFactory::createXMLReader());
    parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, true);
    parser->setFeature(xercesc::XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(xercesc::XMLUni::fgXercesSchema, true);
    parser->setFeature(xercesc::XMLUni::fgXercesSchemaFullChecking, true);
    parser->setErrorHandler(&errors);
    try
    {
      auto* grammar = parser->loadGrammar(argv[1], xercesc::Grammar::SchemaGrammarType, true);
      if (grammar && !errors.failed)
      {
        parser->setFeature(xercesc::XMLUni::fgXercesUseCachedGrammarInParse, true);
        // Compare the explicitly supplied bundled grammar. Older fixtures may
        // reference remote schema versions that are outside the bundled corpus.
        parser->setFeature(xercesc::XMLUni::fgXercesLoadSchema, false);
        if (argc == 3) { parser->parse(argv[2]); }
        result = errors.failed ? 1 : 0;
      }
    }
    catch (...) { result = 1; }
  }
  xercesc::XMLPlatformUtils::Terminate();
  return result;
}
