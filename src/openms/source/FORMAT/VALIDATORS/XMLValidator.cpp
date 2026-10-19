// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Init.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>

#include <libxml/xmlschemas.h>
#include <memory>

using namespace std;

namespace OpenMS
{
  class XMLValidatorErrorHandler_
  {
  public:
    static void error(void* ctx, const xmlError* error)
    {
      // Context-local callback: independent validators do not share error state.
      auto& validator = *static_cast<XMLValidator*>(ctx);
      validator.logError_(error->level == XML_ERR_WARNING ? "warning" : "error",
                          error->message ? error->message : "libxml2 validation error",
                          error->line, error->int2);
    }
  };

  XMLValidator::XMLValidator() : valid_(true), os_(nullptr) {}

  bool XMLValidator::isValid(const std::string& filename, const std::string& schema, std::ostream& os)
  {
    filename_ = filename;
    os_ = &os;
    valid_ = true;
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    unique_ptr<xmlSchemaParserCtxt, decltype(&xmlSchemaFreeParserCtxt)> parser(
      xmlSchemaNewParserCtxt(schema.c_str()), xmlSchemaFreeParserCtxt);
    if (!parser) { throw std::bad_alloc(); }
    xmlSchemaSetParserStructuredErrors(parser.get(), XMLValidatorErrorHandler_::error, this);
    unique_ptr<xmlSchema, decltype(&xmlSchemaFree)> grammar(xmlSchemaParse(parser.get()), xmlSchemaFree);
    if (!grammar) { return false; }
    unique_ptr<xmlSchemaValidCtxt, decltype(&xmlSchemaFreeValidCtxt)> validator(
      xmlSchemaNewValidCtxt(grammar.get()), xmlSchemaFreeValidCtxt);
    if (!validator) { throw std::bad_alloc(); }
    xmlSchemaSetValidStructuredErrors(validator.get(), XMLValidatorErrorHandler_::error, this);
    // Streaming validation preserves bounded memory for large mzML inputs.
    const int result = xmlSchemaValidateFile(validator.get(), filename.c_str(), XML_PARSE_NONET | XML_PARSE_NO_XXE);
    return valid_ && result == 0;
  }

  void XMLValidator::logError_(const std::string& kind, const std::string& message, unsigned long line, unsigned long column)
  {
    (*os_) << "Validation " << kind << " in file '" << filename_ << "' line " << line
           << " column " << column << ": " << message << '\n';
    valid_ = false;
  }

} // namespace OpenMS
