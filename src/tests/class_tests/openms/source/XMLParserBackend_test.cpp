// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/test_config.h>
#include <fstream>
#include <future>
#include <stdexcept>

using namespace OpenMS;
namespace
{
  struct Reader : Internal::XMLFile
  {
    using XMLFile::parseBuffer_;
    using XMLFile::enforceEncoding_;
  };
  struct Handler : Internal::XMLHandler
  {
    Handler() : XMLHandler("memory", "1.0") {}
    std::u16string text, attribute;
    Size starts = 0;
    bool stop = false, fail = false;
    void onStartElement(const char16_t*, const Internal::XMLAttributes& attrs) override
    {
      ++starts;
      if (stop) { throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
      if (fail) { throw std::runtime_error("callback exception"); }
      if (attrs.value("value"))
      {
        attribute = attrs.value("value");
        if (attrs.index(u"value") < 0 || !attrs.valueByIndex(attrs.index(u"value")))
        {
          throw std::runtime_error("attribute lookup mismatch");
        }
      }
    }
    void onCharacters(const char16_t* value, Size length) override { text.append(value, length); }
    void onEndElement(const char16_t*) override {}
  };
}
START_TEST(XMLParserBackend, "$Id$")
START_SECTION((UTF8 Unicode, entities, CDATA, attributes, and chunk boundaries))
  Reader reader;
  Handler handler;
  const std::string repeated(70000, 'x');
  reader.parseBuffer_("<root value='Grüße &amp; &#x1F9EC;'>&lt;<![CDATA[&]]>" + repeated + "🧬</root>", &handler);
  TEST_EQUAL(handler.attribute == u"Grüße & 🧬", true)
  TEST_EQUAL(handler.text == std::u16string(u"<&") + std::u16string(70000, u'x') + u"🧬", true)
  TEST_EQUAL(handler.starts, 1)
END_SECTION
START_SECTION(Malformed XML and exception propagation)
  Reader reader;
  Handler malformed;
  TEST_EXCEPTION(Exception::ParseError, reader.parseBuffer_("<root><broken></root>", &malformed))
  Handler failure;
  failure.fail = true;
  TEST_EXCEPTION(std::runtime_error, reader.parseBuffer_("<root/>", &failure))
  Handler stopped;
  stopped.stop = true;
  reader.parseBuffer_("<root><ignored/></root>", &stopped);
  TEST_EQUAL(stopped.starts, 1)
END_SECTION
START_SECTION(External entity loading is disabled)
  std::string filename;
  NEW_TMP_FILE(filename)
  { std::ofstream out(filename); out << "external secret"; }
  Reader reader;
  Handler handler;
  const std::string xml = "<!DOCTYPE root [<!ENTITY external SYSTEM 'file://" + filename + "'>]><root>&external;</root>";
  reader.parseBuffer_(xml, &handler);
  TEST_TRUE(handler.text.empty())
END_SECTION
START_SECTION(Validator stream exceptions are rethrown after C callbacks return)
  struct FailingBuffer : std::streambuf
  {
    int_type overflow(int_type) override { return traits_type::eof(); }
  } buffer;
  std::ostream output(&buffer);
  output.exceptions(std::ios::badbit | std::ios::failbit);
  XMLValidator validator;
  TEST_EXCEPTION(std::ios_base::failure,
    validator.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_missing_attribute.xml"),
                      OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd"), output))
END_SECTION
START_SECTION(Encoding overrides are per reader)
  Reader latin, utf8;
  latin.enforceEncoding_("ISO-8859-1");
  utf8.enforceEncoding_("UTF-8");
  Handler latin_handler, utf8_handler;
  latin.parseBuffer_("<root>Gr\xfc\xdf" "e</root>", &latin_handler);
  utf8.parseBuffer_("<root>Grüße</root>", &utf8_handler);
  TEST_EQUAL(latin_handler.text == u"Grüße", true)
  TEST_EQUAL(utf8_handler.text == u"Grüße", true)
END_SECTION
START_SECTION(Indexed decoder supports a binary text node larger than 10 MB)
  const std::string large_payload = std::string(13 * 1024 * 1024, ' ') + "AAAAAAAAAAA=";
  const auto array = [](const std::string& kind, const std::string& payload)
  {
    return "<binaryDataArray><cvParam accession='MS:1000523' name='64-bit float'/>"
           "<cvParam accession='" + kind + "'/><binary>" + payload + "</binary></binaryDataArray>";
  };
  const std::string xml = "<spectrum id='large' defaultArrayLength='1'><binaryDataArrayList count='2'>" +
    array("MS:1000514", large_payload) + array("MS:1000515", "AAAAAAAAAAA=") + "</binaryDataArrayList></spectrum>";
  MzMLSpectrumDecoder decoder;
  MSSpectrum spectrum;
  decoder.domParseSpectrum(xml, spectrum);
  TEST_EQUAL(spectrum.size(), 1)
  TEST_EQUAL(spectrum[0].getMZ(), 0.0)
  TEST_EQUAL(spectrum[0].getIntensity(), 0.0)
END_SECTION
START_SECTION(Indexed decoder retains the entity amplification guard under HUGE)
  std::string dtd = "<!DOCTYPE spectrum [<!ENTITY e0 '0123456789'>";
  for (int i = 1; i <= 7; ++i)
  {
    dtd += "<!ENTITY e" + std::to_string(i) + " '";
    for (int j = 0; j < 10; ++j) { dtd += "&e" + std::to_string(i - 1) + ";"; }
    dtd += "'>";
  }
  const std::string xml = dtd + "]><spectrum defaultArrayLength='0'>&e7;</spectrum>";
  MzMLSpectrumDecoder decoder;
  MSSpectrum spectrum;
  TEST_EXCEPTION(Exception::ParseError, decoder.domParseSpectrum(xml, spectrum))
END_SECTION
START_SECTION(Independent parser contexts on eight threads)
  std::vector<std::future<bool>> workers;
  for (int i = 0; i < 8; ++i)
  {
    workers.emplace_back(std::async(std::launch::async, []
    {
      for (int j = 0; j < 100; ++j)
      {
        Reader reader;
        Handler handler;
        reader.parseBuffer_("<root value='parallel'>text</root>", &handler);
        if (handler.text != u"text" || handler.attribute != u"parallel") { return false; }
      }
      return true;
    }));
  }
  for (auto& worker : workers) { TEST_TRUE(worker.get()) }
END_SECTION
END_TEST
