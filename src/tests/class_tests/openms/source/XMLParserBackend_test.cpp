// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <future>
#include <stdexcept>

using namespace OpenMS;
namespace
{
  struct Reader : Internal::XMLFile
  {
    using XMLFile::parseBuffer_;
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
