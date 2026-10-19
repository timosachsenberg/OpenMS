// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XMLFile.h>

#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/FileNameUtils.h>

#include <OpenMS/FORMAT/HANDLERS/SAX2HandlerAdapter.h>

#include <libxml/parser.h>
#include <libxml/parserInternals.h>
#include <xercesc/util/BinInputStream.hpp>
#include <exception>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>


#include <zlib.h>
#include <bzlib.h>

#include <fstream>
#include <iomanip> // setprecision etc.
#include <streambuf>

#include <memory>
#include <xercesc/util/XMLString.hpp>

using namespace std;

namespace OpenMS::Internal
{

    /// This class ensures that the reset() method of the XMLHandler is called when it goes out of scope.
    /// useful when used in exception handling
    class XMLCleaner_
    {
public:
      explicit XMLCleaner_(XMLHandler * handler) :
        p_(handler)
      {

      }

      ~XMLCleaner_()
      {
        p_->reset();
      }

private:
      XMLHandler * p_;
    };

    XMLFile::XMLFile()
    = default;

    XMLFile::XMLFile(const std::string & schema_location, const std::string & version) :
      schema_location_(schema_location),
      schema_version_(version)
    {
    }

    XMLFile::~XMLFile()
    = default;

    void XMLFile::enforceEncoding_(const std::string& encoding)
    {
      enforced_encoding_ = encoding;
    }

    // POC: keep the existing decompression streams, replace their XML parser.
    // C callbacks capture every exception: no C++ exception crosses libxml2 frames.
    struct LibxmlState
    {
      XMLHandler* handler;
      std::exception_ptr exception;
      static LibxmlState& state(void* ctx)
      {
        return *static_cast<LibxmlState*>(static_cast<xmlParserCtxt*>(ctx)->_private);
      }
      static std::u16string utf16(const xmlChar* text, int length = -1)
      {
        std::u16string result;
        const unsigned char* cursor = text;
        const unsigned char* end = text + (length < 0 ? xmlStrlen(text) : length);
        while (cursor < end)
        {
          unsigned int cp = *cursor++;
          if (cp >= 0xc0)
          {
            const int remaining = cp < 0xe0 ? 1 : cp < 0xf0 ? 2 : 3;
            cp &= (1u << (6 - remaining)) - 1;
            for (int i = 0; i < remaining && cursor < end; ++i) { cp = (cp << 6) | (*cursor++ & 0x3f); }
          }
          if (cp > 0xffff)
          {
            cp -= 0x10000;
            result.push_back(static_cast<char16_t>(0xd800 + (cp >> 10)));
            result.push_back(static_cast<char16_t>(0xdc00 + (cp & 0x3ff)));
          }
          else { result.push_back(static_cast<char16_t>(cp)); }
        }
        return result;
      }
      template <typename Callback> static void callback(void* ctx, Callback&& fn)
      {
        // Stopping the parser can re-enter an error callback. Keep the original
        // handler exception rather than replacing it with a parser-stop error.
        if (state(ctx).exception) { return; }
        try { fn(state(ctx)); }
        catch (...)
        {
          state(ctx).exception = std::current_exception();
          xmlStopParser(static_cast<xmlParserCtxt*>(ctx));
        }
      }
      static void start(void* ctx, const xmlChar* name, const xmlChar** attrs)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          std::vector<std::pair<std::u16string, std::u16string>> values;
          if (attrs)
          {
            for (Size i = 0; attrs[i]; i += 2) { values.emplace_back(utf16(attrs[i]), utf16(attrs[i + 1])); }
          }
          s.handler->onStartElement(utf16(name).c_str(), XMLAttributes(values));
        });
      }
      static void end(void* ctx, const xmlChar* name)
      {
        callback(ctx, [&](LibxmlState& s) { s.handler->onEndElement(utf16(name).c_str()); });
      }
      static void characters(void* ctx, const xmlChar* text, int length)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          if (s.handler->tryOnCharactersUTF8(reinterpret_cast<const char*>(text), static_cast<Size>(length))) { return; }
          const auto value = utf16(text, length);
          s.handler->onCharacters(value.c_str(), value.size());
        });
      }
      static void error(void* ctx, const xmlError* error)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          const std::string message = error->message ? error->message : "libxml2 parse error";
          if (error->level == XML_ERR_WARNING) { s.handler->warning(XMLHandler::LOAD, message, error->line, error->int2); }
          else { s.handler->fatalError(XMLHandler::LOAD, message, error->line, error->int2); }
        });
      }
    };
    void parseLibxmlPOC_(xercesc::InputSource* const source, XMLHandler* handler)
    {
      LibxmlState state{handler, {}};
      xmlSAXHandler sax{};
      xmlSAXVersion(&sax, 1);
      sax.startElement = LibxmlState::start;
      sax.endElement = LibxmlState::end;
      sax.characters = LibxmlState::characters;
      sax.cdataBlock = LibxmlState::characters;
      auto free_parser = [](xmlParserCtxt* ctx)
      {
        if (ctx->myDoc) { xmlFreeDoc(ctx->myDoc); }
        xmlFreeParserCtxt(ctx);
      };
      unique_ptr<xmlParserCtxt, decltype(free_parser)> parser(xmlCreatePushParserCtxt(&sax, nullptr, nullptr, 0, nullptr), free_parser);
      if (!parser) { throw std::bad_alloc(); }
      parser->_private = &state;
      parser->userData = parser.get();
      // SAX1 ignores sax.serror; the context-local API works in both SAX modes.
      xmlCtxtSetErrorHandler(parser.get(), LibxmlState::error, parser.get());
      xmlCtxtUseOptions(parser.get(), XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT);
      if (source->getEncoding())
      {
        const auto encoding = StringManager::convert(reinterpret_cast<const char16_t*>(source->getEncoding()));
        if (xmlSwitchEncodingName(parser.get(), encoding.c_str()) != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Unsupported XML encoding: " + encoding);
        }
      }
      unique_ptr<xercesc::BinInputStream> stream(source->makeStream());
      if (!stream) { throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Cannot open XML input stream"); }
      char buffer[65536];
      int status = 0;
      while (!state.exception)
      {
        const auto count = stream->readBytes(reinterpret_cast<XMLByte*>(buffer), sizeof(buffer));
        status = xmlParseChunk(parser.get(), buffer, static_cast<int>(count), count == 0);
        if (status || count == 0) { break; }
      }
      if (state.exception)
      {
        try { std::rethrow_exception(state.exception); }
        catch (const XMLHandler::EndParsingSoftly&) { return; }
      }
      if (status)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "libxml2 rejected XML input");
      }
    }

    // -----------------------------------------------------------------------
    // Benchmark prototype: runtime-selectable XML adapter.
    //
    // OPENMS_XML_ADAPTER selects how XMLFile feeds an XMLHandler:
    //   libxml2-poc  (default) the unchanged proof-of-concept adapter above
    //   xerces       Xerces SAX2 + SAX2HandlerAdapter, as on develop 77cccda
    //   libxml2-sax2 reworked adapter: libxml2 SAX2 (startElementNs, attribute
    //                values as pointers into the parser buffer), cached UTF-16
    //                element/attribute names, reused UTF-16 value buffers
    //   libxml2-utf8 libxml2-sax2 plus a UTF-8 path to the handler: attribute
    //                helpers read UTF-8 directly, UTF-16 values are converted on
    //                demand only, and MzMLHandler takes the UTF-8 element name
    //   libxml2-fast libxml2-utf8 plus (a) Latin-1 input decoded by the adapter:
    //                OpenMS writes mzML with encoding="ISO-8859-1"; libxml2
    //                would run its byte-wise latin1ToUTF8 over the whole file.
    //                The adapter passes pure-ASCII chunks through unchanged
    //                (Latin-1 == UTF-8 there) and converts other chunks itself,
    //                so decoding stays exact; (b) MzMLHandler's UTF-8 character
    //                hook checks "data discarded" before its ASCII scan.
    // The handler API (UTF-16 callbacks) is unchanged in all modes.
    // -----------------------------------------------------------------------
    enum class XMLAdapterMode { Xerces, LibxmlPOC, LibxmlSAX2, LibxmlUTF8, LibxmlFast };

    XMLAdapterMode xmlAdapterMode_()
    {
      static const XMLAdapterMode mode = []
      {
        const char* value = std::getenv("OPENMS_XML_ADAPTER");
        const std::string name = value ? value : "";
        XMLAdapterMode result;
        if (name.empty() || name == "libxml2-poc") { result = XMLAdapterMode::LibxmlPOC; }
        else if (name == "xerces") { result = XMLAdapterMode::Xerces; }
        else if (name == "libxml2-sax2") { result = XMLAdapterMode::LibxmlSAX2; }
        else if (name == "libxml2-utf8") { result = XMLAdapterMode::LibxmlUTF8; }
        else if (name == "libxml2-fast") { result = XMLAdapterMode::LibxmlFast; }
        else
        {
          std::fprintf(stderr, "Unknown OPENMS_XML_ADAPTER '%s'\n", name.c_str());
          std::abort();
        }
        if (std::getenv("OPENMS_XML_ADAPTER_VERBOSE"))
        {
          std::fprintf(stderr, "OpenMS XML adapter: %s\n", name.empty() ? "libxml2-poc" : name.c_str());
        }
        return result;
      }();
      return mode;
    }

    void parseXerces_(xercesc::InputSource* const source, XMLHandler* handler)
    {
      // Verbatim from develop 77cccda7caf40.
      unique_ptr<xercesc::SAX2XMLReader> parser(xercesc::XMLReaderFactory::createXMLReader());

      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

      // XMLHandler exposes only Xerces-free callbacks; bridge it to the SAX parser.
      SAX2HandlerAdapter adapter(*handler);
      parser->setContentHandler(&adapter);
      parser->setErrorHandler(&adapter);

      // try to parse file
      try
      {
        parser->parse(*source);
      }
      catch (const xercesc::XMLException& toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("XMLException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const xercesc::SAXException& toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("SAXException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const XMLHandler::EndParsingSoftly& /*toCatch*/)
      {
        // nothing to do here, as this exception is used to softly abort the
        // parsing for whatever reason.
      }
      catch (...)
      {
        // re-throw
        throw;
      }
    }

    struct LibxmlSAX2State
    {
      struct Name
      {
        std::string utf8;
        std::u16string utf16;
      };

      struct KeyHash
      {
        size_t operator()(const std::pair<const void*, const void*>& key) const noexcept
        {
          const auto a = reinterpret_cast<std::uintptr_t>(key.first);
          const auto b = reinterpret_cast<std::uintptr_t>(key.second);
          return std::hash<std::uintptr_t>()(a ^ (b * 0x9e3779b97f4a7c15ULL));
        }
      };

      /// RAII: the UTF-8 element-name hint is only valid during onStartElement.
      struct QNameHint
      {
        XMLHandler* handler;
        ~QNameHint() { if (handler) { handler->setNativeQNameUTF8(nullptr); } }
      };

      XMLHandler* handler = nullptr;
      std::exception_ptr exception;
      bool utf8 = false;

      // Element and attribute names come from libxml2's per-parser dictionary, so
      // their pointers are stable during a parse. Each cache hit is still verified
      // against the actual string, so correctness never depends on interning.
      std::unordered_map<std::pair<const void*, const void*>, Name, KeyHash> names;

      // Grow-only attribute buffers, reused for every element.
      std::vector<const char*> names8;
      std::vector<const char*> values8;
      std::vector<Size> values8_length;
      std::vector<const std::u16string*> names16;
      std::vector<std::u16string> values16;
      std::vector<unsigned char> ready;
      std::u16string chars16;

      static inline const char xmlns_marker = 0;

      static LibxmlSAX2State& state(void* ctx)
      {
        return *static_cast<LibxmlSAX2State*>(static_cast<xmlParserCtxt*>(ctx)->_private);
      }

      static bool matches(const std::string& qname, const char* prefix, const char* local)
      {
        const char* q = qname.c_str();
        if (prefix)
        {
          while (*prefix) { if (*q++ != *prefix++) { return false; } }
          if (*q++ != ':') { return false; }
        }
        return std::strcmp(q, local) == 0;
      }

      const Name& lookup(const void* key, const char* prefix, const char* local)
      {
        auto [it, inserted] = names.try_emplace(std::make_pair(key, static_cast<const void*>(local)));
        Name& entry = it->second;
        if (inserted || !matches(entry.utf8, prefix, local))
        {
          entry.utf8.clear();
          if (prefix) { entry.utf8 += prefix; entry.utf8 += ':'; }
          entry.utf8 += local;
          XMLAttributes::assignUTF16(entry.utf16, entry.utf8.data(), entry.utf8.size());
        }
        return entry;
      }

      const Name& name(const xmlChar* prefix, const xmlChar* local)
      {
        return lookup(prefix, reinterpret_cast<const char*>(prefix), reinterpret_cast<const char*>(local));
      }

      /// Namespace declarations re-exposed as "xmlns"/"xmlns:p" attributes, as Xerces does with namespaces off.
      const Name& xmlnsName(const xmlChar* prefix)
      {
        static const char xmlns[] = "xmlns";
        return prefix ? lookup(&xmlns_marker, xmlns, reinterpret_cast<const char*>(prefix))
                      : lookup(&xmlns_marker, nullptr, xmlns);
      }

      void reserve(Size count)
      {
        if (names8.size() < count)
        {
          names8.resize(count);
          values8.resize(count);
          values8_length.resize(count);
          names16.resize(count);
          values16.resize(count);
          ready.resize(count);
        }
      }

      template <typename Callback> static void callback(void* ctx, Callback&& fn)
      {
        if (state(ctx).exception) { return; }
        try { fn(state(ctx)); }
        catch (...)
        {
          state(ctx).exception = std::current_exception();
          xmlStopParser(static_cast<xmlParserCtxt*>(ctx));
        }
      }

      static void start(void* ctx, const xmlChar* localname, const xmlChar* prefix, const xmlChar* /*uri*/,
                        int nb_namespaces, const xmlChar** namespaces,
                        int nb_attributes, int /*nb_defaulted*/, const xmlChar** attributes)
      {
        callback(ctx, [&](LibxmlSAX2State& s)
        {
          const Size count = static_cast<Size>(nb_namespaces) + static_cast<Size>(nb_attributes);
          s.reserve(count);
          Size k = 0;
          for (int i = 0; i < nb_namespaces; ++i, ++k)
          {
            const Name& n = s.xmlnsName(namespaces[2 * i]);
            const xmlChar* uri = namespaces[2 * i + 1];
            s.names8[k] = n.utf8.c_str();
            s.names16[k] = &n.utf16;
            s.values8[k] = uri ? reinterpret_cast<const char*>(uri) : "";
            s.values8_length[k] = uri ? static_cast<Size>(xmlStrlen(uri)) : 0;
          }
          for (int i = 0; i < nb_attributes; ++i, ++k)
          {
            const xmlChar** a = attributes + 5 * i; // localname, prefix, URI, value begin, value end
            const Name& n = s.name(a[1], a[0]);
            s.names8[k] = n.utf8.c_str();
            s.names16[k] = &n.utf16;
            s.values8[k] = reinterpret_cast<const char*>(a[3]);
            s.values8_length[k] = static_cast<Size>(a[4] - a[3]);
          }
          if (s.utf8)
          {
            std::fill(s.ready.begin(), s.ready.begin() + count, static_cast<unsigned char>(0));
          }
          else
          {
            for (Size i = 0; i < count; ++i) { XMLAttributes::assignUTF16(s.values16[i], s.values8[i], s.values8_length[i]); }
          }
          XMLAttributes::NativeView view;
          view.count = count;
          view.names8 = s.names8.data();
          view.values8 = s.values8.data();
          view.values8_length = s.values8_length.data();
          view.names16 = s.names16.data();
          view.values16 = s.values16.data();
          view.values16_ready = s.ready.data();
          view.lazy_utf16 = s.utf8;
          view.expose_utf8 = s.utf8;

          const Name& element = s.name(prefix, localname);
          QNameHint hint{s.utf8 ? s.handler : nullptr};
          if (s.utf8) { s.handler->setNativeQNameUTF8(element.utf8.c_str()); }
          s.handler->onStartElement(element.utf16.c_str(), XMLAttributes(view));
        });
      }

      static void end(void* ctx, const xmlChar* localname, const xmlChar* prefix, const xmlChar* /*uri*/)
      {
        callback(ctx, [&](LibxmlSAX2State& s) { s.handler->onEndElement(s.name(prefix, localname).utf16.c_str()); });
      }

      static void characters(void* ctx, const xmlChar* text, int length)
      {
        callback(ctx, [&](LibxmlSAX2State& s)
        {
          const char* chars = reinterpret_cast<const char*>(text);
          if (s.handler->tryOnCharactersUTF8(chars, static_cast<Size>(length))) { return; }
          XMLAttributes::assignUTF16(s.chars16, chars, static_cast<Size>(length));
          s.handler->onCharacters(s.chars16.c_str(), s.chars16.size());
        });
      }

      static void error(void* ctx, const xmlError* error)
      {
        callback(ctx, [&](LibxmlSAX2State& s)
        {
          const std::string message = error->message ? error->message : "libxml2 parse error";
          if (error->level == XML_ERR_WARNING) { s.handler->warning(XMLHandler::LOAD, message, error->line, error->int2); }
          else { s.handler->fatalError(XMLHandler::LOAD, message, error->line, error->int2); }
        });
      }
    };

    /// Used by MzMLHandler (benchmark prototype): enable mode-specific handler fast paths.
    bool xmlAdapterFastPaths_()
    {
      return xmlAdapterMode_() == XMLAdapterMode::LibxmlFast;
    }

    namespace
    {
      bool isLatin1Name_(std::string name)
      {
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return name == "iso-8859-1" || name == "iso_8859-1" || name == "iso8859-1" || name == "latin1" ||
               name == "latin-1" || name == "iso-latin-1" || name == "l1" || name == "iso-ir-100" ||
               name == "cp819" || name == "ibm819" || name == "csisolatin1" || name == "iso_8859-1:1987";
      }

      /// Encoding named in an ASCII XML declaration at the very start of @p data; empty if absent/unreadable.
      std::string declaredEncoding_(const char* data, std::size_t size)
      {
        const std::string_view head(data, std::min<std::size_t>(size, 1024));
        if (head.substr(0, 5) != "<?xml") { return {}; } // BOM or non-ASCII-compatible encoding: leave to libxml2
        const auto decl_end = head.find("?>");
        if (decl_end == std::string_view::npos) { return {}; }
        const std::string_view decl = head.substr(0, decl_end);
        auto pos = decl.find("encoding");
        if (pos == std::string_view::npos) { return {}; }
        pos += 8;
        auto skip_space = [&] { while (pos < decl.size() && (decl[pos] == ' ' || decl[pos] == '\t' || decl[pos] == '\r' || decl[pos] == '\n')) { ++pos; } };
        skip_space();
        if (pos >= decl.size() || decl[pos] != '=') { return {}; }
        ++pos;
        skip_space();
        if (pos >= decl.size() || (decl[pos] != '"' && decl[pos] != '\'')) { return {}; }
        const char quote = decl[pos++];
        const auto value_end = decl.find(quote, pos);
        if (value_end == std::string_view::npos) { return {}; }
        return std::string(decl.substr(pos, value_end - pos));
      }

      bool isASCII_(const char* data, std::size_t size)
      {
        const auto* bytes = reinterpret_cast<const unsigned char*>(data);
        std::uint64_t accumulated = 0;
        std::size_t i = 0;
        for (; i + 8 <= size; i += 8)
        {
          std::uint64_t word;
          std::memcpy(&word, bytes + i, 8);
          accumulated |= word;
        }
        for (; i < size; ++i) { accumulated |= bytes[i]; }
        return (accumulated & 0x8080808080808080ULL) == 0;
      }

      void latin1ToUTF8_(const char* data, std::size_t size, std::string& out)
      {
        out.clear();
        out.reserve(2 * size);
        for (std::size_t i = 0; i < size; ++i)
        {
          const auto c = static_cast<unsigned char>(data[i]);
          if (c < 0x80) { out.push_back(static_cast<char>(c)); }
          else
          {
            out.push_back(static_cast<char>(0xc0 | (c >> 6)));
            out.push_back(static_cast<char>(0x80 | (c & 0x3f)));
          }
        }
      }
    }

    void parseLibxmlSAX2_(xercesc::InputSource* const source, XMLHandler* handler, bool utf8, bool fast = false)
    {
      LibxmlSAX2State state;
      state.handler = handler;
      state.utf8 = utf8;
      xmlSAXHandler sax{};
      xmlSAXVersion(&sax, 2);
      sax.startElement = nullptr;
      sax.endElement = nullptr;
      sax.startElementNs = LibxmlSAX2State::start;
      sax.endElementNs = LibxmlSAX2State::end;
      sax.characters = LibxmlSAX2State::characters;
      sax.cdataBlock = LibxmlSAX2State::characters;
      // Without a DOM, libxml2 never classifies text as ignorable; sharing the
      // callback lets it skip that per-chunk check.
      sax.ignorableWhitespace = LibxmlSAX2State::characters;
      auto free_parser = [](xmlParserCtxt* ctx)
      {
        if (ctx->myDoc) { xmlFreeDoc(ctx->myDoc); }
        xmlFreeParserCtxt(ctx);
      };
      unique_ptr<xmlParserCtxt, decltype(free_parser)> parser(xmlCreatePushParserCtxt(&sax, nullptr, nullptr, 0, nullptr), free_parser);
      if (!parser) { throw std::bad_alloc(); }
      parser->_private = &state;
      parser->userData = parser.get();
      xmlCtxtSetErrorHandler(parser.get(), LibxmlSAX2State::error, parser.get());
      // Same options and chunked feeding as the proof of concept.
      const int options = XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT;
      xmlCtxtUseOptions(parser.get(), options);
      // libxml2-fast: decode Latin-1 in the adapter (see mode description above).
      bool latin1 = false;
      if (source->getEncoding())
      {
        const auto encoding = StringManager::convert(reinterpret_cast<const char16_t*>(source->getEncoding()));
        if (fast && isLatin1Name_(encoding))
        {
          latin1 = true;
          xmlCtxtUseOptions(parser.get(), options | XML_PARSE_IGNORE_ENC);
        }
        else if (xmlSwitchEncodingName(parser.get(), encoding.c_str()) != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Unsupported XML encoding: " + encoding);
        }
      }
      const bool sniff_declaration = fast && !source->getEncoding();
      unique_ptr<xercesc::BinInputStream> stream(source->makeStream());
      if (!stream) { throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Cannot open XML input stream"); }
      char buffer[65536];
      std::string converted;
      bool first = true;
      int status = 0;
      while (!state.exception)
      {
        const auto count = stream->readBytes(reinterpret_cast<XMLByte*>(buffer), sizeof(buffer));
        if (first && sniff_declaration && isLatin1Name_(declaredEncoding_(buffer, count)))
        {
          latin1 = true;
          xmlCtxtUseOptions(parser.get(), options | XML_PARSE_IGNORE_ENC);
        }
        first = false;
        const char* data = buffer;
        std::size_t size = count;
        if (latin1 && count > 0 && !isASCII_(buffer, count))
        {
          latin1ToUTF8_(buffer, count, converted);
          data = converted.data();
          size = converted.size();
        }
        status = xmlParseChunk(parser.get(), data, static_cast<int>(size), count == 0);
        if (status || count == 0) { break; }
      }
      if (state.exception)
      {
        try { std::rethrow_exception(state.exception); }
        catch (const XMLHandler::EndParsingSoftly&) { return; }
      }
      if (status)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "libxml2 rejected XML input");
      }
    }

    void parse(xercesc::InputSource* const source, XMLHandler* handler)
    {
      switch (xmlAdapterMode_())
      {
        case XMLAdapterMode::Xerces: parseXerces_(source, handler); break;
        case XMLAdapterMode::LibxmlPOC: parseLibxmlPOC_(source, handler); break;
        case XMLAdapterMode::LibxmlSAX2: parseLibxmlSAX2_(source, handler, false); break;
        case XMLAdapterMode::LibxmlUTF8: parseLibxmlSAX2_(source, handler, true); break;
        case XMLAdapterMode::LibxmlFast: parseLibxmlSAX2_(source, handler, true, true); break;
      }
    }

    void XMLFile::parse_(const std::string & filename, XMLHandler * handler)
    {
      // ensure handler->reset() is called to save memory (in case the XMLFile
      // reader, e.g. FeatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      StringManager sm;
      //try to open file
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      // initialize parser
      try
      {
        Internal::xercesInitialize(); // serialised: Xerces requires it, and parseBuffer_() runs on several threads at once
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "",std::string("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }


      // peak ahead into the file: is it bzip2 or gzip compressed?
      std::string bz;
      {
        std::ifstream file(filename.c_str());
        char tmp_bz[3];
        file.read(tmp_bz, 2);
        tmp_bz[2] = '\0';
        bz =std::string(tmp_bz);
      }

      unique_ptr<xercesc::InputSource> source;

      char g1 = 0x1f;
      char g2 = 0;
      g2 |= 1 << 7;
      g2 |= 1 << 3;
      g2 |= 1 << 1;
      g2 |= 1 << 0;
      //g2 = static_cast<char>(0x8b); // can make troubles if it is casted to 0x7F which is the biggest number signed char can save
      if ((bz[0] == 'B' && bz[1] == 'Z') || (bz[0] == g1 && bz[1] == g2) || (bz[0] == 'P' && bz[1] == 'K'))
      {
        source.reset(new CompressedInputSource(sm.convert(filename).c_str(), bz));
      }
      else
      {
        source.reset(new xercesc::LocalFileInputSource(sm.convert(filename).c_str()));
      }
      // what if no encoding given http://xerces.apache.org/xerces-c/apiDocs-3/classInputSource.html
      if (!enforced_encoding_.empty())
      {
        const auto encoding = sm.convert(enforced_encoding_);
        source->setEncoding(reinterpret_cast<const XMLCh*>(encoding.c_str()));
      }
      
      parse(source.get(), handler);
    }

    void XMLFile::parseBuffer_(const std::string & buffer, XMLHandler * handler)
    {
      // ensure handler->reset() is called to save memory (in case the XMLFile
      // reader, e.g. FeatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      StringManager sm;

      // initialize parser
      try
      {
        Internal::xercesInitialize(); // serialised: Xerces requires it, and parseBuffer_() runs on several threads at once
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "",std::string("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }

      // TODO: handle non-plain text
      // peak ahead into the file: is it bzip2 or gzip compressed?
      // std::string bz = StringUtils::substr(buffer, 0, 2);

      unique_ptr<xercesc::InputSource> source;
      {
        auto fake_id = sm.convert("inMemory");
        source.reset(new xercesc::MemBufInputSource(reinterpret_cast<const unsigned char *>(buffer.c_str()), buffer.size(), fake_id.c_str()));
      }
      // what if no encoding given http://xerces.apache.org/xerces-c/apiDocs-3/classInputSource.html
      if (!enforced_encoding_.empty())
      {
        const auto encoding = sm.convert(enforced_encoding_);
        source->setEncoding(reinterpret_cast<const XMLCh*>(encoding.c_str()));
      }
      
      parse(source.get(), handler);
    }

    void XMLFile::save_(const std::string & filename, XMLHandler * handler) const
    {
      // Detect compression from filename extension, in any letter case like on input (and in TOPPBase's check)
      const FileTypes::Type compression = FileNameUtils::compressionType(filename);
      const bool use_gzip = compression == FileTypes::GZ;
      const bool use_bzip2 = compression == FileTypes::BZ2;

      if (use_gzip)
      {
        // Minimal std::streambuf that writes gzip-compressed data via zlib
        class GzipOStreambuf : public std::streambuf
        {
        public:
          enum { BUF_SIZE = 65536 };

          explicit GzipOStreambuf(const char* fn) : gz_(gzopen(fn, "wb"))
          {
            if (gz_) setp(buf_, buf_ + BUF_SIZE);
          }

          ~GzipOStreambuf()
          {
            flush_buffer();
            if (gz_) gzclose(gz_);
          }

          bool is_open() const { return gz_ != nullptr; }

        protected:
          int_type overflow(int_type c) override
          {
            if (!gz_) return traits_type::eof();
            if (pptr() == epptr() && !flush_buffer()) return traits_type::eof();
            if (!traits_type::eq_int_type(c, traits_type::eof()))
            {
              *pptr() = traits_type::to_char_type(c);
              pbump(1);
            }
            return traits_type::not_eof(c);
          }

          int sync() override { return flush_buffer() ? 0 : -1; }

        private:
          bool flush_buffer()
          {
            int n = static_cast<int>(pptr() - pbase());
            if (n > 0)
            {
              if (gzwrite(gz_, pbase(), static_cast<unsigned>(n)) != n) return false;
              setp(buf_, buf_ + BUF_SIZE);
            }
            return true;
          }

          gzFile gz_;
          char buf_[BUF_SIZE];
        };

        GzipOStreambuf sbuf(filename.c_str());
        if (!sbuf.is_open())
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        std::ostream os(&sbuf);
        os.precision(writtenDigits(double()));
        handler->writeTo(os);
      }
      else if (use_bzip2)
      {
        // Minimal std::streambuf that writes bzip2-compressed data via bzlib
        class Bzip2OStreambuf : public std::streambuf
        {
        public:
          enum { BUF_SIZE = 65536 };

          explicit Bzip2OStreambuf(const char* fn) : file_(nullptr), bz_(nullptr)
          {
            file_ = fopen(fn, "wb");
            if (!file_) return;
            int err = BZ_OK;
            bz_ = BZ2_bzWriteOpen(&err, file_, 9, 0, 0);
            if (err != BZ_OK)
            {
              BZ2_bzWriteClose(&err, bz_, 1, nullptr, nullptr);
              bz_ = nullptr;
              fclose(file_);
              file_ = nullptr;
              return;
            }
            setp(buf_, buf_ + BUF_SIZE);
          }

          ~Bzip2OStreambuf()
          {
            flush_buffer();
            if (bz_)
            {
              int err = BZ_OK;
              BZ2_bzWriteClose(&err, bz_, 0, nullptr, nullptr);
            }
            if (file_) fclose(file_);
          }

          bool is_open() const { return bz_ != nullptr; }

        protected:
          int_type overflow(int_type c) override
          {
            if (!bz_) return traits_type::eof();
            if (pptr() == epptr() && !flush_buffer()) return traits_type::eof();
            if (!traits_type::eq_int_type(c, traits_type::eof()))
            {
              *pptr() = traits_type::to_char_type(c);
              pbump(1);
            }
            return traits_type::not_eof(c);
          }

          int sync() override { return flush_buffer() ? 0 : -1; }

        private:
          bool flush_buffer()
          {
            if (!bz_) return false;
            int n = static_cast<int>(pptr() - pbase());
            if (n > 0)
            {
              int err = BZ_OK;
              BZ2_bzWrite(&err, bz_, pbase(), n);
              if (err != BZ_OK && err != BZ_RUN_OK) return false;
              setp(buf_, buf_ + BUF_SIZE);
            }
            return true;
          }

          FILE* file_;
          BZFILE* bz_;
          char buf_[BUF_SIZE];
        };

        Bzip2OStreambuf sbuf(filename.c_str());
        if (!sbuf.is_open())
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        std::ostream os(&sbuf);
        os.precision(writtenDigits(double()));
        handler->writeTo(os);
      }
      else
      {
        // Uncompressed: open in binary mode to avoid any line ending conversions
        std::ofstream os(filename.c_str(), std::ios::out | std::ios::binary);
        os.precision(writtenDigits(double()));
        if (!os)
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        handler->writeTo(os);
        os.close();
      }
    }

    std::string encodeTab(const std::string& to_encode)
    {
      if (!StringUtils::has(to_encode, '\t'))
      {
        return to_encode;
      }
      else
      {
        std::string result = to_encode;
        StringUtils::substitute(result, "\t", "&#x9;");
        return result;
      }
    }

    bool XMLFile::isValid(const std::string & filename, std::ostream & os)
    {
      if (schema_location_.empty())
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      std::string current_location = File::find(schema_location_);
      return XMLValidator().isValid(filename, current_location, os);
    }

    const std::string & XMLFile::getVersion() const
    {
      return schema_version_;
    }

} // namespace OpenMS  // namespace Internal
