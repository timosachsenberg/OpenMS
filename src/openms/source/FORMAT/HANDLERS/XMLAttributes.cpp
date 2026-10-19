// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>

#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/TransService.hpp>
#include <cstring>

namespace OpenMS::Internal
{
  // XMLCh is guaranteed to match char16_t in size; the attribute values are
  // consumed purely as UTF-16 code units, so the reinterpret_cast at this
  // boundary is well-defined for our use.
  static_assert(sizeof(::XMLCh) == sizeof(char16_t), "XMLCh is not sized correctly for UTF-16.");

  namespace
  {
    inline const xercesc::Attributes* attrs(const void* handle)
    {
      return static_cast<const xercesc::Attributes*>(handle);
    }

    bool equalsASCII(const char16_t* wide, const char* ascii)
    {
      while (*ascii && *wide == static_cast<unsigned char>(*ascii))
      {
        ++wide;
        ++ascii;
      }
      return *ascii == 0 && *wide == 0;
    }

    bool equals16(const char16_t* a, const char16_t* b)
    {
      while (*a && *a == *b)
      {
        ++a;
        ++b;
      }
      return *a == 0 && *b == 0;
    }
  }

  void XMLAttributes::assignUTF16(std::u16string& out, const char* in, Size length)
  {
    // UTF-16 never needs more code units than UTF-8 needs bytes. The buffer
    // keeps its capacity across calls, so steady-state use does not allocate.
    out.resize(length);
    const unsigned char* cursor = reinterpret_cast<const unsigned char*>(in);
    const unsigned char* const end = cursor + length;
    char16_t* target = out.data();
    while (cursor < end && *cursor < 0x80) { *target++ = static_cast<char16_t>(*cursor++); }
    while (cursor < end)
    {
      // Same decoding as the proof-of-concept adapter; the parser has already validated UTF-8.
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
        *target++ = static_cast<char16_t>(0xd800 + (cp >> 10));
        *target++ = static_cast<char16_t>(0xdc00 + (cp & 0x3ff));
      }
      else { *target++ = static_cast<char16_t>(cp); }
    }
    out.resize(static_cast<Size>(target - out.data()));
  }

  SignedSize XMLAttributes::nativeIndex_(const char16_t* qname) const
  {
    for (Size i = 0; i < native_view_->count; ++i)
    {
      if (equals16(native_view_->names16[i]->c_str(), qname)) { return static_cast<SignedSize>(i); }
    }
    return -1;
  }

  const char16_t* XMLAttributes::nativeValue16_(Size i) const
  {
    const NativeView& v = *native_view_;
    if (v.lazy_utf16 && !v.values16_ready[i])
    {
      assignUTF16(v.values16[i], v.values8[i], v.values8_length[i]);
      v.values16_ready[i] = 1;
    }
    return v.values16[i].c_str();
  }

  const char* XMLAttributes::valueUTF8(const char16_t* qname, Size& length) const
  {
    if (!hasUTF8()) { return nullptr; }
    const SignedSize i = nativeIndex_(qname);
    if (i < 0) { return nullptr; }
    length = native_view_->values8_length[i];
    return native_view_->values8[i];
  }

  const char16_t* XMLAttributes::value(const char16_t* qname) const
  {
    if (native_view_)
    {
      const SignedSize i = nativeIndex_(qname);
      return i < 0 ? nullptr : nativeValue16_(static_cast<Size>(i));
    }
    if (native_attributes_)
    {
      const SignedSize i = index(qname);
      return i < 0 ? nullptr : (*native_attributes_)[i].second.c_str();
    }
    return reinterpret_cast<const char16_t*>(attrs(attributes_)->getValue(reinterpret_cast<const ::XMLCh*>(qname)));
  }

  const char16_t* XMLAttributes::value(const char* qname) const
  {
    if (native_view_)
    {
      for (Size i = 0; i < native_view_->count; ++i)
      {
        if (std::strcmp(native_view_->names8[i], qname) == 0) { return nativeValue16_(i); }
      }
      return nullptr;
    }
    if (native_attributes_)
    {
      bool ascii = true;
      for (const auto* ch = qname; *ch; ++ch)
      {
        if (static_cast<unsigned char>(*ch) > 127) { ascii = false; break; }
      }
      if (ascii)
      {
        // XML vocabulary names are usually ASCII; avoid allocating a transcoder.
        for (const auto& item : *native_attributes_)
        {
          if (equalsASCII(item.first.c_str(), qname)) return item.second.c_str();
        }
        return nullptr;
      }
      const xercesc::TranscodeFromStr transcoded(reinterpret_cast<const XMLByte*>(qname), std::strlen(qname), "UTF-8");
      return value(reinterpret_cast<const char16_t*>(transcoded.str()));
    }
    ::XMLCh* transcoded = xercesc::XMLString::transcode(qname);
    const char16_t* val = reinterpret_cast<const char16_t*>(attrs(attributes_)->getValue(transcoded));
    xercesc::XMLString::release(&transcoded);
    return val;
  }

  SignedSize XMLAttributes::index(const char16_t* qname) const
  {
    if (native_view_) { return nativeIndex_(qname); }
    if (native_attributes_)
    {
      for (Size i = 0; i < native_attributes_->size(); ++i)
      {
        if ((*native_attributes_)[i].first == qname) { return static_cast<SignedSize>(i); }
      }
      return -1;
    }
    return static_cast<SignedSize>(attrs(attributes_)->getIndex(reinterpret_cast<const ::XMLCh*>(qname)));
  }

  const char16_t* XMLAttributes::valueByIndex(Size i) const
  {
    if (native_view_) { return i < native_view_->count ? nativeValue16_(i) : nullptr; }
    if (native_attributes_) { return i < native_attributes_->size() ? (*native_attributes_)[i].second.c_str() : nullptr; }
    return reinterpret_cast<const char16_t*>(attrs(attributes_)->getValue(static_cast<XMLSize_t>(i)));
  }

} // namespace OpenMS::Internal
