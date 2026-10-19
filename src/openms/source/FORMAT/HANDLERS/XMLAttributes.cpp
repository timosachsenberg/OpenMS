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
  }

  const char16_t* XMLAttributes::value(const char16_t* qname) const
  {
    if (native_attributes_)
    {
      const SignedSize i = index(qname);
      return i < 0 ? nullptr : (*native_attributes_)[i].second.c_str();
    }
    return reinterpret_cast<const char16_t*>(attrs(attributes_)->getValue(reinterpret_cast<const ::XMLCh*>(qname)));
  }

  const char16_t* XMLAttributes::value(const char* qname) const
  {
    if (native_attributes_)
    {
      std::u16string name;
      while (*qname) { name.push_back(static_cast<unsigned char>(*qname++)); }
      return value(name.c_str());
    }
    ::XMLCh* transcoded = xercesc::XMLString::transcode(qname);
    const ::XMLCh* val = attrs(attributes_)->getValue(transcoded);
    xercesc::XMLString::release(&transcoded);
    return reinterpret_cast<const char16_t*>(val);
  }

  SignedSize XMLAttributes::index(const char16_t* qname) const
  {
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
    if (native_attributes_) { return i < native_attributes_->size() ? (*native_attributes_)[i].second.c_str() : nullptr; }
    return reinterpret_cast<const char16_t*>(attrs(attributes_)->getValue(static_cast<XMLSize_t>(i)));
  }

} // namespace OpenMS::Internal
