// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS::Internal
{
  /**
    @brief Xerces-free, lightweight view over the attribute list of an XML element.

    Wraps a live parser attribute list (owned by the SAX adapter and valid only
    for the duration of the @c startElement callback). It exposes only the small
    surface the OpenMS handlers actually use (value-by-name, index-by-name,
    value-by-index), so that Xerces types never appear in a public header.

    The underlying attribute list is held as an opaque pointer, so this header
    pulls in no Xerces declarations — consumers of libOpenMS need Xerces neither
    on their include path nor at link time. The accessor implementations
    (in @c XMLAttributes.cpp) recover the concrete Xerces type internally.

    The wrapper is a non-owning view: it stores a pointer to the parser-owned
    attribute list and must not outlive the callback it was created for.
  */
  class OPENMS_DLLAPI XMLAttributes
  {
  public:
    XMLAttributes(const XMLAttributes&) = default;
    XMLAttributes& operator=(const XMLAttributes&) = default;

    /**
      @brief Wrap a live attribute list (non-owning).

      Takes an opaque handle so this header need not name any Xerces type. The
      only caller is the internal SAX2HandlerAdapter, which passes the address of
      the live @c xercesc::Attributes (recovered internally in @c XMLAttributes.cpp).
    */
    explicit XMLAttributes(const void* attributes) noexcept :
      attributes_(attributes)
    {
    }

    /// Non-owning UTF-16 attribute view for the libxml2 proof of concept.
    explicit XMLAttributes(const std::vector<std::pair<std::u16string, std::u16string>>& attributes) noexcept :
      attributes_(nullptr), native_attributes_(&attributes)
    {
    }

    /**
      @brief Adapter-owned attribute storage of the reworked libxml2 SAX2 adapter (benchmark prototype).

      All arrays hold @c count entries and are valid only during the startElement
      callback. UTF-8 names are NUL-terminated, UTF-8 values are not. The adapter
      reuses all buffers across elements, so a steady-state element needs no
      heap allocation. With @c lazy_utf16 set, UTF-16 values are converted on
      first access only.
    */
    struct NativeView
    {
      Size count = 0;
      const char* const* names8 = nullptr;              ///< UTF-8 qualified names
      const char* const* values8 = nullptr;             ///< UTF-8 values (not NUL-terminated)
      const Size* values8_length = nullptr;             ///< value lengths in bytes
      const std::u16string* const* names16 = nullptr;   ///< UTF-16 qualified names (adapter name cache)
      std::u16string* values16 = nullptr;               ///< UTF-16 value buffers (reused)
      unsigned char* values16_ready = nullptr;          ///< conversion flags (lazy mode)
      bool lazy_utf16 = false;                          ///< convert UTF-16 values on first access
      bool expose_utf8 = false;                         ///< let handlers read UTF-8 values directly
    };

    /// Non-owning view over the reworked adapter's attribute storage.
    explicit XMLAttributes(const NativeView& view) noexcept :
      attributes_(nullptr), native_view_(&view)
    {
    }

    /// True if values can be read as UTF-8 without conversion (see valueUTF8()).
    bool hasUTF8() const noexcept { return native_view_ != nullptr && native_view_->expose_utf8; }

    /**
      @brief UTF-8 value of the attribute named @p qname, or @c nullptr if absent.

      The result is not NUL-terminated; @p length receives its size in bytes.
      Only meaningful if hasUTF8() returns true (otherwise always @c nullptr).
    */
    const char* valueUTF8(const char16_t* qname, Size& length) const;

    /// Replace @p out by the UTF-16 form of the (valid) UTF-8 sequence @p in of @p length bytes.
    static void assignUTF16(std::u16string& out, const char* in, Size length);

    /// Value of the attribute named @p qname (UTF-16 code units), or @c nullptr if absent.
    const char16_t* value(const char16_t* qname) const;

    /// Value of the attribute named @p qname (narrow name, transcoded internally), or @c nullptr if absent.
    const char16_t* value(const char* qname) const;

    /// Index of the attribute named @p qname, or a negative value if absent.
    SignedSize index(const char16_t* qname) const;

    /// Value of the attribute at position @p i (UTF-16 code units).
    const char16_t* valueByIndex(Size i) const;

    /// Opaque handle to the underlying attribute list (libOpenMS-internal use).
    const void* handle() const noexcept { return attributes_; }

  private:
    /// Index of @p qname in the native view, or -1.
    SignedSize nativeIndex_(const char16_t* qname) const;
    /// UTF-16 value at position @p i of the native view (converted on demand in lazy mode).
    const char16_t* nativeValue16_(Size i) const;

    const void* attributes_;
    const std::vector<std::pair<std::u16string, std::u16string>>* native_attributes_ = nullptr;
    const NativeView* native_view_ = nullptr;
  };

} // namespace OpenMS::Internal
