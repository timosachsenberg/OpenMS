// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <vector>

#include <QtGui/QColor>

namespace OpenMS
{
  /** @brief An annotation item which represents a vertical line.
          @see Annotation1DItem
  */
  class Annotation1DVerticalLineItem :
      public Annotation1DItem
  {

  public:
    /// Constructor
    Annotation1DVerticalLineItem(const double& x, const QColor& color, const QString & tex="");
    /// Copy constructor
    Annotation1DVerticalLineItem(const Annotation1DVerticalLineItem & rhs);
    /// Destructor
    ~Annotation1DVerticalLineItem() override;
    // Docu in base class
    void ensureWithinDataRange(Spectrum1DCanvas * const canvas) override;
    // Docu in base class
    void draw(Spectrum1DCanvas * const canvas, QPainter & painter, bool flipped = false) override;
    // Docu in base class
    void move(const PointType & delta) override;

    /// Sets the uppermost position of the line
    void setPosition(const double & x);
    /// Returns the position
    const double & getPosition() const;

  protected:
    /// The position of the vertical line
    double x_;

    /// The color of the line
    QColor color_;

  };
} // namespace OpenMS
