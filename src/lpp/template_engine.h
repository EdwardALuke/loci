//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef TEMPLATE_ENGINE_H
#define TEMPLATE_ENGINE_H

#include "template_ast.h"

#include <map>
#include <string>

namespace Loci {

  class TemplateEngine {
    std::map<std::string, TemplateNodeList> templates_ ;

  public:
    static constexpr int kMaxPartialDepth = 64 ;

    void define(std::string const & name, std::string const & text) ;
    void define_from_file(std::string const & name, std::string const & path) ;
    TemplateNodeList const * find(std::string const & name) const ;
    std::string render(
      std::string const & name, TemplateValue const & context
    ) const ;
    std::string render(
      std::string const & name,
      TemplateValue const & context,
      TemplateNewlineHandler newline_handler
    ) const ;
  } ;

} // end: namespace Loci

#endif
