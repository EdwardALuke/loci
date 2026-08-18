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

#include "template_engine.h"
#include "template_parser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace Loci {

  void TemplateEngine::define(
    std::string const & name, std::string const & text
  ) {
    try {
      TemplateParser parser(text) ;
      templates_[name] = parser.parse() ;
    } catch(std::exception const & ex) {
      throw std::runtime_error(
        "template \"" + name + "\": " + ex.what()
      ) ;
    }
  }

  void TemplateEngine::define_from_file(
    std::string const & name, std::string const & path
  ) {
    std::ifstream in(path, std::ios::binary) ;
    if(!in) {
      throw std::runtime_error(
        "template \"" + name + "\": cannot read template file: " + path
      ) ;
    }
    std::ostringstream ss ;
    ss << in.rdbuf() ;
    define(name, ss.str()) ;
  }

  TemplateNodeList const * TemplateEngine::find(
    std::string const & name
  ) const {
    std::map<std::string, TemplateNodeList>::const_iterator iter =
      templates_.find(name) ;
    if(iter == templates_.end()) {
      return nullptr ;
    }
    return &iter->second ;
  }

  std::string TemplateEngine::render(
    std::string const & name, TemplateValue const & context
  ) const {
    return render(name, context, TemplateNewlineHandler{}) ;
  }

  std::string TemplateEngine::render(
    std::string const & name,
    TemplateValue const & context,
    TemplateNewlineHandler newline_handler
  ) const {
    TemplateNodeList const * nodes = find(name) ;
    if(!nodes) {
      throw std::runtime_error("undefined template: " + name) ;
    }

    std::ostringstream ss ;
    TemplateScope scope(context, *this) ;
    if(newline_handler) {
      scope.newline_handler = &newline_handler ;
    }
    TemplateNodeList::const_iterator niter = nodes->begin() ;
    while(niter != nodes->end()) {
      (*niter)->evaluate(ss, scope) ;
      ++niter ;
    }
    return ss.str() ;
  }

} // end: namespace Loci
