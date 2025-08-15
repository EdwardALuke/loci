//#############################################################################
//#
//# Copyright 2025, Loci Consulting Services, LLC
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
#ifndef DOCREPORT_H
#define DOCREPORT_H

#include <rule.h>
#include <iostream>

namespace Loci {
  /** 
   * @brief Function that will extract default and optional rules and output
   * a description of the inputs by extracting the comments from the rules.
   * 
   * @param rdb  [-] Input Rule Database
   * @param o    [-] Output stream to write the summary (usually cout)
   */
  void describeInputs(const rule_db &rdb, std::ostream &o) ;
  /**
   * @brief Convenience function that defaults to writing to cout
   *
   * @param rdb  [-] Input Rule Database
   */
  inline void describeInputs(const rule_db &rdb) { describeInputs(rdb,std::cout) ; }
  /** 
   * @brief Function that will extract and collate documentation on the rules
   * in the provided rule database that compute the provided variables.  If 
   * the variable set input is empty then the function will write out 
   * documentation for all of the rules.
   *
   * @param rdb  [-] Input Rule Database
   * @param vars [-] Set of variables to document
   * @param o    [-] Output stream to write documentation (usually cout)
   */
  void ruleDocumentation(const rule_db &rdb, variableSet vars,
                         std::ostream &o) ;
  /** 
   * @brief Convienience function that defaults to writing rule documentation 
   * to cout based on input rule database and requested variables.
   *
   * @param rdb  [-] Input Rule Database
   * @param vars [-] Set of variables to document
   */
  inline void ruleDocumentation(const rule_db &rdb, variableSet vars) {
    ruleDocumentation(rdb,vars,std::cout) ;
  }
}

#endif
