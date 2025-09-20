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
#ifndef GLOBALS_H
#define GLOBALS_H
#include "defines.h"
class Globals{
public:
  static double fold;
  static double tolerance;
  static int levels;
  static double factor;
  //balance options:
  //0: no edge's depth is greater than 1
  //1: 0 and no cell has more than half of its face split
  //2: 0 and 1 and no cell has two opposite faces split
  static int balance_option;
  static vect3d split;
  static vect3d nosplit;
};

#endif
