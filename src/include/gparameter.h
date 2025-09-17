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

//template class paramRepI and gParamRepI are circular referenced
//in order to resolve the compiling problem, for both classes,
//declaration are separated from definition,
//and whereever either class is used, both declaration files are included first
//followed by the both definition files

#ifndef GPARAMETER_H
#define GPARAMETER_H
#include <parameter.h>
#include <gparameter_def.h>
#include <gparameter_impl.h>
#endif
