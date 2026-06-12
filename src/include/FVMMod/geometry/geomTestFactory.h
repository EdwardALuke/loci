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

#ifndef GEOM_TEST_FACTORY_H
#define GEOM_TEST_FACTORY_H

#include "leftPlane.h"
#include "inBox.h"
#include "inSphere.h"
#include "inCylinder.h"
#include "inCone.h"
#include "nearRectangle.h"
/** ****************************************************************************
 * @brief Function to check that the user input options are appropriate for the
 *        geomTest class.
 * @param name Name of user input (and associated class) geometry.
 * @param ol   List of options provided by user for the specific geometry class.
 * @return Loci::CPTR<geomTest> Pointer to newly created object.
 ******************************************************************************/
inline Loci::CPTR<geomTest> geomTestFactory(string name, const options_list ol)
{
  Loci::CPTR<geomTest> gp = 0;
  if(     name == "sphere")     gp = new inSphere(ol);
  else if(name == "box")        gp = new inBox(ol);
  else if(name == "cylinder")   gp = new inCylinder(ol);
  else if(name == "cone")       gp = new inCone(ol);
  else if(name == "leftPlane")  gp = new leftPlane(ol);
  else if(name == "hex")        gp = new hex(ol);
  else if(name == "") 
  else
    printf("ERROR::geomTestFactory: Don't know what to do with '%s'\n", name.c_str());

  return gp;
} // End of geomTestFactory()


/** ****************************************************************************
 * @brief Function to check that the user input options are appropriate for the
 *        geomTest class.
 * @param name Name of user input (and associated class) geometry.
 * @param ol   List of options provided by user for the specific geometry class.
 * @return Loci::CPTR<geomTest> Pointer to newly created object.
 ******************************************************************************/
inline Loci::CPTR<nearRectangle> geomAreaTest(string name, const options_list ol)
{
  int                       iError = 0;
  Loci::CPTR<nearRectangle> gp     = 0;

  if(name == "nearRectangle") gp     = new nearRectangle(ol);
  else                        iError = 1;

  if(iError)
  {
    printf("ERROR::geomAreaTest: '%s' is not a recognized area class.\n", name.c_str());
    Loci::Abort();
  }

  return gp;
} // End of geomAreaTest()

#endif // GEOM_TEST_FACTORY_H
