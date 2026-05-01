/** ****************************************************************************
 * @file      geomTestFactory.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief     Contains the geometry test prototypes.
 * @date      2022-11-08
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/

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
  if(     name == "inSphere")   gp = new inSphere(ol);
  else if(name == "inBox")      gp = new inBox(ol);
  else if(name == "inCylinder") gp = new inCylinder(ol);
  else if(name == "inCone")     gp = new inCone(ol);
  else if(name == "leftPlane")  gp = new leftPlane(ol);
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
