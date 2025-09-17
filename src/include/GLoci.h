//#############################################################################
//#
//# Copyright 2008-2025,  Mississippi State University
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
#ifndef GLOCI_NOH
#define GLOCI_NOH

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <string>
#include <Loci_version.h>
#include <Loci_types.h>
#include <data_traits.h>
#include <gmap.h>
#include <gmapvec.h>
#include <gmultimap.h>
#include <gstore.h>
//#include <storeVec.h>
#include <gparameter.h>
#include <gblackbox.h>
#include <gconstraint.h>
//#include <rule.h>
//#include <Tools/options_list.h>
#include <Tools/fpe.h>
#include <scheduler.h>
#include <distribute.h>
#include <distribute_io.h>
#include <distribute_container.h>
#include <gmultistore.h>
#include <mod_db.h>
#include <gstore.h>
#include <fact_db.h>
#include <gkeyspace.h>
#include <gkey_manager.h>
//#include <blackbox.h>

#include <LociGridReaders.h>
#include <pnn.h>

#define HAS_MODULE_SEARCH_DIR
namespace Loci {
  extern void AddModuleSearchDir(std::string dirname);
  extern double random() ;
  extern int irandom() ;
}
using Loci::fact_db;
using Loci::gParam;
using Loci::gConstraint;
using Loci::gMap;
using Loci::gStore;
using Loci::gMultiMap;
using Loci::gMapVec;
//using Loci::gStoreVec;
using Loci::gMultiStore;
using Loci::gKeySpace;
using Loci::gKeySpaceP;
using Loci::gEntitySet;
using Loci::GEMPTY;
#endif





