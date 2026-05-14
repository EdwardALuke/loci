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

#ifndef SCHEDULER_H
#define SCHEDULER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <rule.h>
#include <fact_db.h>
#include <sched_db.h>
#include <execute.h>


namespace Loci {

  // function prototype for internal queries
  executeP create_internal_execution_schedule(rule_db& par_rdb,
                                              fact_db &facts,
                                              sched_db &scheds,
                                              const variableSet& target,
                                              int nth=1) ;
  bool internalQuery(rule_db& par_rdb, fact_db& facts,
                     const variableSet& query) ;
    

  extern executeP create_execution_schedule(const rule_db &rdb,
                                            fact_db &facts,
                                            sched_db &scheds,
                                            const variableSet& target,
                                            int nth=1) ;

  extern bool makeQuery(const rule_db &rdb, fact_db &facts,
                        const std::string& target) ;
}

#endif
