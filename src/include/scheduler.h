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

  /// @brief This function is used to create an execution schedule for
  ///   `internalQuery()`.
  ///
  /// This function and the `internalQuery()` function are mainly intended to
  /// be used by the Loci scheduler to compute intermediate values.
  ///
  /// NOTE: the passed in rule database is the expanded parametric rule
  /// database since the expanded rule base is what we needed and this
  /// expansion process only needs to be performed once.
  ///
  /// @param[in] par_rdb Expanded parametric rule database used to build the
  ///   internal schedule.
  /// @param[in,out] facts Fact database used while constructing the internal
  ///   execution schedule.
  /// @param[in,out] scheds Scheduler database populated while constructing the
  ///   internal execution schedule.
  /// @param[in] target Variables whose internal execution schedule is being
  ///   created.
  /// @param[in] nth Currently unused legacy parameter retained for API
  ///   compatibility.
  /// @return Execution schedule for the requested internal query, or a null
  ///   schedule when no target vertices are reachable.
  executeP create_internal_execution_schedule(rule_db& par_rdb,
                                              fact_db &facts,
                                              sched_db &scheds,
                                              const variableSet& target,
                                              int nth=1) ;

  /// @brief This function is used by the Loci scheduler to issue queries for
  ///   intermediate relations.
  ///
  /// It is basically a reduced version of the user function `makeQuery()`.
  ///
  /// NOTE: the passed in rule database is the expanded parametric rule
  /// database since the expanded rule base is what we needed and this
  /// expansion process only needs to be performed once.
  ///
  /// @param[in] par_rdb Expanded parametric rule database used to satisfy the
  ///   internal query.
  /// @param[in,out] facts Fact database read during query planning and updated
  ///   with the queried intensional facts on success.
  /// @param[in] query Intermediate facts requested by the scheduler.
  /// @return True when the internal query produces an execution schedule and
  ///   publishes the queried facts, or false when no schedule can be created.
  bool internalQuery(rule_db& par_rdb, fact_db& facts,
                     const variableSet& query) ;


  /// @brief Creates an execution schedule that derives the target facts from
  ///   the given facts using `rdb`.
  ///
  /// This routine expands the parametric rule database for the requested
  /// target facts, generates the dependency graph, decomposes that graph,
  /// performs existential analysis, and then creates the execution schedule.
  /// The returned schedule is constructed here but not executed here.
  ///
  /// Schedule construction updates both the fact database and the schedule
  /// database. In particular, distributed facts may be localized/frozen for
  /// scheduling, and `scheds` is initialized and populated for the requested
  /// schedule.
  ///
  /// @param[in] rdb Rule database used to derive the target facts.
  /// @param[in,out] facts Fact database containing the current given facts and
  ///   updated as needed during schedule construction.
  /// @param[in,out] scheds Schedule database initialized from `facts` and
  ///   populated for the requested execution schedule.
  /// @param[in] target Target facts to be derived.
  /// @param[in] nth Currently unused legacy parameter retained for API
  ///   compatibility.
  /// @return Execution schedule for the requested target facts, or a null
  ///   schedule when no reachable dependency graph can be constructed.
  extern executeP create_execution_schedule(const rule_db &rdb,
                                            fact_db &facts,
                                            sched_db &scheds,
                                            const variableSet& target,
                                            int nth=1) ;

  /// @brief Issues a user query that derives the queried facts from the given
  ///   known facts using `rdb`.
  ///
  /// This is the user-facing query entry point. It installs default facts,
  /// parses `query` into the target queried facts, ignores queried facts that
  /// are already extensional, and for the remaining targets clears prior
  /// intensional facts, creates and executes an execution schedule, and stores
  /// the queried results back into `facts` as intensional facts.
  ///
  /// Queried results are therefore maintained in `facts` until the next
  /// non-extensional query.
  ///
  /// NOTE: in the current implementation, scheduling or execution failures
  /// abort rather than returning `false`.
  ///
  /// @param[in] rdb Rule database used to derive the queried facts.
  /// @param[in,out] facts Fact database holding the given known facts and
  ///   updated with default facts and queried intensional results.
  /// @param[in] query Query expression naming the target queried facts.
  /// @return True when the queried facts are already extensional or are
  ///   derived successfully.
  extern bool makeQuery(const rule_db &rdb, fact_db &facts,
                        const std::string& query) ;
}

#endif
