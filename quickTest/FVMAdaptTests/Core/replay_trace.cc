//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
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

#include "fvmadapt_core.h"

namespace {

/// Cell-family selector for plan replay.
enum TraceKind {
  TRACE_HEX,
  TRACE_PRISM,
  TRACE_GENERAL
} ;

/// Queue state for one node in the replayed cell tree.
struct TraceState {
  int id ;
  int parent_id ;
  int child_slot ;
  int depth ;
  int local_fold ;
  std::string path ;
} ;

/// Output row for one breadth-first replay step.
struct TraceRow {
  int id ;
  int parent_id ;
  int child_slot ;
  int depth ;
  int local_fold ;
  std::string path ;
  int plan_index ;
  bool explicit_code ;
  int code ;
  int child_count ;
  int first_child_id ;
  int leaf_id ;
} ;

/// Complete replay result for one cell plan.
struct PlanReplayTrace {
  std::vector<TraceRow> rows ;
  int consumed_plan_entries ;
  int leaf_count ;
} ;

/// Format the unused suffix of a plan vector.
std::string plan_tail_string(const std::vector<char>& plan, int start) {
  if(start >= int(plan.size())) {
    return "[]" ;
  }
  std::vector<char> tail ;
  for(size_t i = size_t(start); i < plan.size(); ++i) {
    tail.push_back(plan[i]) ;
  }
  return plan_string(tail) ;
}

/// Return the number of children created by one HexCell split code.
int hex_trace_child_count(int code) {
  switch(code) {
  case 0:
    return 0 ;
  case 1:
  case 2:
  case 4:
    return 2 ;
  case 3:
  case 5:
  case 6:
    return 4 ;
  case 7:
    return 8 ;
  default:
    break ;
  }

  std::ostringstream msg ;
  msg << "unsupported hex plan code " << code ;
  throw std::runtime_error(msg.str()) ;
}

/// Return the number of children created by one Prism split code.
int prism_trace_child_count(int code, int nfold) {
  switch(code) {
  case 0:
    return 0 ;
  case 1:
    return 2 ;
  case 2:
    return nfold ;
  case 3:
    return 2*nfold ;
  default:
    break ;
  }

  std::ostringstream msg ;
  msg << "unsupported prism plan code " << code ;
  throw std::runtime_error(msg.str()) ;
}

/// Return the number of children created by one general Cell or DiamondCell code.
int general_trace_child_count(int code,
                              int local_fold,
                              int depth,
                              const std::vector<int>& root_child_folds) {
  switch(code) {
  case 0:
    return 0 ;
  case 1:
    return depth == 0 ? int(root_child_folds.size()) : 2*local_fold + 2 ;
  default:
    break ;
  }

  std::ostringstream msg ;
  msg << "unsupported general-cell plan code " << code ;
  throw std::runtime_error(msg.str()) ;
}

/// Dispatch child-count rules by cell family.
int trace_child_count(TraceKind kind,
                      int code,
                      int local_fold,
                      int depth,
                      const std::vector<int>& root_child_folds) {
  switch(kind) {
  case TRACE_HEX:
    return hex_trace_child_count(code) ;
  case TRACE_PRISM:
    return prism_trace_child_count(code, local_fold) ;
  case TRACE_GENERAL:
    return general_trace_child_count(code, local_fold, depth, root_child_folds) ;
  }
  throw std::runtime_error("unsupported trace kind") ;
}

/// Return the local fold value for each replayed child.
std::vector<int> trace_child_folds(TraceKind kind,
                                   int code,
                                   int local_fold,
                                   int depth,
                                   const std::vector<int>& root_child_folds) {
  const int count =
    trace_child_count(kind, code, local_fold, depth, root_child_folds) ;
  std::vector<int> folds(count, -1) ;

  if(kind == TRACE_HEX || code == 0) {
    return folds ;
  }

  if(kind == TRACE_PRISM) {
    if(code == 1) {
      for(int i = 0; i < count; ++i) {
        folds[i] = local_fold ;
      }
    } else {
      for(int i = 0; i < count; ++i) {
        folds[i] = 4 ;
      }
    }
    return folds ;
  }

  if(depth == 0) {
    return root_child_folds ;
  }

  if(count >= 1) {
    folds[0] = local_fold ;
  }
  if(count >= 2) {
    folds[1] = local_fold ;
  }
  for(int i = 2; i < count; ++i) {
    folds[i] = 3 ;
  }
  return folds ;
}

/// Name a split code in the replay trace output.
std::string code_meaning(TraceKind kind, int code, int depth) {
  switch(kind) {
  case TRACE_HEX:
    switch(code) {
    case 0: return "leaf" ;
    case 1: return "split_z" ;
    case 2: return "split_y" ;
    case 3: return "split_yz" ;
    case 4: return "split_x" ;
    case 5: return "split_xz" ;
    case 6: return "split_xy" ;
    case 7: return "split_xyz" ;
    default: return "unknown" ;
    }
  case TRACE_PRISM:
    switch(code) {
    case 0: return "leaf" ;
    case 1: return "split_z" ;
    case 2: return "split_xy_ring" ;
    case 3: return "split_xy_and_z" ;
    default: return "unknown" ;
    }
  case TRACE_GENERAL:
    switch(code) {
    case 0: return "leaf" ;
    case 1: return depth == 0 ? "split_cell_to_diamonds"
                              : "split_diamond_isotropic" ;
    default: return "unknown" ;
    }
  }
  return "unknown" ;
}

/// Replay a breadth-first split plan without calling the FVMAdapt cell classes.
PlanReplayTrace make_plan_replay_trace(TraceKind kind,
                                       const std::vector<char>& plan,
                                       int root_fold,
                                       const std::vector<int>& root_child_folds) {
  std::vector<TraceState> queue ;
  queue.push_back({0, -1, -1, 0, root_fold, "R"}) ;

  PlanReplayTrace trace ;
  trace.consumed_plan_entries = 0 ;
  trace.leaf_count = 0 ;

  for(size_t cursor = 0; cursor < queue.size(); ++cursor) {
    const TraceState state = queue[cursor] ;
    const bool explicit_code = trace.consumed_plan_entries < int(plan.size()) ;
    const int plan_index = explicit_code ? trace.consumed_plan_entries : -1 ;
    const int code = explicit_code ? int(plan[trace.consumed_plan_entries]) : 0 ;
    if(explicit_code) {
      ++trace.consumed_plan_entries ;
    }

    const int child_count = trace_child_count(kind,
                                              code,
                                              state.local_fold,
                                              state.depth,
                                              root_child_folds) ;
    TraceRow row ;
    row.id = state.id ;
    row.parent_id = state.parent_id ;
    row.child_slot = state.child_slot ;
    row.depth = state.depth ;
    row.local_fold = state.local_fold ;
    row.path = state.path ;
    row.plan_index = plan_index ;
    row.explicit_code = explicit_code ;
    row.code = code ;
    row.child_count = child_count ;
    row.first_child_id = child_count == 0 ? -1 : int(queue.size()) ;
    row.leaf_id = child_count == 0 ? ++trace.leaf_count : 0 ;
    trace.rows.push_back(row) ;

    const std::vector<int> child_folds =
      trace_child_folds(kind, code, state.local_fold, state.depth, root_child_folds) ;
    for(int child = 0; child < child_count; ++child) {
      std::ostringstream path ;
      path << state.path << "." << child ;
      TraceState child_state ;
      child_state.id = int(queue.size()) ;
      child_state.parent_id = state.id ;
      child_state.child_slot = child ;
      child_state.depth = state.depth + 1 ;
      child_state.local_fold = child_folds[child] ;
      child_state.path = path.str() ;
      queue.push_back(child_state) ;
    }
  }

  return trace ;
}

/// Map a CaseResult shape name to the replay model.
TraceKind trace_kind_for_shape(const std::string& shape) {
  if(shape == "hex") {
    return TRACE_HEX ;
  }
  if(shape == "prism") {
    return TRACE_PRISM ;
  }
  return TRACE_GENERAL ;
}

/// Return the root fold used by prism and general-cell replay.
int root_fold_for_shape(const std::string& shape) {
  if(shape == "prism") {
    return 3 ;
  }
  if(shape == "general_tet") {
    return 4 ;
  }
  if(shape == "general_pyramid") {
    return 5 ;
  }
  return -1 ;
}

/// Return the DiamondCell folds created by a general-cell root split.
std::vector<int> root_child_folds_for_shape(const std::string& shape) {
  if(shape == "general_tet") {
    return {3, 3, 3, 3} ;
  }
  if(shape == "general_pyramid") {
    return {3, 3, 3, 3, 4} ;
  }
  return std::vector<int>() ;
}

/// Build a replay trace and check its leaf count against the plan result.
PlanReplayTrace checked_plan_replay_trace(const CaseResult& result) {
  const TraceKind kind = trace_kind_for_shape(result.shape) ;
  const PlanReplayTrace trace =
    make_plan_replay_trace(kind,
                           result.input_plan,
                           root_fold_for_shape(result.shape),
                           root_child_folds_for_shape(result.shape)) ;

  if(trace.leaf_count != result.counted_leaves) {
    std::ostringstream msg ;
    msg << "plan replay trace for " << result.shape << " " << result.case_name
        << " has " << trace.leaf_count << " leaves, counted "
        << result.counted_leaves ;
    throw std::runtime_error(msg.str()) ;
  }

  return trace ;
}

/// Write the replay trace rows for one plan case.
void write_one_plan_replay(std::ofstream& out, const CaseResult& result) {
  const TraceKind kind = trace_kind_for_shape(result.shape) ;
  const PlanReplayTrace trace = checked_plan_replay_trace(result) ;

  out << "case " << result.shape << " " << result.case_name << "\n" ;
  out << "  input_cell_plan " << plan_string(result.input_plan) << "\n" ;
  out << "  canonical_cell_plan " << plan_string(result.canonical_plan) << "\n" ;
  out << "  consumed_plan_entries " << trace.consumed_plan_entries
      << " of " << result.input_plan.size() << "\n" ;
  if(trace.consumed_plan_entries < int(result.input_plan.size())) {
    out << "  unused_plan_tail "
        << plan_tail_string(result.input_plan, trace.consumed_plan_entries)
        << "\n" ;
  }
  out << "  replay_rows\n" ;
  out << "    # id parent child depth path local_fold plan_index source code "
      << "child_count first_child leaf meaning\n" ;
  for(size_t i = 0; i < trace.rows.size(); ++i) {
    const TraceRow& row = trace.rows[i] ;
    out << "    node "
        << row.id << " "
        << row.parent_id << " "
        << row.child_slot << " "
        << row.depth << " "
        << row.path << " "
        << row.local_fold << " "
        << row.plan_index << " "
        << (row.explicit_code ? "explicit" : "default") << " "
        << row.code << " "
        << row.child_count << " "
        << row.first_child_id << " "
        << row.leaf_id << " "
        << code_meaning(kind, row.code, row.depth) << "\n" ;
  }
}

} // namespace

/// Write breadth-first replay traces for all plan cases.
void write_plan_replay(const std::string& output_dir,
                       const std::vector<CaseResult>& results) {
  const std::string path = join_path(output_dir, "core_plan_replay.dat") ;
  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open plan replay output file") ;
  }

  out << "# FVMAdapt core plan replay trace\n" ;
  out << "# Each row is the breadth-first visit of one cell-tree node.\n" ;
  out << "# source=default means the plan vector was exhausted and replay used code 0.\n" ;
  out << "# local_fold is prism nfold, general Cell node count at depth 0, or\n" ;
  out << "# general DiamondCell nfold below depth 0. Hex rows use -1.\n" ;
  for(size_t i = 0; i < results.size(); ++i) {
    write_one_plan_replay(out, results[i]) ;
    out << "\n" ;
  }
}

/// Check replay-trace leaf counts for all plan cases.
void check_plan_replay_counts(const std::vector<CaseResult>& results) {
  for(size_t i = 0; i < results.size(); ++i) {
    checked_plan_replay_trace(results[i]) ;
  }
}
