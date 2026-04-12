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
#ifndef SCHED_DB_H
#define SCHED_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <fact_db.h>
#include <rule.h>
#include <execute.h>

namespace Loci {

  /// @brief Describes one variable transfer between this processor and one peer.
  ///
  /// Each record names the variable being communicated and the peer involved.
  /// `send_set` holds the local entities to pack and send to `processor`.
  /// `recv_set` holds the local entities to fill from `processor`, in the order
  /// required by unpacking.
  ///
  /// These records are produced during scheduler analysis, cached by `sched_db`,
  /// and later assembled into communication blocks for barrier, reduction, loop,
  /// and recursion schedules.
  struct comm_info {
    variable v ;

    /// @brief Peer processor for this transfer record.
    /// For `send_set`, this is the destination rank.
    /// For `recv_set`, this is the source rank.
    int processor ;

    /// @brief Local entities to pack and send to `processor`.
    entitySet send_set ;

    /// @brief Local entities to fill from `processor`, in unpack order.
    sequence recv_set ;
  } ;


  /// @brief Send-side variable transfer grouped under a peer processor.
  ///
  /// The surrounding communication map or vector stores the peer rank; this
  /// record only names the variable and the local entities to pack for that peer.
  struct send_var_info {
    variable v ;

    /// @brief Local entities to pack for the associated peer.
    entitySet set ;
    send_var_info(variable iv, const entitySet &iset) : v(iv),set(iset) {}
  } ;


  /// @brief Receive-side variable transfer grouped under a peer processor.
  ///
  /// The surrounding communication map or vector stores the peer rank; this
  /// record only names the variable and the ordered local entities to fill from
  /// that peer.
  struct recv_var_info {
    variable v ;

    /// @brief Local entities to unpack into, in message order.
    sequence seq ;
    recv_var_info(variable iv, const sequence &iseq) : v(iv),seq(iseq) {}
  } ;


  /// @brief Compiler-side database of derived scheduling state.
  ///
  /// `sched_db` is initialized from a `fact_db` and then updated by scheduler
  /// passes as they analyze rule existence, propagate variable requests, and
  /// prepare distributed communication.
  ///
  /// It records per-variable scheduling metadata such as known existence,
  /// requested entities, rule-by-rule existential contributions, alias and
  /// synonym relationships, cached map-query results, communication records for
  /// distributed scheduling phases, and duplication-analysis data.
  ///
  /// Other scheduler components use this database to answer questions about what
  /// can be computed, what must be requested or communicated, and how a rule or
  /// variable should be scheduled.
  ///
  /// `sched_db` owns scheduling metadata only; the underlying fact stores and
  /// typed variables remain owned by `fact_db`.
  class sched_db {
  private:

    /// @brief Shared metadata for variables that point at the same
    ///   scheduling-data slot.
    ///
    /// Alias variables share one `sched_data` record so they can reuse map
    /// classification, cached map queries, and loop-rotation bookkeeping.
    struct sched_data {
      /// @brief Variables that share this `sched_data` slot.
      variableSet aliases ;

      /// @brief Variables linked by loop-rotation bookkeeping.
      variableSet rotations ;

      /// @brief Original-side names recorded when aliases are introduced into
      ///   this slot.
      variableSet antialiases ;

      /// @brief True when the backing store representation is a `MAP` or
      ///   `GPUMAP`.
      bool ismap ;

      /// @brief Cached `MapRep` used by `image()` and `preimage()` when
      ///   `ismap` is true.
      MapRepP minfo ;

      /// @brief Memoized `MapRep::image()` results keyed by the queried
      ///   domain.
      std::map<entitySet,entitySet> imageMap ;

      /// @brief Memoized `MapRep::preimage()` results keyed by the queried
      ///   codomain.
      /// The stored pair preserves both preimage variants returned by the map
      /// representation.
      std::map<entitySet,std::pair<entitySet,entitySet> > preimageMap ;

      sched_data() {}
      sched_data(variable v, storeRepP &st) {
        aliases += v ;
        ismap = (st->RepType() == Loci::MAP || st->RepType() == Loci::GPUMAP) ;
        if(ismap) minfo = MapRepP(st->getRep()) ;
      }
    } ;


    /// @brief One rule's contribution to the known existence of a variable.
    struct existential_info {
      /// @brief Variable instance associated with this rule contribution.
      /// Priority decoration on this name is used when overlapping producers
      /// are compared.
      variable v ;

      /// @brief Entities that remain attributed to the rule after any
      ///   priority-based conflict handling.
      entitySet exists ;

      existential_info() {}
      existential_info(const variable &vin,const entitySet &e) :
        v(vin), exists(e) {}
    } ;

    /// @brief Canonical per-variable scheduler state.
    ///
    /// Each variable handled directly by `vmap` owns one of these records.
    /// Shared alias/map metadata lives in `sched_infov`; this struct stores the
    /// existential, request, reduction, and duplication state attached to the
    /// variable itself.
    struct sched_info {
      /// @brief Index into `sched_infov` for shared alias, rotation, and map
      ///   metadata.
      int sched_info_ref ;

      // Legacy field kept commented out because it was initialized in `init()`
      // but not otherwise used by the current scheduler implementation.
      // entitySet fact_installed ;

      /// @brief Names that resolve to this `sched_info` entry through synonym
      ///   lookup.
      variableSet synonyms ;

      /// @brief Rule-keyed existential contributions merged into `existence`.
      std::map<rule,existential_info> exist_map ;

      /// @brief Union of directly known and rule-derived entities for the
      ///   variable.
      entitySet existence ;

      /// @brief Union of entities requested by downstream scheduling passes.
      entitySet requested ;

      /// @brief Distributed shadow-region entities for output-mapped apply or
      ///   reduction scheduling.
      /// Later passes use this set to exchange partial results with owner
      /// processors.
      entitySet shadow ;

      /// @brief Additional UNIT-rule requests created when mapped apply logic
      ///   needs target entities beyond the normal request set.
      entitySet extra_unit_request;

      /// @brief Maps each rule to the entities it can compute under the broader
      ///   duplication-analysis context.
      std::map<rule, entitySet> proc_able_map;

      /// @brief Maps each rule to the entities it can compute using only
      ///   processor-local context.
      std::map<rule, entitySet> my_proc_able_map;

      /// @brief Bitmask of enabled `duplicate_policy` selectors.
      unsigned int policy;

      /// @brief True when duplication has been selected for this variable.
      bool duplicate_variable;

      /// @brief Reduction-specific computable set used when deciding whether
      ///   owner communication can be avoided.
      entitySet reduce_proc_able_entities;

      /// @brief True when any producer rule for the variable maps its output
      ///   during reduction handling.
      bool reduction_outputmap;

      /// @brief Accumulated model-based computation-time estimates for the
      ///   original and duplicated schedules.
      double original_computation_time, duplication_computation_time;

      /// @brief Accumulated model-based communication-time estimates for the
      ///   original and duplicated schedules.
      double original_communication_time, duplication_communication_time;

      sched_info(int ref = -1) {
        sched_info_ref = ref ;
        policy = 0;
        duplicate_variable = false;
        reduce_proc_able_entities = ~EMPTY;
        reduction_outputmap = false;

        original_computation_time = 0;
        duplication_computation_time = 0;
        original_communication_time = 0;
        duplication_communication_time = 0;
      }
    } ;

    /// @brief Stores the two coefficients used by the scheduler's timing
    ///   models for duplication decisions.
    ///
    /// `ts` is the fixed-cost term and `tw` is the variable-cost term.
    /// Communication costs are estimated as
    /// `ts * processor_count + tw * packed_size`, while computation costs are
    /// estimated as `ts + tw * context_size`.
    ///
    /// The units come from the external model data loaded by the scheduler, so
    /// these values are treated as generic timing coefficients rather than
    /// hard-coded seconds.
    struct model {
      /// @brief Fixed-cost term in the timing model.
      double ts ;

      /// @brief Variable-cost term in the timing model.
      double tw ;

      /// @brief Returns true when `val` is a populated coefficient.
      /// @param[in] val Coefficient value to validate.
      /// @return True when `val` is not the invalid sentinel.
      static bool is_valid_val(double val) { return (val > -99999) ; }

      /// @brief Copies the stored coefficients into output references.
      /// @param[out] t0 Receives the fixed-cost term.
      /// @param[out] tc Receives the variable-cost term.
      void get_parameters(double &t0, double &tc) { t0 = ts; tc = tw ; }

      /// @brief Replaces the stored coefficients.
      /// @param[in] t0 Fixed-cost term to store.
      /// @param[in] tc Variable-cost term to store.
      void set_parameters(double t0, double tc) { ts = t0; tw = tc ; }

      model(double t0, double tc) {
        ts = t0 ;
        tw = tc ;
      }

      model() { ts = -100000; tw = -100000 ; }
    } ;

    /// @brief Variables already queued for later duplicate-work analysis.
    variableSet possible_duplicate_vars;

    /// @brief Per-rule computation timing models used by model-based
    ///   duplication decisions.
    std::map<rule, model> comp_model;

    /// @brief Communication timing model used by model-based duplication
    ///   decisions.
    model comm_model;

    /// @brief Registers a variable in the scheduler's canonical lookup tables.
    void register_variable(variable v) ;

    /// @brief Every variable name currently known to the scheduler, including
    ///   aliases and synonyms added after initialization.
    variableSet all_vars ;

    /// @brief Maps non-canonical synonym names to the canonical variables
    ///   stored in `vmap`.
    std::map<variable,variable> synonyms ;

    /// @brief Lookup-table type for variables that carry independent scheduler
    ///   state.
    typedef std::map<variable, sched_info> vmap_type ;

    /// @brief Per-variable scheduler state addressed after synonym resolution.
    vmap_type vmap ;

    /// @brief Shared `sched_data` table referenced by `sched_info_ref`.
    /// Alias variables point into the same entry in this table.
    std::vector<sched_data> sched_infov ;

    // Legacy field kept commented out because it is neither initialized nor
    // queried in the current scheduler implementation.
    //  intervalSet free_set ;

    /// @brief Sticky flag raised when scheduler analysis reports a conflict or
    ///   other invalid state.
    bool detected_errors ;

    /// @brief Phase-specific communication and execution caches populated
    ///   during request analysis.
    /// Send-entity caches hold mapped-output send domains. Communication-list
    /// caches hold variable-keyed `comm_info` records for barrier, reduction,
    /// loop-advance, and recursion phases. `exec_seq_map` stores rule-keyed
    /// execution domains.

    /// @brief Barrier-phase send sets used to build mapped-output
    ///   precommunication.
    std::map<variable,entitySet> barrier_send_entities_map ;

    /// @brief Recursive pre-phase send sets used to build mapped-output
    ///   precommunication.
    std::map<variable,entitySet> recurse_pre_send_entities_map ;

    /// @brief Barrier-phase clone/fill communication records derived from
    ///   variable requests.
    std::map<variable, std::list<comm_info> > barrier_clist_map;

    /// @brief Barrier-phase precommunication records built for mapped outputs.
    std::map<variable, std::list<comm_info> > barrier_plist_map;

    /// @brief Reduction-phase communication records that return shadow values
    ///   for owner-side joining.
    std::map<variable, std::list<comm_info> > reduce_rlist_map;

    /// @brief Reduction-phase clone/fill communication records derived from
    ///   variable requests.
    std::map<variable, std::list<comm_info> > reduce_clist_map;

    /// @brief Loop-advance communication records cached before time-level names
    ///   are restored.
    std::map<variable, std::list<comm_info> > loop_advance_list_map;

    /// @brief General recursion communication cache keyed by variable.
    std::map<variable, std::list<comm_info> > recurse_clist_map;

    /// @brief Recursive pre-phase request communication records.
    std::map<variable, std::list<comm_info> > recurse_pre_clist_map;

    /// @brief Recursive post-phase communication records issued after recursive
    ///   iterations.
    std::map<variable, std::list<comm_info> > recurse_post_clist_map;

    /// @brief Recursive pre-phase precommunication records for mapped outputs.
    std::map<variable, std::list<comm_info> > recurse_pre_plist_map;

    /// @brief Rule-keyed execution domains cached during request processing.
    /// The first write wins; later writes only warn.
    std::map<rule, entitySet> exec_seq_map;

  public:
    /// @brief Resolves a scheduler synonym to its canonical variable name.
    ///
    /// Synonyms share one canonical `sched_info` entry. Alias variables are not
    /// resolved by this helper; aliases keep separate names and share
    /// `sched_data`.
    /// @param[in] v Variable name to canonicalize.
    /// @return The canonical variable for `v`, or `v` if no synonym is
    ///   registered.
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      if((mi=synonyms.find(v)) != synonyms.end())
        return mi->second ;
      return v ;
    }

    /// @brief Selects a duplication-policy bit for a variable.
    ///
    /// `add_policy()` and `is_policy()` translate these selectors into bits in
    /// `sched_info::policy`: `NEVER` -> 1, `ALWAYS` -> 2, `MODEL_BASED` -> 4.
    enum duplicate_policy{NEVER, ALWAYS, MODEL_BASED};

    /// @brief Selects one of the cached `comm_info` lists by scheduler phase.
    ///
    /// The chosen value determines which communication cache
    /// `get_comm_info_list()` and `update_comm_info_list()` read or update for
    /// barrier, reduction, loop-advance, or recursive scheduling.
    enum list_type{BARRIER_CLIST, BARRIER_PLIST,
                   REDUCE_RLIST, REDUCE_CLIST,
                   LOOP_ADVANCE_LIST,RECURSE_CLIST,
                   RECURSE_PRE_CLIST, RECURSE_POST_CLIST, RECURSE_PRE_PLIST};

    /// @brief Selects which cached mapped-output send-entity set to query or
    ///   update.
    ///
    /// These caches are populated during existence analysis and later used
    /// while building precommunication for barrier or recursive pre phases.
    enum send_entities_type{BARRIER, RECURSE_PRE};

    /// @brief Creates an empty scheduler database with no registered
    ///   variables.
    ///
    /// Call `init()` before using variable-dependent scheduler state.
    sched_db() ;

    /// @brief Destroys the scheduler metadata owned by this database.
    ~sched_db() ;

    /// @brief Creates a scheduler database initialized from the typed facts in
    ///   `facts`.
    ///
    /// This convenience constructor clears the error flag and then seeds the
    /// scheduler tables from the variables currently typed in `facts`.
    /// @param[in] facts Fact database whose typed variables seed the initial
    ///   scheduler state.
    sched_db(fact_db &facts) ;

    /// @brief Seeds scheduler state from the typed variables in `facts`.
    ///
    /// For each typed variable, this records its current store domain as known
    /// existence, creates its canonical `sched_info`, initializes its shared
    /// `sched_data`, and records the variable in `all_vars`.
    ///
    /// This routine is intended for a fresh `sched_db` or one that has been
    /// explicitly prepared for reinitialization.
    /// @param[in] facts Fact database whose typed variables seed scheduler
    ///   state.
    void init(fact_db &facts) ;

    /// @brief Returns whether scheduler analysis has reported a sticky error.
    ///
    /// Top-level schedule builders check this flag after schedule construction
    /// and abort if any MPI rank reports a failure.
    bool errors_found() {return detected_errors ;}

    /// @brief Clears the sticky scheduler error flag.
    ///
    /// This resets only the recorded flag; it does not repair the underlying
    /// state that originally triggered the error.
    void clear_errors() {detected_errors = false ;}

    /// @brief Marks scheduler analysis as failed.
    ///
    /// Use this when a scheduler pass detects an invalid or unsupported state
    /// but continues far enough for the caller to report diagnostics.
    void set_error() { detected_errors = true ; }

    /// @brief Adds shared scheduling data for a new variable and registers `v`.
    ///
    /// The supplied `sched_data` is appended to `sched_infov`, and `v` then
    /// receives a `sched_info` record pointing at that new shared-data slot.
    /// @param[in] v Variable to register.
    /// @param[in] data Shared scheduling data to append for `v`.
    void install_sched_data(variable v, sched_data data) ;

    /// @brief Installs scheduler state for `v` in the variable registry.
    ///
    /// Use this when `v` should reuse an existing `sched_data` slot described
    /// by `info.sched_info_ref`, such as when registering an alias.
    /// @param[in] v Variable to register.
    /// @param[in] info Scheduler state to install for `v`.
    /// @note Aborts if `v` is already registered.
    void install_sched_info(variable v, sched_info info) ;

    /// @brief Returns the synonym set recorded for `v`'s canonical scheduler
    ///   entry.
    ///
    /// This covers names linked through `synonym_variable()`. It does not
    /// include alias or rotation relationships.
    /// @param[in] v Variable whose synonym set is requested.
    /// @note `v` must already be known to the scheduler; use
    ///   `try_get_synonyms()` for unknown-safe lookup.
    /// @return The synonym set for `v`'s canonical scheduler entry.
    variableSet get_synonyms(variable v) const {
      return get_sched_info(v).synonyms ;
    }

    /// @brief Returns the synonym set recorded for `v`'s canonical scheduler
    ///   entry.
    /// @param[in] v Variable whose synonym set is requested.
    /// @return `EMPTY` when `v` is not known to the scheduler.
    variableSet try_get_synonyms(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end())
        return variableSet(EMPTY) ;
      else
        return (mi->second).synonyms ;
    }

    /// @brief Returns the alias group that shares `sched_data` with `v`.
    ///
    /// Alias relationships come from `alias_variable()` and are distinct from
    /// synonym resolution.
    /// @param[in] v Variable whose alias group is requested.
    /// @note `v` must already be known to the scheduler; use
    ///   `try_get_aliases()` for unknown-safe lookup.
    /// @return The alias group that shares `sched_data` with `v`.
    variableSet get_aliases(variable v) const {
      return get_sched_data(v).aliases ;
    }

    /// @brief Returns the alias group that shares `sched_data` with `v`.
    /// @param[in] v Variable whose alias group is requested.
    /// @return `EMPTY` when `v` is not known to the scheduler.
    variableSet try_get_aliases(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end())
        return variableSet(EMPTY) ;
      else
        return sched_infov[(mi->second).sched_info_ref].aliases ;
    }

    /// @brief Returns the anti-alias names recorded in `v`'s shared alias
    ///   metadata.
    ///
    /// These names are added when aliases are introduced and are stored in the
    /// shared `sched_data` slot rather than in per-variable synonym state.
    /// @param[in] v Variable whose anti-alias names are requested.
    /// @note `v` must already be known to the scheduler; use
    ///   `try_get_antialiases()` for unknown-safe lookup.
    /// @return The anti-alias names recorded for `v`.
    variableSet get_antialiases(variable v) const {
      return get_sched_data(v).antialiases ;
    }

    /// @brief Returns the anti-alias names recorded in `v`'s shared alias
    ///   metadata.
    /// @param[in] v Variable whose anti-alias names are requested.
    /// @return `EMPTY` when `v` is not known to the scheduler.
    variableSet try_get_antialiases(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end())
        return variableSet(EMPTY) ;
      else
        return sched_infov[(mi->second).sched_info_ref].antialiases ;
    }

    /// @brief Returns the loop-rotation group recorded for `v`.
    ///
    /// Rotation groups are registered by `set_variable_rotations()` so later
    /// scheduling passes can treat those variables together.
    /// @param[in] v Variable whose rotation group is requested.
    /// @note `v` must already be known to the scheduler; use
    ///   `try_get_rotations()` for unknown-safe lookup.
    /// @return The loop-rotation group recorded for `v`.
    variableSet get_rotations(variable v) const {
      return get_sched_data(v).rotations ;
    }

    /// @brief Returns the loop-rotation group recorded for `v`.
    /// @param[in] v Variable whose rotation group is requested.
    /// @return `EMPTY` when `v` is not known to the scheduler.
    variableSet try_get_rotations(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end())
        return variableSet(EMPTY) ;
      else
        return sched_infov[(mi->second).sched_info_ref].rotations ;
    }

    /// @brief Ensures that `v` is typed in both the scheduler and `facts`.
    ///
    /// This creates an intensional fact for `v` in `facts`. If the scheduler
    /// is not already tracking `v`, it also installs a fresh `sched_data`
    /// record seeded from `st`. Existing scheduler state for `v` is left in
    /// place.
    /// @param[in] v Variable being typed.
    /// @param[in] st Store representation to associate with `v`.
    /// @param[in,out] facts Fact database updated with the corresponding
    ///   intensional fact.
    void set_variable_type(variable v, storeRepP st, fact_db &facts) ;
    void set_variable_type(std::string vname,const storeRepP st,
                           fact_db &facts) {
      set_variable_type(variable(vname),st, facts) ;
    }

    void set_variable_type(variable v, store_instance &si, fact_db &facts) {
      set_variable_type(v,si.Rep(), facts) ;
    }
    void set_variable_type(std::string vname, store_instance &si,
                           fact_db &facts) {
      set_variable_type(variable(vname),si, facts) ;
    }

    /// @brief Alternate overload set for callers that do not pass a
    ///   `fact_db`.
    /// @param[in] v Variable being typed.
    /// @param[in] st Store representation to associate with `v`.
    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname,const storeRepP st) {
      set_variable_type(variable(vname),st) ;
    }

    void set_variable_type(variable v, store_instance &si) {
      set_variable_type(v,si.Rep()) ;
    }
    void set_variable_type(std::string vname, store_instance &si) {
      set_variable_type(variable(vname),si) ;
    }

    // Legacy `variable_is_fact_at()` overloads remain commented out because the
    // scheduler no longer uses them.
    //  void variable_is_fact_at(variable v,entitySet s, fact_db &facts) ;
    //     void variable_is_fact_at(std::string vname, const entitySet s, fact_db &facts)
    //       { variable_is_fact_at(variable(vname),s, facts) ; }

    //     void variable_is_fact_at(variable v,entitySet s) ;
    //     void variable_is_fact_at(std::string vname, const entitySet s)
    //       { variable_is_fact_at(variable(vname),s) ; }

    /// @brief Stores `rot` as the complete rotation group for every variable in
    ///   the set.
    ///
    /// Each member must already be known to the scheduler. The group is written
    /// into shared scheduling data so later loop-scheduling passes can recover
    /// the same rotation set from any member. Aliases that share the same
    /// `sched_data` record observe the same rotation metadata.
    /// @param[in] rot Complete rotation group to record for each member.
    void set_variable_rotations(variableSet rot) ;

    /// @brief Registers `alias` as a new variable that shares `v`'s
    ///   `sched_data`, and mirrors the relationship into `facts`.
    ///
    /// Aliases share type- and map-level scheduling metadata through the same
    /// `sched_data` slot, but each alias keeps its own `sched_info` fields such
    /// as existence, requests, and rule ownership.
    /// @param[in] v Existing scheduler variable that provides the shared
    ///   `sched_data`.
    /// @param[in] alias New alias name to register.
    /// @param[in,out] facts Fact database updated to mirror the alias
    ///   relationship.
    /// @note If exactly one of `v` and `alias` already exists, the existing
    ///   variable is used as the source regardless of argument order.
    void alias_variable(variable v, variable alias, fact_db &facts) ;
    void alias_variable(std::string vname, std::string alias, fact_db &facts) {
      alias_variable(variable(vname),variable(alias), facts) ;
    }

    /// @brief Registers `alias` as a new variable that shares `v`'s
    ///   `sched_data`.
    ///
    /// This is the scheduler-only form of `alias_variable(..., fact_db&)`.
    /// @param[in] v Existing scheduler variable that provides the shared
    ///   `sched_data`.
    /// @param[in] alias New alias name to register.
    void alias_variable(variable v, variable alias) ;
    void alias_variable(std::string vname, std::string alias) {
      alias_variable(variable(vname),variable(alias)) ;
    }

    /// @brief Registers `synonym` as an alternate name for canonical variable
    ///   `v`, and mirrors the relationship into `facts`.
    ///
    /// After registration, lookups for `synonym` resolve through
    /// `remove_synonym()` to `v`, so the synonym shares the same `sched_info`
    /// state instead of getting a separate scheduler record.
    /// @param[in] v Canonical variable name.
    /// @param[in] synonym Alternate name to resolve through `v`.
    /// @param[in,out] facts Fact database updated to mirror the synonym
    ///   relationship.
    void synonym_variable(variable v, variable synonym, fact_db &facts) ;
    void synonym_variable(std::string vname, std::string synonym,
                          fact_db &facts) {
      synonym_variable(variable(vname),variable(synonym), facts) ;
    }

    /// @brief Registers `synonym` as an alternate name for canonical variable
    ///   `v`.
    ///
    /// This is the scheduler-only form of `synonym_variable(..., fact_db&)`.
    /// @param[in] v Canonical variable name.
    /// @param[in] synonym Alternate name to resolve through `v`.
    void synonym_variable(variable v, variable synonym) ;
    void synonym_variable(std::string vname, std::string synonym) {
      synonym_variable(variable(vname),variable(synonym)) ;
    }

    /// @brief Records or extends the subset of `v` produced by rule `f`.
    /// @param[in] v Variable instance, including any priority decoration.
    /// @param[in] f Rule that produces the entities in `x`.
    /// @param[in] x Produced entities to merge into `f`'s existential record.
    ///
    /// Repeated calls for the same rule accumulate into one entry. If the new
    /// entities overlap entities already owned by a different rule, variable
    /// priority annotations are used to resolve ownership when possible.
    /// Unresolved overlaps are reported as conflicts and mark the scheduler as
    /// having detected errors. Calls with `EMPTY` still register `f` in the
    /// existential map.
    void set_existential_info(variable v,rule f,entitySet x) ;

    /// @brief Returns the rules currently present in `v`'s existential map.
    ///
    /// This reports registered rule entries, even if a rule's current entity
    /// contribution is empty.
    /// @param[in] v Variable whose existential producers are requested.
    /// @return The rules currently present in `v`'s existential map.
    ruleSet get_existential_rules(variable v) ;

    /// @brief Returns `v`'s canonical `sched_info` after synonym resolution.
    /// @param[in] v Variable whose scheduler state is requested.
    /// @note Alias variables are not collapsed here; aliases keep separate
    ///   `sched_info` records.
    /// @note Aborts if the canonical variable is not known to the scheduler.
    /// @return The canonical read-only scheduler record for `v`.
    const sched_info & get_sched_info(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end()) {
        cerr << "const get_sched_info: variable " << v << " does not exist" << endl ;
        abort() ;
      }
      return mi->second ;
    }

    /// @brief Returns the shared `sched_data` record referenced by `v`'s
    ///   canonical scheduler entry.
    ///
    /// Variables in the same alias group share this record.
    /// @param[in] v Variable whose shared scheduling data is requested.
    /// @return The shared read-only scheduling data for `v`.
    const sched_data & get_sched_data(variable v) const {
      return sched_infov[get_sched_info(v).sched_info_ref] ;
    }

    /// @brief Returns `v`'s canonical `sched_info` after synonym resolution.
    /// @param[in] v Variable whose scheduler state is requested.
    /// @note Alias variables are not collapsed here; aliases keep separate
    ///   `sched_info` records.
    /// @note Unknown variables are routed through the scheduler's `EMPTY`
    ///   placeholder record instead of aborting.
    /// @return The canonical mutable scheduler record for `v`.
    sched_info & get_sched_info(variable v);

    /// @brief Returns the shared `sched_data` record referenced by `v`'s
    ///   canonical scheduler entry.
    ///
    /// Variables in the same alias group share this record, so mutations affect
    /// the whole alias group.
    /// @param[in] v Variable whose shared scheduling data is requested.
    /// @return The shared mutable scheduling data for `v`.
    sched_data & get_sched_data(variable v) {
      return sched_infov[get_sched_info(v).sched_info_ref] ;
    }

    /// @brief Returns true when `v`'s shared scheduling data is backed by a map
    ///   representation.
    ///
    /// Alias variables report the same result because they share `sched_data`.
    /// @param[in] v Variable to inspect.
    /// @return True when `v` is backed by a map representation.
    bool is_a_Map(variable v) {
      return get_sched_data(v).ismap ;
    }

    /// @brief Returns the entities currently attributed to rule `f` for `v`.
    ///
    /// This queries `v`'s per-variable existential map, so synonyms resolve to
    /// the same record while aliases keep separate rule-attribution state.
    /// @param[in] v Variable whose existential map is queried.
    /// @param[in] f Rule whose attributed entities are requested.
    /// @return The stored existential subset for `f`, or `EMPTY` when `f` has
    ///   no entry for `v`.
    entitySet get_existential_info(variable v, rule f) {
      sched_info &finfo = get_sched_info(v) ;
      std::map<rule,existential_info>::const_iterator mi ;
      mi = finfo.exist_map.find(f) ;
      if(mi!=finfo.exist_map.end()) {
        return mi->second.exists ;
      } else
        return EMPTY ;
    }

    /// @brief Returns the accumulated existence currently known for `v`.
    ///
    /// This includes direct existence updates and entities recorded through
    /// `set_existential_info()`. Because existence lives in `sched_info`,
    /// synonyms share it while aliases keep separate existence sets.
    /// @param[in] v Variable whose known existence is requested.
    /// @return The accumulated existence set, or `~EMPTY` for time variables.
    entitySet variable_existence(variable v) ;

    /// @brief Adds `e` to the pending request set for `v`.
    ///
    /// Requests accumulate rather than replace the current set. Time variables
    /// ignore request updates because they are treated as existing everywhere.
    /// @param[in] v Variable whose request set is updated.
    /// @param[in] e Requested entities to add.
    void variable_request(variable v, entitySet e) ;

    /// @brief Clears pending requests and extra UNIT-rule requests for every
    ///   tracked `sched_info` entry.
    ///
    /// Synonyms reset together because they resolve to the same scheduler
    /// entry, while aliases are cleared independently.
    void clear_variable_request() {
      for(vmap_type::iterator mi = vmap.begin();
          mi != vmap.end(); mi++){
        mi->second.requested=EMPTY;
        mi->second.extra_unit_request = EMPTY;
      }
    }

    // void reset_variable_request(variable v, entitySet e=EMPTY) ;

    /// @brief Replaces the distributed shadow set tracked for `v`.
    ///
    /// This communication-oriented set is used by later distributed scheduling
    /// passes, especially output-mapped and reduction cases. Because shadow
    /// lives in `sched_info`, synonyms share it while aliases keep separate
    /// shadow sets.
    /// @param[in] v Variable whose shadow set is replaced.
    /// @param[in] e New shadow set to store.
    void set_variable_shadow(variable v, entitySet e) {
      get_sched_info(v).shadow = e ;
    }

    /// @brief Adds entities to the distributed shadow set tracked for `v`.
    /// @param[in] v Variable whose shadow set is updated.
    /// @param[in] e Shadow entities to add.
    void variable_shadow(variable v, entitySet e) {
      get_sched_info(v).shadow += e ;
    }

    /// @brief Returns the distributed shadow set currently tracked for `v`.
    ///
    /// Shadow state lives in `sched_info`, so synonyms share it while aliases
    /// keep separate shadow sets.
    /// @param[in] v Variable whose shadow set is requested.
    /// @return The distributed shadow set currently recorded for `v`.
    entitySet get_variable_shadow(variable v) const {
      return get_sched_info(v).shadow ;
    }

    /// @brief Returns the requested portion of `v` that overlaps rule `f`'s
    ///   recorded existential contribution.
    /// @param[in] f Rule whose requested overlap is requested.
    /// @param[in] v Variable whose request set is queried.
    /// @return The intersection of `f`'s existential contribution and the
    ///   current request set for `v`, or `EMPTY` when `f` has no existential
    ///   entry for `v`.
    entitySet get_variable_request(rule f, variable v) ;

    /// @brief Returns the full pending request set currently stored for `v`.
    ///
    /// Requests live in `sched_info`, so synonyms share them while aliases keep
    /// separate request sets.
    /// @param[in] v Variable whose request set is requested.
    /// @return The full pending request set for `v`.
    entitySet get_variable_requests(variable v) const {
      return  get_sched_info(v).requested ;
    }

    /// @brief Returns every variable name currently registered with the
    ///   scheduler.
    ///
    /// This is the scheduler's `all_vars` registry, so it includes canonical
    /// variables plus any aliases or synonyms that were added later.
    /// @return Every variable name currently registered with the scheduler.
    variableSet get_typed_variables() const { return all_vars ; }

    /// @brief Returns the memoized image of `e` through map variable `v`.
    ///
    /// The first lookup delegates to the stored `MapRep`; later lookups for
    /// the same entity set reuse the cached result. Because this cache lives in
    /// `sched_data`, aliases share it, and synonyms resolve to the same cache
    /// through their canonical scheduler entry.
    /// @param[in] v Map variable to query.
    /// @param[in] e Domain entities whose image is requested.
    /// @return `EMPTY` when `v` is not a map variable.
    entitySet image(variable v, entitySet e) ;

    /// @brief Returns the memoized preimage information of `e` through map
    ///   variable `v`.
    ///
    /// The stored pair is whatever the underlying `MapRep::preimage()`
    /// produces for `e`: for simple maps the two sets match, while multi-valued
    /// maps use the pair to distinguish fully-contained from merely touching
    /// source entities. Because this cache lives in `sched_data`, aliases share
    /// it, and synonyms resolve to the same cache through their canonical
    /// scheduler entry.
    /// @param[in] v Map variable to query.
    /// @param[in] e Codomain entities whose preimage is requested.
    /// @return The cached `MapRep::preimage()` result, or `(EMPTY, EMPTY)` when
    ///   `v` is not a map variable.
    std::pair<entitySet,entitySet> preimage(variable v, entitySet e) ;

    /// @brief Replaces the raw duplication-policy bitmask stored for `v`.
    ///
    /// This is the low-level setter behind `add_policy()`. Because policy
    /// lives in `sched_info`, synonyms share it while aliases keep separate
    /// policy state. The mask itself is only an input to later duplication
    /// selection; it does not set `duplicate_variable`.
    /// @param[in] v Variable whose raw policy mask is replaced.
    /// @param[in] p New raw policy mask.
    void set_policy(variable v, unsigned int p) { get_sched_info(v).policy = p;}

    /// @brief Returns the raw duplication-policy bitmask stored for `v`.
    ///
    /// Use `duplicate_policy`, `add_policy()`, and `is_policy()` when you want
    /// the symbolic interpretation of the bits.
    /// @param[in] v Variable whose policy mask is requested.
    /// @return The raw duplication-policy bitmask stored for `v`.
    unsigned int get_policy(variable v) { return get_sched_info(v).policy; }

    /// @brief Adds duplication policy `p` to `v`'s current policy mask.
    ///
    /// This accumulates policy bits rather than replacing the existing mask.
    /// The update is applied across `v`'s synonym set, but not across aliases.
    /// Later duplication analysis interprets that mask and records the final
    /// choice with `set_duplicate_variable()`.
    /// @param[in] v Variable whose policy mask is updated.
    /// @param[in] p Policy bit to add.
    void add_policy(variable v, duplicate_policy p) ;

    /// @brief Tests whether duplication policy `p` is active for `v`.
    ///
    /// This checks the symbolic `duplicate_policy` bit in `v`'s raw policy
    /// mask. It does not mean duplication has already been selected.
    /// @param[in] v Variable whose policy mask is queried.
    /// @param[in] p Policy bit to test.
    /// @return True when policy `p` is active for `v`.
    bool is_policy(variable v, duplicate_policy p) ;

    /// @brief Returns true when duplication has been selected for `v`.
    ///
    /// This is the scheduler's stored duplication decision, not just a policy
    /// recommendation. Later request-minimization and communication code uses
    /// this flag to switch from owner-based requests to duplicated local
    /// computation.
    /// @param[in] v Variable whose duplication decision is queried.
    /// @return True when duplication has been selected for `v`.
    bool is_duplicate_variable(variable v) { return get_sched_info(v).duplicate_variable; }

    /// @brief Sets the duplicate-variable flag for `v`'s synonym set.
    ///
    /// This records the scheduler's final duplication choice for those names.
    /// Aliases keep their own duplicate-variable flag because they have
    /// separate `sched_info` entries.
    /// @param[in] v Variable whose synonym set is updated.
    /// @param[in] p New duplicate-variable flag.
    void set_duplicate_variable(variable v, bool p) {
      variableSet syns = get_synonyms(v) ;
      for(variableSet::const_iterator vi = syns.begin(); vi != syns.end(); vi++)
        get_sched_info(*vi).duplicate_variable = p ;
    }

    /// @brief Returns the broader duplicate-work compute set recorded for rule
    ///   `f` and variable `v`.
    ///
    /// This per-rule set is used when `v` has been marked as a duplicate
    /// variable. Unlike `get_my_proc_able_entities()`, it is recorded before
    /// the `my_entities` restriction is applied, so it may include non-owned
    /// entities this processor can recompute under work duplication. Because it
    /// lives in `sched_info`, synonyms share it while aliases keep separate
    /// entries.
    /// @param[in] v Variable whose duplicate-work compute set is requested.
    /// @param[in] f Rule whose compute set is requested.
    /// @return The duplication-expanded compute set recorded for (`v`, `f`).
    entitySet get_proc_able_entities(variable v, rule f) {
      sched_info &finfo = get_sched_info(v) ;
      std::map<rule, entitySet>::const_iterator mi ;
      mi = finfo.proc_able_map.find(f) ;
      if(mi != finfo.proc_able_map.end())
        return mi->second ;
      else
        return EMPTY ;
    }

    /// @brief Adds entities to the duplicate-work compute set recorded for rule
    ///   `f` and variable `v`.
    ///
    /// Repeated calls accumulate into the cached set rather than replacing it.
    /// @param[in] v Variable whose duplicate-work compute set is updated.
    /// @param[in] f Rule whose compute set is updated.
    /// @param[in] x Entities to add.
    void set_proc_able_entities(variable v, rule f, entitySet x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.proc_able_map[f] += x ;
    }

    /// @brief Returns the processor-local compute set recorded for rule `f`
    ///   and variable `v`.
    ///
    /// This narrower per-rule set is recorded after restricting the rule's
    /// context to this processor's local entities. Later request processing
    /// uses it when `v` is not being duplicated.
    /// @param[in] v Variable whose local compute set is requested.
    /// @param[in] f Rule whose local compute set is requested.
    /// @return The `my_entities`-restricted compute set for (`v`, `f`).
    entitySet get_my_proc_able_entities(variable v, rule f) {
      sched_info &finfo = get_sched_info(v) ;
      std::map<rule, entitySet>::const_iterator mi ;
      mi = finfo.my_proc_able_map.find(f) ;
      if(mi != finfo.my_proc_able_map.end())
        return mi->second ;
      else
        return EMPTY ;
    }

    /// @brief Adds entities to the processor-local compute set recorded for
    ///   rule `f` and variable `v`.
    ///
    /// Repeated calls accumulate into the cached set rather than replacing it.
    /// @param[in] v Variable whose local compute set is updated.
    /// @param[in] f Rule whose local compute set is updated.
    /// @param[in] x Entities to add.
    void set_my_proc_able_entities(variable v, rule f, entitySet x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.my_proc_able_map[f] += x ;
    }

    /// @brief Returns the reduction-specific computable summary recorded for
    ///   `v`.
    ///
    /// This is a later summary used by duplicated reduction handling, not a
    /// per-rule map like `get_proc_able_entities()`.
    /// @param[in] v Reduction variable whose summary is requested.
    /// @return The reduction-specific computable summary for `v`.
    entitySet get_reduce_proc_able_entities(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.reduce_proc_able_entities) ;
    }

    /// @brief Replaces the reduction-specific computable summary recorded for
    ///   `v`.
    /// @param[in] v Reduction variable whose summary is replaced.
    /// @param[in] x New reduction-specific computable summary.
    void set_reduce_proc_able_entities(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.reduce_proc_able_entities = x ;
    }

    /// @brief Returns true when duplicated reduction handling for `v` must use
    ///   the output-mapped case.
    ///
    /// This flag is set during reduction analysis when at least one producer of
    /// `v` maps its output. In that case later duplicated-reduction logic uses
    /// the conservative `reduce_filter` path instead of widening by
    /// `get_reduce_proc_able_entities(v)`.
    /// @param[in] v Reduction variable to query.
    /// @return True when `v` must use the output-mapped reduction path.
    bool is_reduction_outputmap(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.reduction_outputmap) ;
    }

    /// @brief Records whether duplicated reduction handling for `v` should use
    ///   the output-mapped case.
    ///
    /// This value is computed by the reduction-analysis pass.
    /// @param[in] v Reduction variable to update.
    /// @param[in] x New reduction-output-mapping flag.
    void set_reduction_outputmap(variable v, bool x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.reduction_outputmap = x ;
    }

    /// @brief Adds extra UNIT-rule requests for `v`.
    ///
    /// This accumulates request entities that mapped apply/reduction handling
    /// needs UNIT rules to cover in addition to the normal request set. These
    /// extras are folded in only for UNIT-rule request processing and are
    /// cleared by `clear_variable_request()`.
    /// @param[in] v Variable whose extra UNIT requests are updated.
    /// @param[in] x Extra UNIT-rule requests to add.
    void add_extra_unit_request(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.extra_unit_request += x ;
    }
    // void reset_extra_unit_request(variable v, entitySet x=EMPTY) {
//       sched_info &finfo = get_sched_info(v) ;
//       finfo.extra_unit_request = x ;
//     }

    /// @brief Returns the extra UNIT-rule requests currently recorded for `v`.
    ///
    /// These requests are only added into UNIT-rule processing, not every rule
    /// request for `v`.
    /// @param[in] v Variable whose extra UNIT requests are requested.
    /// @return The extra UNIT-rule requests currently recorded for `v`.
    entitySet get_extra_unit_request(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.extra_unit_request) ;
    }

    /// @brief Adds known existence for `v` without recording a producing rule.
    ///
    /// Despite the name, this extends the current existence set rather than
    /// replacing it, and it does not create an `exist_map` producer entry. Use
    /// it when existence is known directly, especially for map-related cases
    /// that do not have a scheduler rule to attribute.
    /// @param[in] v Variable whose known existence is extended.
    /// @param[in] x Directly known existence to add.
    void set_variable_existence(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.existence += x ;
    }

    /// @brief Loads the timing models used by `MODEL_BASED` duplication.
    ///
    /// `comm_ts` and `comm_tw` populate the global communication model used to
    /// estimate communication as startup plus size-based cost. `comp_info`
    /// populates the per-rule computation models used to estimate rule work as
    /// a fixed term plus a context-size term.
    /// @param[in] comm_ts Fixed communication startup coefficient.
    /// @param[in] comm_tw Variable communication coefficient.
    /// @param[in] comp_info Per-rule computation-model coefficients.
    void add_model_info(double comm_ts, double comm_tw, const std::map<rule,
                        std::pair<double, double> > &comp_info) ;

    // These helpers accumulate model estimates computed while comparing the
    // original schedule against the duplicated-work alternative for `v`.
    // "original" means the owner-communication schedule before duplication is
    // selected, while "duplication" means the alternative that recomputes
    // locally available entities and communicates only the remaining requests.
    /// @brief Adds to `v`'s original-schedule computation estimate.
    /// @param[in] v Variable whose estimate is updated.
    /// @param[in] add Computation time to add.
    void add_original_computation_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.original_computation_time += add ;
    }

    /// @brief Adds to `v`'s duplicated-schedule computation estimate.
    /// @param[in] v Variable whose estimate is updated.
    /// @param[in] add Computation time to add.
    void add_duplication_computation_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.duplication_computation_time += add ;
    }

    /// @brief Adds to `v`'s original-schedule communication estimate.
    /// @param[in] v Variable whose estimate is updated.
    /// @param[in] add Communication time to add.
    void add_original_communication_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.original_communication_time += add ;
    }

    /// @brief Adds to `v`'s duplicated-schedule communication estimate.
    /// @param[in] v Variable whose estimate is updated.
    /// @param[in] add Communication time to add.
    void add_duplication_communication_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v) ;
      finfo.duplication_communication_time += add ;
    }

    /// @brief Returns `v`'s accumulated original-schedule computation estimate.
    /// @param[in] v Variable whose estimate is requested.
    /// @return The accumulated original-schedule computation estimate for `v`.
    double get_precalculated_original_computation_time(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.original_computation_time) ;
    }

    /// @brief Returns `v`'s accumulated duplicated-schedule computation
    ///   estimate.
    /// @param[in] v Variable whose estimate is requested.
    /// @return The accumulated duplicated-schedule computation estimate for `v`.
    double get_precalculated_duplication_computation_time(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.duplication_computation_time) ;
    }

    /// @brief Returns `v`'s accumulated original-schedule communication
    ///   estimate.
    /// @param[in] v Variable whose estimate is requested.
    /// @return The accumulated original-schedule communication estimate for `v`.
    double get_precalculated_original_communication_time(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.original_communication_time) ;
    }

    /// @brief Returns `v`'s accumulated duplicated-schedule communication
    ///   estimate.
    /// @param[in] v Variable whose estimate is requested.
    /// @return The accumulated duplicated-schedule communication estimate for `v`.
    double get_precalculated_duplication_communication_time(variable v) {
      sched_info &finfo = get_sched_info(v) ;
      return(finfo.duplication_communication_time) ;
    }

    /// @brief Returns the computation timing model for rule `r`, or an invalid
    ///   sentinel model.
    ///
    /// Missing rules return `model(-100000, -100000)`. Check the coefficients
    /// with `model::is_valid_val()` before using them in timing math.
    /// @param[in] r Rule whose computation model is requested.
    /// @return The computation timing model for `r`, or an invalid sentinel.
    model get_comp_model(rule r) const {
      std::map<rule, model>::const_iterator mi = comp_model.find(r) ;
      if(mi != comp_model.end())
        return mi->second ;
      else {
        return model(-100000, -100000) ;
      }
    }

    /// @brief Returns the global communication timing model installed by
    ///   `add_model_info()`.
    ///
    /// Before model data is loaded, the coefficients remain invalid sentinel
    /// values.
    /// @return The global communication timing model.
    model get_comm_model() const { return comm_model ; }

    /// @brief Returns the duplicate-analysis worklist for future target sets.
    ///
    /// This is not the final duplication decision; use
    /// `is_duplicate_variable()` for that. The worklist is carried forward for
    /// variables whose duplication should be reconsidered when they appear in a
    /// later target set.
    /// @return The current duplicate-analysis worklist.
    variableSet get_possible_duplicate_vars() const { return possible_duplicate_vars ; }

    /// @brief Adds variables and their synonyms to the future duplicate-analysis
    ///   worklist.
    ///
    /// The set is append-only. Synonyms are inserted for each variable, while
    /// aliases are not expanded here.
    /// @param[in] vars Variables to add to the future duplicate-analysis
    ///   worklist.
    void add_possible_duplicate_vars(variableSet vars) {
      for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); vi++) {
        variableSet syns = get_synonyms(*vi) ;
        possible_duplicate_vars += syns ;
      }
    }

    /// @brief Writes a readable summary of the scheduler's current existential
    ///   state.
    ///
    /// The summary includes each tracked variable's container kind, synonym
    /// set, known existence, and current request set.
    /// @param[in] facts Fact database used to identify each variable's
    ///   container kind.
    /// @param[in,out] s Output stream that receives the summary text.
    /// @return `s`.
    std::ostream &print_summary(fact_db &facts, std::ostream &s) ;

    /// @brief Returns cached send-entity sets for variables in `eset` whose
    ///   existential producers map output entities.
    ///
    /// Only variables with at least one output-mapped producer are returned.
    /// Missing cache entries are skipped with a warning to `debugout`.
    /// @param[in] eset Variables whose send-entity caches are requested.
    /// @param[in] e Phase-specific send-entity cache to query.
    /// @return Cached send-entity sets for the applicable variables in `eset`.
    std::vector<std::pair<variable,entitySet> > get_send_entities(variableSet eset, send_entities_type e) ;

    /// @brief Stores send-entity cache entries for a specific scheduler phase.
    ///
    /// Existing entries are overwritten, with a warning when the same variable
    /// is written more than once for that phase.
    /// @param[in] evec Send-entity entries to store.
    /// @param[in] e Phase-specific send-entity cache to update.
    void  update_send_entities( const std::vector<std::pair<variable,entitySet> >& evec, send_entities_type e) ;

    /// @brief Returns one of the cached communication lists for variables in
    ///   `eset`.
    ///
    /// Retrieval filters the cache to the requested variables, then repacks the
    /// result so all sends appear first in processor order, followed by all
    /// receives in processor order.
    /// @param[in] eset Variables whose cached communication records are
    ///   requested.
    /// @param[in] facts Fact database associated with the schedule.
    /// @param[in] e Phase-specific communication cache to query.
    /// @return The filtered and reordered communication list.
    std::list<comm_info> get_comm_info_list(variableSet eset, fact_db& facts, list_type e) const ;

    /// @brief Appends communication records into one of the phase-specific
    ///   communication caches.
    ///
    /// Records are stored per variable and are not cleared or deduplicated by
    /// this helper.
    /// @param[in] elist Communication records to append.
    /// @param[in] e Phase-specific communication cache to update.
    void update_comm_info_list(const std::list<comm_info>& elist, list_type e) ;

    /// @brief Returns the cached execution sequence for `r`.
    /// @param[in] r Rule whose cached execution sequence is requested.
    /// @return The cached sequence, or `EMPTY` when no sequence has been
    ///   stored.
    ///
    /// Reading before the sequence is populated emits a warning to `debugout`.
    entitySet get_exec_seq(rule r) {
      entitySet re = EMPTY ;
      std::map<rule,entitySet>::const_iterator mi = exec_seq_map.find(r) ;
      if(mi != exec_seq_map.end()){
        re += mi->second ;
      }else{
        debugout << "WARNING: exec_seq_map is read before write at rule "
                 << r << endl ;
      }
      return re ;
    }

    /// @brief Stores the execution sequence for `r` the first time it is seen.
    ///
    /// Later writes for the same rule are ignored and reported as warnings to
    /// `debugout`.
    /// @param[in] r Rule whose execution sequence is being cached.
    /// @param[in] eset Execution sequence to cache for `r`.
    void update_exec_seq(rule r, entitySet eset) {
      std::map<rule,entitySet>::const_iterator mi = exec_seq_map.find(r) ;
      if(mi == exec_seq_map.end()) {
        exec_seq_map[r] = eset ;
      }else{
        debugout << "WARNING: exec_seq_map is written more than once at rule "
                 << r << endl ;
      }
    }
  } ;
}

#endif
