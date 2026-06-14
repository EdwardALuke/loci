# FVMAdapt Refactor Resolution Ledger {#fvmadapt_refactor_ledger}

This ledger tracks cleanup items that should be resolved before the current
FVMAdapt documentation and `$rule` translation work is considered finished.
Keep entries concrete and remove them only after the code builds and the
relevant smoke tests have been run.

## Open Items

1. Keep reusable apply-rule reducers out of the FVMAdapt documentation PR.
   - Boolean reducers such as generic logical-and/logical-or operators are
     useful beyond FVMAdapt and should be proposed in a separate framework PR
     if they are still wanted.
   - This branch should keep translated FVMAdapt rules on existing local
     reducers unless that separate PR is already merged.
   - Known local reducers:
     `src/include/FVMAdapt/defines.h` defines `logicalAnd`, and
     `src/FVMAdapt/check_folded_face.loci` defines `logicalOr`.

2. Finish the temporary side-by-side rule translations.
   - `src/FVMAdapt/balance_general_cell.loci` currently preserves old C++ rule
     classes in block comments below the translated `$rule` forms.
   - Once the translated rules are reviewed and accepted, remove the commented
     old implementations so the file does not carry duplicate source truth.

3. Finish the `.lh` declaration split.
   - Standard mesh facts should come from `FVM.lh`; do not redeclare
     `lower`, `upper`, `boundary_map`, `face2node`, `face2edge`, `edge2node`,
     `pos`, `geom_cells`, `faces`, `interior_faces`, or `fileNumber(X)` in
     converted FVMAdapt rule files.
   - `face2edge`, `edge2node`, `faces`, `interior_faces`, `cells`, and
     `ghost_cells` have been promoted to `FVM.lh` because they are generated
     and consumed by shared FVM grid/setup infrastructure.
   - FVMAdapt-specific facts used only by converted internal rules should live
     beside those rules in `src/FVMAdapt/fvmadapt_internal.lh`, not under
     `src/include/FVMAdapt`.
   - `src/FVMAdapt/Makefile` adds the source directory to `INCLUDES` so
     configured `OBJ/FVMAdapt` builds can resolve this private include through
     the symlinked Makefile.
   - Before finishing, scan converted files for leftover local `$type`
     declarations that should be promoted or removed, including temporary
     FVMAdapt-only facts such as `balancedPointSet1`.

4. Run final verification from a clean configured tree.
   - At minimum, run `OBJ/lpp/lpp` on translated files and a targeted
     `make -C OBJ/FVMAdapt LOCI_BASE=..` build.
   - Run the existing FVMAdapt quick test path documented in
     `src/FVMAdapt/doc/testing.md` before claiming behavior is unchanged.

5. Check documentation consistency after code cleanup.
   - Keep `src/FVMAdapt/doc/concepts.md` aligned with the final translation
     patterns.
   - Remove ledger items that were resolved, and leave only genuinely open
     follow-up work.
