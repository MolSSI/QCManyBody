# Parallel Execution Remediation Roadmap

This plan tracks the follow-up work required to harden the `parallel_execution` branch before merging to `main`. Each task includes a problem summary, desired outcome, concrete deliverables, and dependencies.

## Summary Dashboard

| ID | Priority | Area | Status | Owner | Target | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | High | QCEngine integration | ☑ Done | TBD | 2025-10-15 | Spec mapping + driver propagation covered by new executor + tests |
| T2 | High | Multiprocessing | ☑ Done | TBD | 2025-10-18 | Spawn-safe worker refactor with regression test + docs |
| T3 | Medium | Metrics | ☑ Done | TBD | 2025-10-10 | Real serial baselines, per-level timings, histogram support |
| T4 | Medium | Dependency validation | ☐ Not started | TBD | 2025-10-22 | Extend level batching tests (supersystem/VMFC) |
| T5 | Low | Config UX | ☐ Not started | TBD | 2025-10-05 | Allow `ParallelConfig` without qcengine install when disabled |

Legend: ☐ Not started ▣ In progress ☑ Done

---

## T1 — Align QCEngine Inputs With Model Chemistry Abstractions

**Problem**  
`ParallelManyBodyExecutor.execute_fragment` currently lowers the model-chemistry label into `AtomicInput.model["method"]` and hard-codes `driver="energy"`. Real jobs use higher-level aliases (e.g., `"nwc-ccsd/tz"`) and may request gradients or Hessians. This mapping will fail QCEngine validation and silently ignore requested drivers/keywords.

**Impact**  
Parallel runs on production specs will error or return incorrect physics, blocking branch adoption.

**Objectives**
- Reuse the existing `ManyBodyInput` / `AtomicSpecification` machinery to translate each fragment into a valid `AtomicInput`.
- Preserve per-fragment driver selection (energy/gradient/hessian) and keywords.
- Ensure program selection follows the original compute map (e.g., psi4 vs nwchem) rather than a global config default.

**Deliverables**
1. Adapter that consumes the compute plan entries and returns `(AtomicInput, program)` pairs for the executor. *(Implemented in `ParallelManyBodyExecutor._build_fragment_tasks` with specification resolution.)*  
2. Unit tests for gradient/hessian fragments proving driver propagation. *(See `test_parallel_executor.py::test_fragment_driver_propagation` and related cases.)*  
3. End-to-end validation (serial vs parallel) on a mixed-level example, using real QCEngine when available. *(Covered by the water4 integration exercised in `test_multiprocessing_serialization.py`; enable QCEngine to run with real programs.)*

**Dependencies**
- ManyBodyCore compute map introspection.  
- Regression tester updates (T4) to provide mixed-model fixtures.

---

## T2 — Make Multiprocessing Execution Portable

**Problem**  
`ProcessPoolExecutor` pickles `ParallelManyBodyExecutor.execute_fragment`, which drags in complex state (qcengine, Molecule). This fails on spawn-based multiprocessing (Windows, macOS) and raises serialization errors—the current branch only works under fork.

**Impact**  
Multiprocessing mode is unreliable, undermining the key deliverable of the branch.

**Objectives**
- Move fragment execution into a top-level helper (or module-level function) that accepts a minimal payload (serialized `AtomicInput`, execution settings) to avoid pickling the executor object.
- Cache heavy objects (e.g., programs, options) via `initializer` and `initargs` for performance while keeping workers stateless across pickles.
- Add automated regression that runs a multiprocessing test under Python’s `spawn` start method.

**Deliverables**
1. Refactored executor wiring using a pickle-safe worker function. *(Completed via `_execute_fragment_worker` + initializer in `parallel.py`.)*  
2. New test in `test_multiprocessing_serialization.py` asserting no pickling errors under `multiprocessing.set_start_method("spawn", force=True)`. *(Implemented and runs against the real water4 example.)*  
3. Documentation update explaining worker start modes and environment expectations. *(Added guidance in `parallel-execution-project/docs/usage/quick-start.md`.)*

**Dependencies**
- Completion of T1 (payload now an `AtomicInput`).

---

## T3 — Report Meaningful Performance Metrics

**Problem**  
`execution_stats["speedup_factor"]` currently divides total parallel time by itself, always yielding 1.0. The dashboard cannot detect regressions or wins.

**Impact**  
Benchmarks and validation artifacts mislead stakeholders, obscuring performance regressions.

**Objectives**
- Measure a real serial baseline: either reuse the serial executor or accumulate per-fragment timing before batching.
- Record per-level timings to surface load imbalance in reports.
- Update benchmark scripts to visualize speedup and per-level critical path.

**Deliverables**
1. Enhanced statistics struct with serial baseline, per-level stats, optional histograms. *(Implemented in `ParallelManyBodyExecutor` with level timing capture, serial baseline synthesis, and configurable histograms.)*  
2. Updated `scripts/benchmark_parallel_performance.py` producing plots/tables using the new data. *(Script now consumes measured speedup/serial metrics and reports histogram-ready data.)*  
3. Unit tests asserting non-trivial speedup values when parallel execution is faster than serial. *(See `TestParallelManyBodyExecutor::test_speedup_measurement_threading` and associated histogram/disable checks.)*

**Dependencies**
- None, but T1 ensures serial reference uses the same core path.

---

## T4 — Strengthen Dependency-Level Validation

**Problem**  
`iterate_molecules_by_level()` assumes every level can execute independently. Supersystem-only tasks, VMFC bookkeeping, or mixed model-chemistry levels may break this assumption, but current tests don’t cover those cases.

**Impact**  
Potential silent correctness bugs once users enable advanced BSSE options or supersystem levels.

**Objectives**
- Add fixtures covering supersystem entries, VMFC-only workloads, and mixed-level (outer supersystem, inner multi-level) scenarios.
- Ensure dependency graph output matches original `iterate_molecules()` order and respects required sequencing (e.g., supersystem last, no missing fragments).
- Update validation scripts to run these fixtures in both serial and parallel modes.

**Deliverables**
1. New tests in `test_dependency_graph.py` and `test_integration_dependency_graph.py` covering the above scenarios.  
2. Regression data under `qcmanybody/testing` that exercises the new cases.  
3. Documentation note in `docs/validation/ultra-strict-validation.md` describing coverage.

**Dependencies**
- T1 (ensures fragment generation handles program/driver combos used in tests).

---

## T5 — Relax QCEngine Dependency Requirement When Disabled

**Problem**  
`ParallelConfig.__post_init__` raises if QCEngine isn’t importable, even when `use_qcengine=False`. Scripts that rely on placeholder execution (e.g., quick demos, validation) require qcengine to be installed unnecessarily.

**Impact**  
Onboards face a hard dependency on qcengine for docs/tests, complicating installation and CI configuration.

**Objectives**
- Defer the import guard until `use_qcengine=True` is requested.  
- Document optional dependencies clearly and add CI smoke tests for placeholder mode.

**Deliverables**
1. Guarded import logic in `ParallelConfig` with targeted unit tests ensuring `use_qcengine=False` works without the package.  
2. README/docs update clarifying optional installs.  
3. Lightweight CI job running placeholder-mode demos without qcengine present.

**Dependencies**
- None.

---

## Cross-Cutting Milestones

1. **Design Review (2025-10-07)** — Confirm solution approach for T1/T2 with maintainers; update this plan based on feedback.  
2. **Implementation Window (2025-10-08 → 2025-10-22)** — Execute T1–T5, landing changes behind feature flags where prudent.  
3. **Validation & Benchmarks (2025-10-24)** — Re-run `scripts/validate_parallel_execution.py` and performance suite with real qcengine workloads; publish updated results in `PARALLEL_SYSTEM_STATUS.md`.  
4. **Merge Readiness (2025-10-28)** — Final review, ensure CI green, resolve documentation deltas.

---

## Tracking & Updates

- Update the status checkbox and notes column as work progresses.  
- Record meeting outcomes or blockers directly beneath each task section.  
- On completion, move deliverables to the repository and link PR numbers in the associated section.
