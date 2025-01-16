# Changelog

<!--
## vX.Y.0 / 2025-MM-DD (Unreleased)

#### Breaking Changes

#### New Features

#### Enhancements

#### Bug Fixes

#### Misc.

#### MUST (Unmerged)

#### WIP (Unmerged)
-->


## v0.4.0 / 2025-01-16

#### Breaking Changes
 * [\#36](https://github.com/MolSSI/QCManyBody/pull/36) Feature -- as the embedded point charges aren't fully validated
   and (in the QCEngine computer function of the high-level interface) only work with Psi4 anyways, they are now hidden.
   Set environment variable `QCMANYBODY_EMBEDDING_CHARGES=1` to access functionality. @loriab

#### New Features

#### Enhancements
 * [\#36](https://github.com/MolSSI/QCManyBody/pull/36) Utils -- when adding up results from many molecular species,
   now the compensated sums `math.fsum` or `numpy.sum` are used.

#### Bug Fixes

#### Misc.
 * Maint -- pinned to QCElemental <0.70 to use only QCSchema v1.


## v0.3.0 / 2024-07-21

#### Breaking Changes

 * [\#28](https://github.com/MolSSI/QCManyBody/pull/28) Intf -- low-level "core" interface renamed from
   `ManyBodyCalculator` to `ManyBodyCore`. The old name will continue to work for a few months. Also, its file changed
   from `manybody.py` to `core.py` but it was already a top-level import. @loriab
 * [\#30](https://github.com/MolSSI/QCManyBody/pull/30) Intf -- low-level "core" interface now requires named arguments
   beyond the first recognizable ones (mol, bsse_type, levels). @loriab
 * [\#32](https://github.com/MolSSI/QCManyBody/pull/32) Intf -- "high-level" interface now no longer stores QCVariables
   (or any other results dicts) in extras @loriab
 * [\#32](https://github.com/MolSSI/QCManyBody/pull/32) Utils -- `qcmanybody.utils.collect_vars` now returns with keys
    from ManyBodyResultProperties rather than QCVariables. @loriab
 * [\#32](https://github.com/MolSSI/QCManyBody/pull/32) Utils -- arguments rearranged in
   `qcmanybody.tests.utils.run_qcengine` (use serial backend for core interface) to align with `ManyBodyCore` init
   arguments. @loriab

#### New Features

 * [\#32](https://github.com/MolSSI/QCManyBody/pull/32) Schema -- a new function
   `ManyBodyResultProperties.to_qcvariables()` returns a translation map to QCVariables keys. @loriab
 * [\#32](https://github.com/MolSSI/QCManyBody/pull/32) Schema -- a new function
   `qcmanybody.utils.translate_qcvariables(map)` switches between QCVariable and QCSchema keys. @loriab
 * [\#33](https://github.com/MolSSI/QCManyBody/pull/33) Schema -- `ManyBodySpecification.extras` added. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Schema -- Add schema_version to `AtomicSpecification`,
   `ManyBodySpecification`, `ManyBodyKeywords`, `ManyBodyInput`, and `ManyBodyResultProperties`. @loriab

#### Enhancements

 * [\#28](https://github.com/MolSSI/QCManyBody/pull/28) Intf -- high-level interface is now importable from the top level
   module. @loriab
 * [\#29](https://github.com/MolSSI/QCManyBody/pull/29) Maint -- QCEngine is needed only for the continuous running
   function of the high-level interface, so making it an optional dependency. @loriab
 * [\#30](https://github.com/MolSSI/QCManyBody/pull/30) Intf -- low-level "core" interface now accepts a molecule in
   partial schema dictionary format rather than requiring a constructed `qcelemental.Molecule` object. If the molecule
   is a single large fragment, an error is thrown. @loriab
 * [\#30](https://github.com/MolSSI/QCManyBody/pull/30) Docs -- add end-to-end demos in test_examples. @loriab
 * [\#31](https://github.com/MolSSI/QCManyBody/pull/31) Schema -- add "none" as a bsse_type alias to "nocp". @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Schema -- Allow environment variable QCMANYBODY_MAX_NBODY to
   influence the body-level to which `ManyBodyResultProperties` is defined added. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Schema -- added discriminator to input for
   `GeneralizedOptimizationInput` and `GeneralizedOptimizationResult` models to allow input from dicts (rather than
   models) in OptKing. Further specialized QCElemental. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Schema -- `ManyBodyResultProperties` is still only explicitly
   enumerated up to tetramers, but now it allows through higher-body fields when they match a pattern. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Maint -- start testing optimizations through QCEngine. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Util -- add `labeler(..., opaque=False)` option to produce
   eye-friendly `(1)@(1, 2)` style lablels as well as the semi-opaque internal style. Also always convent single ints
   to tuples now. Function `delabeler` can decode the new style. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Intf -- sort "core" `nbodies_per_mc_level` dictionary so model
   chemistries are in a predictable 1b, 2b, ..., supersystem order. Check that high-level (different data structure) agrees.
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Util -- add short ordinal model chemistry level (e.g., §A) to the
   `format_calc_plan` and `print_nbody_energy` summaries. @loriab
 * [\#34](https://github.com/MolSSI/QCManyBody/pull/34) Util -- add function `modelchem_labels` to associate n-body
   level, model chemistry level, one-char ordinal modelchem label, and n-bodies-covered modelchem label. @loriab


## v0.2.1 / 2024-05-14

#### Enhancements

* [\#27](https://github.com/MolSSI/QCManyBody/pull/27) Intf -- move high-level interface from
  `qcmb.qcng_computer.ManyBodyComputerQCNG` to `qcmb.computer.ManyBodyComputer`. Suppressed remaining printing. @loriab


## v0.2.0 / 2024-05-13

#### New Features

* [\#25](https://github.com/MolSSI/QCManyBody/pull/25) Schema -- added a new field `ManyBodyResults.component_results`
  to store subsystem `AtomicResult`s (or other Result if layered computation). By default this is not stored, but it can
  be with `ManyBodyInput.protocols.component_results = "all"`. @loriab
* [\#25](https://github.com/MolSSI/QCManyBody/pull/25) Schema -- added `GeneralizedOptimizationInput` and
  `GeneralizedOptimizationResult` models as temporary extensions of `OptimizationInput/Result` when the optimizer can
  run `ManyBodyInput`s through QCEngine, not just `AtomicInput`s. Needs special QCElemental, QCEngine, and OptKing to
  work for now. @loriab

#### Enhancements

* [\#22](https://github.com/MolSSI/QCManyBody/pull/22) Intf -- move high-level interface to main directory and remove
  unused functions. Most common route into interface is now: `qcmb.ManyBodyComputerQCNG.from_manybodyinput`. @loriab
* [\#25](https://github.com/MolSSI/QCManyBody/pull/25) Schema -- `AtomicSpecification.protocols` set in
  `ManyBodySpecification.specification` will now be observed. @loriab
* [\#26](https://github.com/MolSSI/QCManyBody/pull/26) Schema -- add `AtomicSpecification.extras`. @loriab

#### Bug Fixes

* [\#25](https://github.com/MolSSI/QCManyBody/pull/25) Schema -- `ManyBodyKeywords.embedding_charges` now default to None
  rather than empty dict. @loriab

#### Misc.

* [\#21](https://github.com/MolSSI/QCManyBody/pull/21) Docs -- installation and molecule. @bennybp
* [\#24](https://github.com/MolSSI/QCManyBody/pull/24) Intf -- release "frozen" pydantic on high-level
  ManyBodyComputerQCNG, so core-interface ManyBodyCalculator can live on the class. @loriab
* [\#26](https://github.com/MolSSI/QCManyBody/pull/26) Cleanup -- remove most debug printing. @loriab


## v0.1.0 / 2024-04-24

#### New Features

* Runs CP, NoCP, VMFC energies, gradients, and Hessians. @bennybp @loriab
