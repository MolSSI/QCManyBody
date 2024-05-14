# Changelog

<!--
## vX.Y.0 / 2024-MM-DD (Unreleased)

#### Breaking Changes

#### New Features

#### Enhancements

#### Bug Fixes

#### Misc.

#### MUST (Unmerged)

#### WIP (Unmerged)
-->


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

