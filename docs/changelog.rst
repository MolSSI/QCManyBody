=========
Changelog
=========

v0.5.2 / 2026-02-16
===================

v0.5.1 / 2025-06-14
===================

Bug Fixes
---------

* :pr:`40` Core -- Fix bug where an input ``fix_symmetry='c1'`` molecule didn't apply no symmetry to its fragments.


v0.5.0 / 2025-06-13
===================

Breaking Changes
----------------

New Features
------------

* :pr:`38` Feature -- alpha ManyBody QCSchema v2 added accessible through ``from qcmanybody.models.v2 import ManyBodyInput`` etc. Changes are:

  * In v2, ``ManyBodyResult`` gained ``native_files`` and ``molecule`` fields. The latter is unchanged by QCManyBody from the input. Also, ``cluster_results`` became required so the protocol can process it (can always be empty dict).
  * In v2, ``ManyBodyInput`` gained ``provenance`` and ``id`` fields.
  * In v2, ``ManyBodySpecification`` gained a ``program`` field, empty by default since prog may be supplied by argument.
  * ``ManyBodyResult.component_results`` in v1 is now in v2 ``ManyBodyResult.cluster_results``.
  * ``ManyBodyProtocols.component_results`` in v1 is now in v2 ``ManyBodyProtocols.cluster_results`` with the same default.
  * In v2, ``ManyBodyInput.extras`` field was removed. Extras should be on ``ManyBodySpecification``.
  * In v2, ``ManyBodyProperties``, ``ManyBodyKeywords``, and ``ManyBodySpecification`` lost their ``schema_version`` field.
  * All v2 models got their ``schema_name`` standardized to ``qcschema_many_body_<type>``. ``ManyBodyProtocols`` got a ``schema_name`` for the first time.
  * ``ManyBodyResult.success`` now in v2 is fixed True.
  * ``ManyBodyResult.properties`` composed of ``ManyBodyResultProperties`` in v1 is now in v2 composed of ``ManyBodyProperties``.
  * ``ManyBodySpecification.protocols`` composed of ``AtomicResultProtocols`` in v1 is now in v2 composed of ``AtomicProtocols``.
  * ``ManyBodyResult.component_properties`` composed of ``AtomicResultProperties`` in v1 is now in v2 ``ManyBodyResult.cluster_properties`` composed of ``AtomicProperties``.

Enhancements
------------

* :pr:`38` Utils -- updated the precise math function for arrays according to the NumPy deprecation warning.
* :pr:`39` Docs -- Added high-level interface overview.

Bug Fixes
---------

* :pr:`38` Docs -- Fix typos in core docs page.

Misc.
-----

* :pr:`38` Maint -- Pydantic package version must be v2 or v1 >=1.10.17. This ensures the v2 API will be available for optional QCSchema v2, while QCManyBody remains on v1 API (importable from package v2).
* :pr:`38` Maint -- Updated license and package spec.


v0.4.0 / 2025-01-16
===================

Breaking Changes
----------------

* :pr:`36` Feature -- as the embedded point charges aren't fully validated and (in the QCEngine computer function of the high-level interface) only work with Psi4 anyways, they are now hidden. Set environment variable ``QCMANYBODY_EMBEDDING_CHARGES=1`` to access functionality.

New Features
------------

Enhancements
------------

* :pr:`36` Utils -- when adding up results from many molecular species, now the compensated sums ``math.fsum`` or ``numpy.sum`` are used.

Bug Fixes
---------

Misc.
-----

* Maint -- pinned to QCElemental <0.70 to use only QCSchema v1.


v0.3.0 / 2024-07-21
===================

Breaking Changes
----------------

* :pr:`28` Intf -- low-level "core" interface renamed from ``ManyBodyCalculator`` to ``ManyBodyCore``. The old name will continue to work for a few months. Also, its file changed from ``manybody.py`` to ``core.py`` but it was already a top-level import.
* :pr:`30` Intf -- low-level "core" interface now requires named arguments beyond the first recognizable ones (mol, bsse_type, levels).
* :pr:`32` Intf -- "high-level" interface now no longer stores QCVariables (or any other results dicts) in extras.
* :pr:`32` Utils -- ``qcmanybody.utils.collect_vars`` now returns with keys from ManyBodyResultProperties rather than QCVariables.
* :pr:`32` Utils -- arguments rearranged in ``qcmanybody.tests.utils.run_qcengine`` (use serial backend for core interface) to align with ``ManyBodyCore`` init arguments.

New Features
------------

* :pr:`32` Schema -- a new function ``ManyBodyResultProperties.to_qcvariables()`` returns a translation map to QCVariables keys.
* :pr:`32` Schema -- a new function ``qcmanybody.utils.translate_qcvariables(map)`` switches between QCVariable and QCSchema keys.
* :pr:`33` Schema -- ``ManyBodySpecification.extras`` added.
* :pr:`34` Schema -- Add schema_version to ``AtomicSpecification``, ``ManyBodySpecification``, ``ManyBodyKeywords``, ``ManyBodyInput``, and ``ManyBodyResultProperties``.

Enhancements
------------

* :pr:`28` Intf -- high-level interface is now importable from the top level module.
* :pr:`29` Maint -- QCEngine is needed only for the continuous running function of the high-level interface, so making it an optional dependency.
* :pr:`30` Intf -- low-level "core" interface now accepts a molecule in partial schema dictionary format rather than requiring a constructed ``qcelemental.Molecule`` object. If the molecule is a single large fragment, an error is thrown.
* :pr:`30` Docs -- add end-to-end demos in test_examples.
* :pr:`31` Schema -- add "none" as a bsse_type alias to "nocp".
* :pr:`34` Schema -- Allow environment variable QCMANYBODY_MAX_NBODY to influence the body-level to which ``ManyBodyResultProperties`` is defined added.
* :pr:`34` Schema -- added discriminator to input for ``GeneralizedOptimizationInput`` and ``GeneralizedOptimizationResult`` models to allow input from dicts (rather than models) in OptKing. Further specialized QCElemental.
* :pr:`34` Schema -- ``ManyBodyResultProperties`` is still only explicitly enumerated up to tetramers, but now it allows through higher-body fields when they match a pattern.
* :pr:`34` Maint -- start testing optimizations through QCEngine.
* :pr:`34` Util -- add ``labeler(..., opaque=False)`` option to produce eye-friendly ``(1)@(1, 2)`` style labels as well as the semi-opaque internal style. Also always convert single ints to tuples now. Function ``delabeler`` can decode the new style.
* :pr:`34` Intf -- sort "core" ``nbodies_per_mc_level`` dictionary so model chemistries are in a predictable 1b, 2b, ..., supersystem order. Check that high-level (different data structure) agrees.
* :pr:`34` Util -- add short ordinal model chemistry level (e.g., §A) to the ``format_calc_plan`` and ``print_nbody_energy`` summaries.
* :pr:`34` Util -- add function ``modelchem_labels`` to associate n-body level, model chemistry level, one-char ordinal modelchem label, and n-bodies-covered modelchem label.


v0.2.1 / 2024-05-14
===================

Enhancements
------------

* :pr:`27` Intf -- move high-level interface from ``qcmb.qcng_computer.ManyBodyComputerQCNG`` to ``qcmb.computer.ManyBodyComputer``. Suppressed remaining printing.


v0.2.0 / 2024-05-13
===================

New Features
------------

* :pr:`25` Schema -- added a new field ``ManyBodyResults.component_results`` to store subsystem ``AtomicResult`` s (or other Result if layered computation). By default this is not stored, but it can be with ``ManyBodyInput.protocols.component_results = "all"``.
* :pr:`25` Schema -- added ``GeneralizedOptimizationInput`` and ``GeneralizedOptimizationResult`` models as temporary extensions of ``OptimizationInput/Result`` when the optimizer can run ``ManyBodyInput`` s through QCEngine, not just ``AtomicInput`` s. Needs special QCElemental, QCEngine, and OptKing to work for now.

Enhancements
------------

* :pr:`22` Intf -- move high-level interface to main directory and remove unused functions. Most common route into interface is now: ``qcmb.ManyBodyComputerQCNG.from_manybodyinput``.
* :pr:`25` Schema -- ``AtomicSpecification.protocols`` set in ``ManyBodySpecification.specification`` will now be observed.
* :pr:`26` Schema -- add ``AtomicSpecification.extras``.

Bug Fixes
---------

* :pr:`25` Schema -- ``ManyBodyKeywords.embedding_charges`` now default to None rather than empty dict.

Misc.
-----

* :pr:`21` Docs -- installation and molecule.
* :pr:`24` Intf -- release "frozen" pydantic on high-level ManyBodyComputerQCNG, so core-interface ManyBodyCalculator can live on the class.
* :pr:`26` Cleanup -- remove most debug printing.


v0.1.0 / 2024-04-24
===================

New Features
------------

* Runs CP, NoCP, VMFC energies, gradients, and Hessians.
