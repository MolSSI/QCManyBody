#!/usr/bin/env python3
"""Convert HMBE fragmentation inputs into QCManyBody-ready scripts."""

# =============================================================================
# Usage Guide
# -----------------------------------------------------------------------------
# Purpose
# -------
# This CLI tool translates an HMBE-style JSON (*.inp) fragmentation input into
# a standalone Python helper that constructs a ``ManyBodyInput`` model and a
# ``ParallelConfig`` for QCManyBody's parallel executor.
#
# Basic Invocation
# ----------------
#     python hmbe_to_qcmanybody.py water_000.inp
#
# Output Location
# ---------------
# By default the script writes ``<stem>_qcmanybody.py`` next to the source
# input. Provide ``--output-file`` to control the destination or overwrite an
# existing helper.
#
# Key Flags
# ---------
# --output-file PATH        : Override the generated script path.
# --execution-mode MODE     : Choose ``multiprocessing`` (default),
#                             ``threading``, or ``serial`` for the
#                             parallel executor.
# --max-workers N           : Cap the worker pool size for parallel execution.
# --memory-mb MB            : Per-worker memory budget passed to QCEngine.
# --timeout-seconds SEC     : Per-fragment timeout safeguard in seconds.
# --no-qcengine             : Emit a config that skips QCEngine dispatch.
# --qc-program NAME         : Force the QC backend (default follows
#                             input file).
# --basis-set NAME          : Override the basis set encoded in the input file.
# --bsse-type FLAG          : Repeatable; request additional BSSE treatments.
# --no-return-total-data    : Disable total energy aggregation in results.
# --supersystem-ie-only     : Restrict the expansion to the supersystem IE.
#
# Tips
# ----
# * Run with ``--help`` for the full option list and detailed descriptions.
# * Generated helper scripts are importable; they expose
#   ``create_manybody_input`` and ``create_parallel_config`` for reuse.
# * Use a temporary ``--output-file`` path when you only need a quick preview.
# =============================================================================

from __future__ import annotations

import argparse
import json
import sys
import textwrap
from pathlib import Path
from string import Template
from typing import Dict, Iterable, List, Optional, Tuple

import pprint


def _format_literal(value: object) -> str:
    """Return a literal string suitable for embedding in generated code."""

    literal = pprint.pformat(value, width=79, compact=False)
    if "\n" not in literal:
        return literal

    indented = textwrap.indent(literal, " " * 4).lstrip()
    return f"(\n    {indented}\n)"


def _collect_system_data(
    hmbe_data: Dict,
) -> Tuple[
    List[str],
    List[List[float]],
    List[List[int]],
    int,
    int,
]:
    """Flatten the fragment tree into QCManyBody-friendly lists."""

    chemical_system = hmbe_data["chemical_system"]
    fragments: List[List[int]] = []
    symbols: List[str] = []
    geometry: List[List[float]] = []

    atom_index = 0
    for parent_fragment in chemical_system.get("fragments", []):
        for subfragment in parent_fragment.get("fragments", []):
            sub_symbols: List[str] = subfragment["symbols"]
            raw_geometry: List[float] = subfragment["geometry"]

            fragment_atom_indices: List[int] = []
            for atom_offset, symbol in enumerate(sub_symbols):
                symbols.append(symbol)
                xyz = raw_geometry[atom_offset * 3: atom_offset * 3 + 3]
                geometry.append([float(coord) for coord in xyz])
                fragment_atom_indices.append(atom_index)
                atom_index += 1

            fragments.append(fragment_atom_indices)

    total_charge = chemical_system.get("molecular_charge")
    if total_charge is None:
        fragment_charges = (
            parent.get("molecular_charge", 0)
            for parent in chemical_system.get("fragments", [])
        )
        total_charge = sum(fragment_charges)

    total_multiplicity = chemical_system.get("molecular_multiplicity")
    if total_multiplicity is None:
        for parent in chemical_system.get("fragments", []):
            total_multiplicity = parent.get("molecular_multiplicity")
            if total_multiplicity is not None:
                break
        else:
            total_multiplicity = 1

    return (
        symbols,
        geometry,
        fragments,
        int(total_charge),
        int(total_multiplicity),
    )


def convert_hmbe_to_qcmanybody(
    hmbe_file: str,
    output_file: Optional[str] = None,
    *,
    bsse_types: Iterable[str],
    return_total_data: bool,
    supersystem_ie_only: bool,
    execution_mode: str,
    max_workers: int,
    use_qcengine: bool,
    qc_program_override: Optional[str],
    basis_override: Optional[str],
    memory_limit_mb: int,
    timeout_seconds: int,
) -> str:
    """Convert an HMBE JSON input file to a QCManyBody helper script."""

    with open(hmbe_file, "r", encoding="utf-8") as handle:
        hmbe_data = json.load(handle)

    runtime_params = hmbe_data["runtime_params"]

    (
        symbols,
        geometry,
        fragments,
        total_charge,
        total_multiplicity,
    ) = _collect_system_data(hmbe_data)

    mbe_orders: List[int] = runtime_params.get("mbe_orders", [])
    if mbe_orders:
        level_orders = sorted({int(level) for level in mbe_orders})
    else:
        level_orders = list(range(1, len(fragments) + 1))
    max_nbody = max(level_orders) if level_orders else len(fragments)

    driver = runtime_params.get("driver", "energy").lower()
    model_params = runtime_params.get("model", {})
    method = model_params.get("method", "hf").lower()
    basis = (basis_override or model_params.get("basis", "sto-3g")).lower()
    base_label = f"{method}_{basis}".replace(" ", "_").replace("-", "_")

    level_to_spec = {
        int(level): f"{base_label}_lvl{int(level)}"
        for level in level_orders
    }
    if not level_to_spec:
        level_to_spec[1] = f"{base_label}_lvl1"

    program_map = {
        "psi4": "psi4",
        "qchem": "qchem",
        "gaussian": "gaussian",
        "nwchem": "nwchem",
        "orca": "orca",
    }

    program_name = (runtime_params.get("program") or "psi4").lower()
    qc_program = qc_program_override or program_map.get(
        program_name,
        program_name,
    )

    spec_keywords = runtime_params.get("keywords", {})

    bsse_type_list = list(
        dict.fromkeys(str(bt).lower() for bt in (bsse_types or ["nocp"]))
    )

    if output_file is None:
        input_path = Path(hmbe_file)
        output_file = input_path.stem + "_qcmanybody.py"

    symbols_literal = _format_literal(symbols)
    geometry_literal = _format_literal(geometry)
    fragments_literal = _format_literal(fragments)
    level_map_literal = _format_literal(level_to_spec)
    spec_keywords_literal = _format_literal(spec_keywords)
    bsse_literal = _format_literal(bsse_type_list)

    script_template = Template(
        textwrap.dedent(
            '''#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from qcelemental.models import Molecule  # noqa: E402
from qcmanybody.models import (  # noqa: E402
    AtomicSpecification,
    ManyBodyInput,
    ManyBodyKeywords,
    ManyBodySpecification,
)
from qcmanybody.parallel import (  # noqa: E402
    ParallelConfig,
    ParallelManyBodyExecutor,
)

SYMBOLS = $symbols_literal
GEOMETRY = $geometry_literal
FRAGMENTS = $fragments_literal
LEVEL_TO_SPEC = $level_map_literal
SPEC_KEYWORDS = $spec_keywords_literal
BSSE_TYPES = $bsse_literal
MAX_NBODY = $max_nbody
RETURN_TOTAL_DATA = $return_total_data
SUPERSYSTEM_IE_ONLY = $supersystem_ie_only
TOTAL_CHARGE = $total_charge
TOTAL_MULTIPLICITY = $total_multiplicity


def create_molecule() -> Molecule:
    """Construct the fragmented molecular system."""

    return Molecule(
        symbols=SYMBOLS,
        geometry=GEOMETRY,
        fragments=FRAGMENTS,
        molecular_charge=TOTAL_CHARGE,
        molecular_multiplicity=TOTAL_MULTIPLICITY,
    )


def create_specifications() -> Dict[str, AtomicSpecification]:
    """Build per-level atomic specifications."""

    return {
        label: AtomicSpecification(
            program="$qc_program",
            driver="$driver",
            model={"method": "$method", "basis": "$basis"},
            keywords=dict(SPEC_KEYWORDS),
            protocols={"stdout": True},
            extras={},
        )
        for label in set(LEVEL_TO_SPEC.values())
    }


def create_manybody_input() -> ManyBodyInput:
    """Create the QCManyBody input model for the calculation."""

    molecule = create_molecule()
    keywords = ManyBodyKeywords(
        bsse_type=BSSE_TYPES,
        max_nbody=MAX_NBODY,
        return_total_data=RETURN_TOTAL_DATA,
        supersystem_ie_only=SUPERSYSTEM_IE_ONLY,
        levels=LEVEL_TO_SPEC,
    )

    specification = ManyBodySpecification(
        driver="$driver",
        keywords=keywords,
        specification=create_specifications(),
    )

    return ManyBodyInput(molecule=molecule, specification=specification)


def create_parallel_config() -> ParallelConfig:
    """Return the parallel execution configuration."""

    return ParallelConfig(
        max_workers=$max_workers,
        execution_mode="$execution_mode",
        use_qcengine=$use_qcengine,
        qc_program="$qc_program",
        basis_set="$basis",
        memory_limit_mb=$memory_limit_mb,
        timeout_seconds=$timeout_seconds,
        qcengine_config={
            "keywords": dict(SPEC_KEYWORDS),
            "protocols": {"stdout": False},
        },
    )


def run_parallel_calculation():
    """Execute the many-body calculation using the parallel executor."""

    manybody_input = create_manybody_input()
    config = create_parallel_config()
    executor = ParallelManyBodyExecutor.from_manybodyinput(
        manybody_input,
        config,
    )

    start = time.time()
    fragment_results = executor.execute_full_calculation()
    stats = executor.get_execution_statistics()
    elapsed = time.time() - start

    print("QCManyBody parallel calculation complete!")
    total_fragments = stats.get("total_fragments")
    if total_fragments is not None:
        print(f"Total fragments executed: {total_fragments}")
    levels_executed = stats.get("levels_executed")
    if levels_executed is not None:
        print(f"Levels processed: {levels_executed}")
    parallel_time = stats.get("parallel_time")
    if parallel_time is not None:
        print(f"Parallel wall time: {parallel_time:.2f} s")
    serial_time = stats.get("serial_time")
    if serial_time:
        print(f"Serial baseline: {serial_time:.2f} s")
    speedup_factor = stats.get("speedup_factor")
    if speedup_factor:
        print(f"Measured speedup: {speedup_factor:.2f}x")
    print(f"Elapsed (wall clock): {elapsed:.2f} s")

    analysis = executor.core.analyze(fragment_results)
    results = analysis.get("results", {}) if isinstance(analysis, dict) else {}
    if results:
        total_key = f"nocp_corrected_total_energy_through_{MAX_NBODY}_body"
        interaction_key = (
            f"nocp_corrected_interaction_energy_through_{MAX_NBODY}_body"
        )
        total_energy = results.get(total_key)
        interaction_energy = results.get(interaction_key)
        if total_energy is not None:
            print(
                f"Total energy (through {MAX_NBODY}-body): "
                f"{float(total_energy):.12f} Eh"
            )
        if interaction_energy is not None:
            print(
                f"Interaction energy (through {MAX_NBODY}-body): "
                f"{float(interaction_energy):.12f} Eh"
            )

    return analysis


if __name__ == "__main__":
    run_parallel_calculation()
'''
        )
    )

    script_content = script_template.substitute(
        symbols_literal=symbols_literal,
        geometry_literal=geometry_literal,
        fragments_literal=fragments_literal,
        level_map_literal=level_map_literal,
        spec_keywords_literal=spec_keywords_literal,
        bsse_literal=bsse_literal,
        max_nbody=max_nbody,
        return_total_data=return_total_data,
        supersystem_ie_only=supersystem_ie_only,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        qc_program=qc_program,
        driver=driver,
        method=method,
        basis=basis,
        max_workers=max_workers,
        execution_mode=execution_mode,
        use_qcengine=use_qcengine,
        memory_limit_mb=memory_limit_mb,
        timeout_seconds=timeout_seconds,
    )

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write(script_content)

    print("Conversion complete!")
    print(f"Input file: {hmbe_file}")
    print(f"Output file: {output_file}")
    print(f"System: {len(symbols)} atoms in {len(fragments)} fragments")
    print(f"Method: {method}/{basis}")
    print(f"Max n-body: {max_nbody}")
    print(
        f"Parallel execution mode: {execution_mode} "
        f"with {max_workers} workers"
    )

    return output_file


def main() -> None:
    """Entry point for command-line usage."""

    parser = argparse.ArgumentParser(
        description="Convert HMBE input files to QCManyBody helper scripts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """Examples::

    python hmbe_to_qcmanybody.py water_000.inp
    python hmbe_to_qcmanybody.py water_000.inp my_output.py --max-workers 8 \
        --execution-mode multiprocessing
            """
        ),
    )

    parser.add_argument("input_file", help="HMBE input file (.inp)")
    parser.add_argument(
        "output_file",
        nargs="?",
        help="Output Python script path",
    )

    parser.add_argument(
        "--bsse-type",
        action="append",
        dest="bsse_types",
        help="BSSE treatments to include (repeatable)",
    )
    parser.add_argument(
        "--no-return-total-data",
        action="store_false",
        dest="return_total_data",
        help="Disable total data reporting in ManyBodyKeywords",
    )
    parser.add_argument(
        "--supersystem-ie-only",
        action="store_true",
        help="Enable supersystem interaction-energy-only mode",
    )
    parser.add_argument(
        "--execution-mode",
        choices=["multiprocessing", "threading", "serial"],
        default="multiprocessing",
        help="Parallel execution mode for ParallelConfig",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Maximum number of parallel workers",
    )
    parser.add_argument(
        "--memory-mb",
        type=int,
        default=1000,
        help="Per-worker memory limit in megabytes",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=7200,
        help="Timeout per fragment calculation (seconds)",
    )
    parser.add_argument(
        "--no-qcengine",
        action="store_false",
        dest="use_qcengine",
        help="Disable QCEngine execution in the generated ParallelConfig",
    )
    parser.set_defaults(return_total_data=True, use_qcengine=True)

    parser.add_argument(
        "--qc-program",
        dest="qc_program_override",
        help="Override QC program in generated script",
    )
    parser.add_argument(
        "--basis-set",
        dest="basis_override",
        help="Override the basis set used in the output script",
    )

    args = parser.parse_args()

    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file '{input_path}' not found!", file=sys.stderr)
        sys.exit(1)

    try:
        convert_hmbe_to_qcmanybody(
            hmbe_file=str(input_path),
            output_file=args.output_file,
            bsse_types=args.bsse_types or ["nocp"],
            return_total_data=args.return_total_data,
            supersystem_ie_only=args.supersystem_ie_only,
            execution_mode=args.execution_mode,
            max_workers=args.max_workers,
            use_qcengine=args.use_qcengine,
            qc_program_override=args.qc_program_override,
            basis_override=args.basis_override,
            memory_limit_mb=args.memory_mb,
            timeout_seconds=args.timeout_seconds,
        )
    except Exception as exc:  # pragma: no cover - surface CLI diagnostics
        print(f"Error during conversion: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
