"""
QCManyBody CLI Plan Command

Shows the execution plan for a calculation without running it.
"""

import logging
from argparse import Namespace

logger = logging.getLogger(__name__)


def handle_plan(args: Namespace) -> int:
    """
    Handle the 'plan' command.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments

    Returns
    -------
    int
        Exit code (0 for success, non-zero for failure)
    """
    logger.info("Generating execution plan")
    logger.info(f"Input file: {args.input}")

    # Step 1: Load and validate input file
    try:
        from qcmanybody.cli.input_parser import parse_input_file

        logger.debug("Parsing input file...")
        cli_input = parse_input_file(args.input)
        logger.info(f"✓ Input file parsed successfully")

    except Exception as e:
        logger.error(f"Failed to parse input file: {e}")
        return 1

    # Step 2: Convert to ManyBodyInput
    try:
        from qcmanybody.cli.converter import convert_to_manybody_input

        logger.debug("Converting to ManyBodyInput...")
        mb_input = convert_to_manybody_input(cli_input, input_file_path=args.input)
        logger.info(f"✓ Input converted successfully")

    except Exception as e:
        logger.error(f"Failed to convert input: {e}")
        return 1

    # Step 3: Build task list using builder (doesn't require QC programs)
    try:
        from qcmanybody.builder import build_nbody_compute_list

        logger.debug("Building task list...")

        # Build list of n-body levels to compute
        max_nbody = mb_input.specification.keywords.max_nbody
        nfragments = len(mb_input.molecule.fragments)
        if max_nbody is None:
            max_nbody = nfragments
        nbodies = list(range(1, max_nbody + 1))

        # Build compute list (this gives us fragment/basis combinations, not full task specs)
        compute_dict = build_nbody_compute_list(
            bsse_type=mb_input.specification.keywords.bsse_type,
            nfragments=nfragments,
            nbodies=nbodies,
            return_total_data=mb_input.specification.keywords.return_total_data,
            supersystem_ie_only=mb_input.specification.keywords.supersystem_ie_only,
        )

        # Count total number of tasks across all BSSE types
        total_tasks = 0
        for bsse_type, level_dict in compute_dict.items():
            if bsse_type in ["all", "cp", "nocp", "vmfc_compute"]:
                for level, task_set in level_dict.items():
                    total_tasks += len(task_set)

        logger.info(f"✓ Execution plan generated ({total_tasks} tasks)")

    except Exception as e:
        logger.error(f"Failed to generate plan: {e}")
        logger.debug("Full error:", exc_info=True)
        return 1

    # Step 4: Display execution plan
    print()
    print("=" * 80)
    print(f"Execution Plan: {args.input}")
    print("=" * 80)
    print()

    # Molecular system info
    mol = mb_input.molecule
    print(f"Molecular System:")
    print(f"  Atoms: {len(mol.symbols)}")
    print(f"  Fragments: {len(mol.fragments)}")
    for i, frag in enumerate(mol.fragments, 1):
        atoms = [mol.symbols[idx] for idx in frag]
        print(f"    Fragment {i}: {atoms} (atoms {frag})")
    print()

    # Calculation info
    spec = mb_input.specification
    print(f"Calculation Settings:")
    print(f"  Driver: {spec.driver.value}")
    print(f"  BSSE types: {[bt.value for bt in spec.keywords.bsse_type]}")
    print(f"  Max n-body: {spec.keywords.max_nbody}")

    # Method/basis info
    if len(spec.specification) == 1 and "(auto)" in spec.specification:
        # Single-level
        atomic_spec = spec.specification["(auto)"]
        print(f"  Method: {atomic_spec.model.method}")
        print(f"  Basis: {atomic_spec.model.basis}")
        print(f"  Program: {atomic_spec.program}")
    else:
        # Multi-level
        print(f"  Levels: {len(spec.specification)}")
        for level_key, atomic_spec in spec.specification.items():
            print(f"    {level_key}: {atomic_spec.model.method}/{atomic_spec.model.basis}")
    print()

    # Task list
    print(f"Computational Tasks:")
    print(f"  Total tasks: {total_tasks}")

    if args.show_tasks:
        print()
        # Display tasks by BSSE type and n-body level
        bsse_types_to_show = mb_input.specification.keywords.bsse_type
        for bsse_type in bsse_types_to_show:
            bsse_key = bsse_type.value
            if bsse_key in compute_dict and compute_dict[bsse_key]:
                print(f"  {bsse_key.upper()} tasks:")
                level_dict = compute_dict[bsse_key]
                for nbody in sorted(level_dict.keys()):
                    task_set = level_dict[nbody]
                    print(f"    {nbody}-body: {len(task_set)} tasks")
                    # Show first few tasks
                    for i, task in enumerate(list(task_set)[:3]):
                        frag, bas = task
                        print(f"      Task {i+1}: fragments {frag} in basis {bas}")
                    if len(task_set) > 3:
                        print(f"      ... and {len(task_set) - 3} more")
                print()

    print("=" * 80)
    print()

    logger.info("✓ Plan command completed successfully")
    return 0
