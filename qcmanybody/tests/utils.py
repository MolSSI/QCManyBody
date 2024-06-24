import json
import math
import os
from typing import Any, Iterable, Literal, Mapping, Optional, Union

import numpy
import zstandard
from qcelemental.models import AtomicInput, Molecule

from qcmanybody import ManyBodyCore, delabeler
from qcmanybody.models import BsseEnum
from qcmanybody.utils import translate_qcvariables

_my_dir = os.path.dirname(os.path.realpath(__file__))


def compare(a, b):
    if isinstance(a, float) and isinstance(b, float):
        if not math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-7):
            abs_diff = abs(a - b)
            rel_diff = abs_diff / min(abs(a), abs(b))
            raise RuntimeError("Not close: {} vs {} rel_diff={} abs_diff={}".format(a, b, rel_diff, abs_diff))
    elif isinstance(a, (numpy.ndarray, list)) and isinstance(b, (numpy.ndarray, list)):
        a = numpy.array(a).ravel()
        b = numpy.array(b).ravel()
        if a.size != b.size:
            raise RuntimeError(f"Size mismatch: {a.size} vs {b.size}")
        if not numpy.allclose(a, b, rtol=1e-7, atol=1e-7):
            raise RuntimeError("Not close: {} vs {}".format(a, b))
    else:
        raise RuntimeError(f"Unknown types: {type(a)} vs {type(b)}")


def jsonify(data):
    # Make numpy arrays into nested lists
    if isinstance(data, dict):
        return {k: jsonify(v) for k, v in data.items()}
    elif isinstance(data, (list, tuple)):
        return [jsonify(x) for x in data]
    elif isinstance(data, numpy.ndarray):
        return jsonify(data.tolist())
    else:
        return data


def unjsonify(data):
    # Make nested lists into numpy arrays
    if isinstance(data, dict):
        return {k: unjsonify(v) for k, v in data.items()}
    elif isinstance(data, (list, tuple, numpy.ndarray)):
        return numpy.array(data)
    else:
        return data


def load_ref_data(file_base):
    filepath = os.path.join(_my_dir, "ref_data", file_base + ".json.zst")
    with zstandard.open(filepath, "rt") as f:
        return unjsonify(json.load(f))


def load_component_data(file_base):
    filepath = os.path.join(_my_dir, "component_data", file_base + ".json.zst")
    with zstandard.open(filepath, "rt") as f:
        return unjsonify(json.load(f))


def generate_component_data(
    mol,
    levels,
    specifications,
    bsse_type,
    return_total_data,
    out_filename,
    supsersytem_ie_only=False,
    embedding_charges=None,
):
    mbc, component_results = run_qcengine(
        specifications, mol, bsse_type, levels, return_total_data, supsersytem_ie_only, embedding_charges
    )

    component_results = jsonify(component_results)
    filepath = os.path.join(_my_dir, "component_data", out_filename + ".json.zst")
    with zstandard.open(filepath, "wt") as f:
        json.dump(component_results, f, indent=2)


def compare_results(qcmb_results, ref_results, levels):
    import pprint

    print("*" * 80)
    pprint.pprint(qcmb_results)
    print("-" * 80)
    pprint.pprint(ref_results)
    print("*" * 80)

    print("COMPARING CURRENT ENERGY")
    compare(qcmb_results["ret_energy"], ref_results["CURRENT ENERGY"])

    levels_no_ss = [k for k in levels.keys() if k != "supersystem"]
    for b in list(BsseEnum):
        for l in levels_no_ss:
            # for "1CP"
            # qcmb is also broken down into dicts
            qcmb_key = f"{l}{b.value}"
            ref_key = qcmb_key.upper()
            print("COMPARING ", qcmb_key, " vs ", ref_key)
            # compare(qcmb_results["energy_body_dict"][b.value][l], ref_results[ref_key])
            compare(qcmb_results["energy_body_dict"][qcmb_key], ref_results[ref_key])

    if "ret_gradient" in qcmb_results and "CURRENT GRADIENT" in ref_results:
        print("COMPARING CURRENT GRADIENT")
        compare(qcmb_results["ret_gradient"], ref_results["CURRENT GRADIENT"])

        for b in list(BsseEnum):
            for l in levels_no_ss:
                # for example, "GRADIENT 1CP"
                # qcmb is broken down into dicts
                ref_key = f"GRADIENT {l}{b.value.upper()}"
                if ref_key in ref_results:
                    print("COMPARING ", ref_key)
                    compare(qcmb_results["gradient_body_dict"][b.value][l], ref_results[ref_key])

    if "ret_hessian" in qcmb_results and "CURRENT HESSIAN" in ref_results:
        print("COMPARING CURRENT HESSIAN")
        compare(qcmb_results["ret_hessian"], ref_results["CURRENT HESSIAN"])

        for b in list(BsseEnum):
            for l in levels_no_ss:
                # for example, "GRADIENT 1CP"
                # qcmb is broken down into dicts
                ref_key = f"HESSIAN {l}{b.value.upper()}"
                if ref_key in ref_results:
                    print("COMPARING ", ref_key)
                    compare(qcmb_results["hessian_body_dict"][b.value][l], ref_results[ref_key])

    ###########################################################

    res = qcmb_results["results"]
    if not res:
        return

    # res (from ManyBodyCore) keys are ManyBodyResultProperties fields, while ref_results keys are
    #   Psi4.core.Wavefunction.variables() keys, so need to translate former for compare() calls below
    res = translate_qcvariables(res)

    if not f"NOCP-CORRECTED TOTAL ENERGY" in ref_results:
        # Psi4 used during the bootstrapping tests does not have data for multi+ss
        return

    for b in list(BsseEnum):
        bstr = b.value.upper()

        print("COMPARING ", f"{bstr}-CORRECTED TOTAL ENERGY")
        compare(res[f"{bstr}-CORRECTED TOTAL ENERGY"], ref_results[f"{bstr}-CORRECTED TOTAL ENERGY"])
        print("COMPARING ", f"{bstr}-CORRECTED INTERACTION ENERGY")
        compare(res[f"{bstr}-CORRECTED INTERACTION ENERGY"], ref_results[f"{bstr}-CORRECTED INTERACTION ENERGY"])

        for l in levels_no_ss:
            print("COMPARING ", f"{bstr}-CORRECTED TOTAL ENERGY THROUGH {l}-BODY")
            compare(
                res[f"{bstr}-CORRECTED TOTAL ENERGY THROUGH {l}-BODY"],
                ref_results[f"{bstr}-CORRECTED TOTAL ENERGY THROUGH {l}-BODY"],
            )

            if l > 1:
                print("COMPARING ", f"{bstr}-CORRECTED INTERACTION ENERGY THROUGH {l}-BODY")
                compare(
                    res[f"{bstr}-CORRECTED INTERACTION ENERGY THROUGH {l}-BODY"],
                    ref_results[f"{bstr}-CORRECTED INTERACTION ENERGY THROUGH {l}-BODY"],
                )

                print("COMPARING ", f"{bstr}-CORRECTED {l}-BODY CONTRIBUTION TO ENERGY")
                compare(
                    res[f"{bstr}-CORRECTED {l}-BODY CONTRIBUTION TO ENERGY"],
                    ref_results[f"{bstr}-CORRECTED {l}-BODY CONTRIBUTION TO ENERGY"],
                )


def run_qcengine(
    specifications: Mapping[str, Mapping[str, Any]],
    molecule: Molecule,
    bsse_type: Iterable[BsseEnum],
    levels: Mapping[Union[int, Literal["supersystem"]], str],
    return_total_data: bool,
    supersystem_ie_only: bool,
    embedding_charges: Optional[Mapping[int, list]],
):
    import qcengine as qcng

    mbc = ManyBodyCore(
        molecule,
        bsse_type,
        levels,
        return_total_data=return_total_data,
        supersystem_ie_only=supersystem_ie_only,
        embedding_charges=embedding_charges,
    )

    component_results = {}

    computation_count = {}
    for chem, label, imol in mbc.iterate_molecules():
        print(label)
        inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])

        _, real, bas = delabeler(label)
        computation_count.setdefault(len(real), 0)
        computation_count[len(real)] += 1

        result = qcng.compute(inp, specifications[chem]["program"])

        if not result.success:
            print(result.error.error_message)
            raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        # pull out stuff
        props = {"energy", "gradient", "hessian"}

        component_results[label] = {}

        for p in props:
            if hasattr(result.properties, f"return_{p}"):
                v = getattr(result.properties, f"return_{p}")
                # print(f"  {label} {p}: {v}")
                if v is not None:
                    component_results[label][p] = v

    return mbc, component_results
