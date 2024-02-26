import math


def compare(a, b):
    if not math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-7):
        abs_diff = abs(a - b)
        rel_diff = abs_diff / min(abs(a), abs(b))
        raise RuntimeError("Not close: {} vs {} rel_diff={} abs_diff={}".format(a, b, rel_diff, abs_diff))
