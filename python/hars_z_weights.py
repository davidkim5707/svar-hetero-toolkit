"""Within-theta estimator weights for HARS-Z output.

The HARS-Z estimator of a posterior summary is

    E_hat[h] = (1/D) * sum_{d=1..D} [ sum_r v_{d,r} h(theta_d, Q^(d,r))
                                       / sum_r v_{d,r} ],

where d indexes the theta blocks, r the rotations retained within a block,
and v is the volume element of the zero-restriction construction. This
module rebuilds the weights of that estimator from a saved ``output`` struct
and matches ``local_block_normalize`` and ``local_assemble_weights`` in
``core/harsz.m``. Take the weight vector of every pooled summary from
:func:`estimator_weights`.

Weight fields in the saved ``.mat``

``output.draw_weight``
    Raw volume element v(Q; theta) per retained draw. Never use raw.
``output.draw_block``
    Theta-block id (the index d above). When the field is absent, blocks
    are recovered from run-lengths of ``output.draw_lh``, because all
    retained draws of one sweep share theta and hence the same lh.
``output.draw_weightA`` / ``output.draw_weight_final``
    Precomputed by the sampler. :func:`estimator_weights` recomputes them
    from the raw fields and, when both are present, asserts agreement, so
    the normalization is verifiable from the raw ``.mat`` alone.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = [
    "block_ids_from_lh",
    "block_normalize",
    "assemble_weights",
    "estimator_weights",
]


def _has(output, name: str) -> bool:
    """Field presence for both scipy.io.loadmat conventions.

    Supports ``struct_as_record=False`` (attribute access) and the default
    record-array loading (``output['name'].item()`` with ``squeeze_me=True``).
    """
    if hasattr(output, "dtype") and getattr(output.dtype, "names", None):
        return name in output.dtype.names
    return hasattr(output, name)


def _get(output, name: str) -> np.ndarray:
    if hasattr(output, "dtype") and getattr(output.dtype, "names", None):
        f = output[name]
        return np.asarray(f.item() if getattr(f, "size", None) == 1 else f)
    return np.asarray(getattr(output, name))


def block_ids_from_lh(lh: np.ndarray) -> np.ndarray:
    """Recover theta-block ids from run-lengths of the stored log posterior.

    Fallback for ``.mat`` files without ``output.draw_block``.
    Retained draws of one Step-3 sweep share theta and hence ``lh``;
    consecutive equal values are one block. Caveat, documented in the
    sampler: if a Metropolis rejection repeats theta across consecutive
    sweeps with an identical stored ``lh``, this merges the two sweeps into
    one block. Output of ``core/harsz.m`` carries the explicit id and never
    hits this path.
    """
    lh = np.asarray(lh, dtype=float).ravel()
    if lh.size == 0:
        return np.zeros(0, dtype=int)
    change = np.r_[False, np.diff(lh) != 0]
    return np.cumsum(change).astype(int)


def block_normalize(v: np.ndarray, blk: np.ndarray) -> np.ndarray:
    """Within-theta normalization w_A(d, r) = v_{d,r} / sum_r' v_{d,r'}.

    Every theta block carries the same total mass, so a weighted mean under
    the result equals the equal-weighted average over theta of the per-block
    ratio estimators (Corollary ratio_identity), and it is invariant to any
    per-theta constant rescaling of v. Blocks whose weights do not sum to a
    positive finite number (volume-element underflow) fall back to uniform
    weights within the block, with a warning, never silently.
    """
    v = np.asarray(v, dtype=float).ravel()
    blk = np.asarray(blk).ravel()
    if v.shape != blk.shape:
        raise ValueError("block id vector must match the weight vector")
    _, gid = np.unique(blk, return_inverse=True)
    bsum = np.bincount(gid, weights=v)
    blen = np.bincount(gid)
    bad = ~(np.isfinite(bsum) & (bsum > 0))
    if bad.any():
        warnings.warn(
            f"{int(bad.sum())} theta block(s) have a non-positive or "
            "non-finite volume-weight sum (underflow/overflow of the volume "
            "element); using uniform weights within those blocks.",
            RuntimeWarning,
            stacklevel=2,
        )
    with np.errstate(invalid="ignore", divide="ignore"):
        wA = v / bsum[gid]
    badd = bad[gid]
    wA[badd] = 1.0 / blen[gid[badd]]
    return wA


def assemble_weights(
    v: np.ndarray, blk: np.ndarray, omega: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(w_final, w_A)`` with ``w_final = block_normalize(v) / omega``.

    The 1/omega factor (AD-RR narrative importance correction, extension
    beyond the paper's formal target) keeps its across-theta variation, so it
    is applied after, and outside of, the within-block normalization of v.
    ``omega is None`` or all ones (``use_nrr_omega = false``, the paper's
    baseline) gives ``w_final == w_A`` exactly.
    """
    wA = block_normalize(v, blk)
    if omega is None:
        return wA.copy(), wA
    omega = np.asarray(omega, dtype=float).ravel()
    if omega.shape != wA.shape:
        raise ValueError("omega vector must match the weight vector")
    return wA / omega, wA


def estimator_weights(output, check_stored: bool = True) -> np.ndarray:
    """Estimator weights of eq. (summary_estimator) for a loaded ``.mat``.

    ``output`` is the ``output`` struct as loaded by ``scipy.io.loadmat``,
    under either loading convention (``struct_as_record=False`` attribute
    access, or the default record array with ``squeeze_me=True``). Uses
    ``output.draw_block`` when present, else the ``draw_lh``
    run-length fallback. When the sampler stored ``draw_weight_final``,
    the recomputed weights are asserted against it (relative tolerance
    1e-10) so a stale or mismatched ``.mat`` fails loudly.
    """
    v = np.asarray(_get(output, "draw_weight"), dtype=float).ravel()
    if _has(output, "draw_block"):
        blk = np.asarray(_get(output, "draw_block")).ravel()
    else:
        warnings.warn(
            "output.draw_block absent; recovering theta blocks "
            "from run-lengths of output.draw_lh.",
            RuntimeWarning,
            stacklevel=2,
        )
        blk = block_ids_from_lh(np.asarray(_get(output, "draw_lh"), dtype=float))
    omega = None
    if _has(output, "draw_omega_nrr"):
        om = np.asarray(_get(output, "draw_omega_nrr"), dtype=float).ravel()
        if om.size == v.size and not np.allclose(om, 1.0):
            omega = om
    w_final, _ = assemble_weights(v, blk, omega)
    if check_stored and _has(output, "draw_weight_final"):
        stored = np.asarray(_get(output, "draw_weight_final"), dtype=float).ravel()
        if not np.allclose(w_final, stored, rtol=1e-10, atol=0.0):
            raise AssertionError(
                "recomputed weights disagree with the stored "
                "output.draw_weight_final; the .mat and this module are out "
                "of sync"
            )
    return w_final