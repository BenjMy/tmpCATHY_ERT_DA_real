# =============================================================================
# ERT (pygimli) frame -> CATHY frame, shared by nb6 (truth) and nb7 (DA)
# =============================================================================
# Factored out of nb5 (choose_ert_frame): same candidates, same alignment gate.
import numpy as np
from nc_model import X, Y, ELEVATION


def choose_ert_frame(elecs, mesh_points, cathy_bounds, z_sign=None, mesh_shift=None, auto_align=True):
    """Return (signs, shift, origin): the transform p_cathy = p_ert * signs + shift.

    Every candidate is tested with check_ert_cathy_alignment on the electrodes AND the pygimli mesh
    nodes. z_sign None = try +1 and -1; mesh_shift None = legacy shift of farming_DA_ERT3d.py
    [X/2-0.5, Y/2-0.5, +ELEVATION], then (auto_align) the shift that centres the electrode array."""
    from pyCATHY.ERT.alignment import (suggest_frame, transform_points, check_ert_cathy_alignment,
                                       ERTAlignmentError)
    z_signs = (1, -1) if z_sign is None else (z_sign,)
    cands = []
    if mesh_shift is not None:
        cands += [((1, 1, sz), tuple(mesh_shift), 'manual') for sz in z_signs]
    else:
        legacy = (X / 2 - 0.5, Y / 2 - 0.5, ELEVATION)
        cands += [((1, 1, sz), legacy, 'legacy (farming_DA_ERT3d.py)') for sz in z_signs]
        if auto_align:
            sg, sh, _ov = suggest_frame(elecs, mesh_points, cathy_bounds, z_signs)
            cands.append((sg, sh, 'auto (array centred on the plot, surface = CATHY top)'))
    log = []
    for signs, shift, origin in cands:
        rep = check_ert_cathy_alignment(
            transform_points(elecs, signs, shift), cathy_bounds,
            pg_points=transform_points(mesh_points, signs, shift), on_fail='ignore')
        if rep['ok']:
            return signs, shift, origin
        log.append(f"  * {origin}, signs={signs}, shift={np.round(shift, 3).tolist()}: "
                   + "; ".join(rep['problems']))
    raise ERTAlignmentError("no ERT->CATHY frame passes the alignment check; tried:\n" + "\n".join(log))
