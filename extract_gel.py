#!/usr/bin/env python3
"""
extract_gel.py - STEP 1: extract the band table from ONE GelAnalyzer gel (.gap)

    python3 extract_gel.py "Primer i5-Population2.gap"

Writes, next to the .gap file:
    bands_<gel>.csv          band table (input of step 2, bin_gel.R)
    check_<gel>.png          gel image with detected bands, Rf curves and lane
                             numbers, for a visual check before going on
    calibration_<gel>.json   calibration details (ladder, Rf curves, residuals)

Options:
    --curves auto     (default) Rf curves fitted automatically to rows of bands;
                      corrects tilted / smiling gels
    --curves drawn    use the Rf curves drawn by hand in GelAnalyzer
    --curves none     no correction (Rf = peak position / lane length)
    --out DIR         output folder (default: folder of the .gap)

Standard library only (Python >= 3.8).
"""

import argparse
import csv
import json
import math
import statistics
import struct
import sys
import zlib
from bisect import bisect_right
from pathlib import Path

# Sanity limits used while scanning the binary file
MAX_LANES = 400
MAX_BANDS = 1000
MAX_LENGTH = 100_000
MAX_POINTS = 20_000

# If band rows drift by more than this (px) from one side of the gel to the
# other, they are tracked again with the overall tilt removed.
TILT_THRESHOLD_PX = 4.0


# =========================================================================== #
# 1. Reading the .gap file (big-endian binary)
# =========================================================================== #
class _Reader:
    def __init__(self, raw, pos=0):
        self.raw, self.pos = raw, pos

    def i32(self):
        if self.pos + 4 > len(self.raw):
            raise EOFError
        v = struct.unpack_from(">i", self.raw, self.pos)[0]
        self.pos += 4
        return v

    def f64(self):
        if self.pos + 8 > len(self.raw):
            raise EOFError
        v = struct.unpack_from(">d", self.raw, self.pos)[0]
        self.pos += 8
        return v


def read_header(raw):
    """Signature, image name and the embedded image bytes."""
    p = 0
    n = struct.unpack_from(">H", raw, p)[0]; p += 2
    signature = raw[p:p + n].decode("latin-1"); p += n
    n = struct.unpack_from(">H", raw, p)[0]; p += 2
    name = raw[p:p + n]; p += n
    try:
        name = name.decode("utf-8")
    except UnicodeDecodeError:
        name = name.decode("latin-1")
    size = struct.unpack_from(">I", raw, p)[0]; p += 4
    image = raw[p:p + size]
    if not signature.startswith("GAP"):
        raise RuntimeError(f"unexpected file signature: {signature!r}")
    return dict(signature=signature, image_name=name, image=image,
                image_end=p + size)


def _read_curves(r):
    n = r.i32()
    if not (0 <= n <= 100):
        raise ValueError
    curves = []
    for _ in range(n):
        m = r.i32()
        if not (1 <= m <= MAX_POINTS):
            raise ValueError
        pts = [(r.i32(), r.i32()) for _ in range(m)]
        rf = r.f64()
        if not (-0.001 <= rf <= 1.001):
            raise ValueError
        if any(b[0] < a[0] for a, b in zip(pts, pts[1:])):
            raise ValueError
        curves.append(dict(rf=rf, pts=pts))
    return curves


def _read_lanes(r):
    n = r.i32()
    if not (1 <= n <= MAX_LANES):
        raise ValueError
    lanes = []
    for _ in range(n):
        if r.i32() != 2:
            raise ValueError
        x, y0, x_bottom = r.i32(), r.i32(), r.i32()
        length, width, _reserved = r.i32(), r.i32(), r.i32()
        if min(x, x_bottom) < 0 or not (0 < length <= MAX_LENGTH) \
                or not (0 < width <= MAX_LENGTH):
            raise ValueError
        nb = r.i32()
        if not (0 <= nb <= MAX_POINTS):
            raise ValueError
        baseline = [(r.i32(), r.i32()) for _ in range(nb)]
        if any(not (-1 <= q <= length + 1) for q, _ in baseline) or \
                any(b[0] < a[0] for a, b in zip(baseline, baseline[1:])):
            raise ValueError
        n_bands = r.i32()
        if not (0 <= n_bands <= MAX_BANDS):
            raise ValueError
        bands = []
        for _ in range(n_bands):
            start, end, peak, flag = r.i32(), r.i32(), r.i32(), r.i32()
            if flag not in (0, 1) or not (-1 <= start <= peak <= end <= length + 1):
                raise ValueError
            mw = r.f64() if flag == 1 else None
            if mw is not None and not (0 < mw < 1e9):
                raise ValueError
            bands.append(dict(start=start, end=end, peak=peak, ladder_mw=mw))
        lanes.append(dict(x=x, y0=y0, x_bottom=x_bottom, length=length,
                          width=width, baseline=baseline, bands=bands))
    if any(b["x"] <= a["x"] for a, b in zip(lanes, lanes[1:])):
        raise ValueError
    return lanes


def parse_gap(path):
    """Return dict(header, curves, lanes).

    The curve block starts shortly after the embedded image at an offset that
    varies between files, so offsets are tried one by one and a parse is
    accepted only if it ends (almost) exactly at the end of the file.
    """
    raw = Path(path).read_bytes()
    header = read_header(raw)
    start = header["image_end"]
    for st in range(start, min(start + 4096, len(raw) - 8)):
        try:
            r = _Reader(raw, st)
            curves = _read_curves(r)
            lanes = _read_lanes(r)
        except (ValueError, EOFError, struct.error):
            continue
        if len(raw) - r.pos <= 32:
            return dict(header=header, curves=curves, lanes=lanes)
    raise RuntimeError("curve/lane structure not recognised in this .gap file")


# =========================================================================== #
# 2. Embedded image (minimal PNG decoder, standard library only)
# =========================================================================== #
def decode_png(data):
    """Return (width, height, rows) as grey-scale intensities; colour images
    are converted to the mean of the RGB channels."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("embedded image is not a PNG")
    p, idat, ihdr, plte = 8, [], None, None
    while p < len(data):
        n = struct.unpack_from(">I", data, p)[0]
        kind = data[p + 4:p + 8]
        body = data[p + 8:p + 8 + n]
        p += 12 + n
        if kind == b"IHDR":
            ihdr = struct.unpack(">IIBBBBB", body)
        elif kind == b"PLTE":
            plte = body
        elif kind == b"IDAT":
            idat.append(body)
        elif kind == b"IEND":
            break
    w, h, depth, colour, _c, _f, interlace = ihdr
    if interlace:
        raise RuntimeError("interlaced PNG not supported")
    channels = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}[colour]
    bpp_bits = depth * channels
    bpp = max(1, bpp_bits // 8)
    stride = (w * bpp_bits + 7) // 8
    raw = zlib.decompress(b"".join(idat))
    rows, prev = [], bytearray(stride)
    q = 0
    for _ in range(h):
        filt = raw[q]; q += 1
        row = bytearray(raw[q:q + stride]); q += stride
        if filt == 1:                                   # Sub
            for i in range(bpp, stride):
                row[i] = (row[i] + row[i - bpp]) & 255
        elif filt == 2:                                 # Up
            for i in range(stride):
                row[i] = (row[i] + prev[i]) & 255
        elif filt == 3:                                 # Average
            for i in range(stride):
                a = row[i - bpp] if i >= bpp else 0
                row[i] = (row[i] + ((a + prev[i]) >> 1)) & 255
        elif filt == 4:                                 # Paeth
            for i in range(stride):
                a = row[i - bpp] if i >= bpp else 0
                b = prev[i]
                c = prev[i - bpp] if i >= bpp else 0
                pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
                pr = a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
                row[i] = (row[i] + pr) & 255
        rows.append(row)
        prev = row
    grey = []
    for row in rows:
        if depth == 8 and channels == 1 and colour == 0:
            grey.append(list(row))
        elif depth == 16 and colour == 0:
            grey.append([(row[2 * i] << 8) | row[2 * i + 1] for i in range(w)])
        elif depth == 8 and colour == 3:
            grey.append([sum(plte[3 * v:3 * v + 3]) / 3 for v in row])
        elif depth == 8:
            k = min(3, channels)
            grey.append([sum(row[channels * i:channels * i + k]) / k for i in range(w)])
        else:
            raise RuntimeError(f"PNG with bit depth {depth} and colour type "
                               f"{colour} not supported")
    return w, h, grey


# =========================================================================== #
# 3. Geometry: lanes and Rf curves
# =========================================================================== #
def _interp(xs, ys, x):
    """Piecewise-linear interpolation, constant outside the range."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    i = bisect_right(xs, x) - 1
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])


def image_position(lane, pos):
    """Image (x, y) of a point `pos` px below the top of the lane."""
    t = pos / lane["length"]
    return (lane["x"] + t * (lane["x_bottom"] - lane["x"]),
            lane["y0"] + t * lane["length"])


def curve_y(curve, x):
    xs = [p[0] for p in curve["pts"]]
    ys = [p[1] for p in curve["pts"]]
    if len(xs) == 1:
        return ys[0]
    return _interp(xs, ys, x)


def rf_from_curves(curves, x, y):
    """Rf of image point (x, y), interpolated between the Rf curves directly
    above and below it (as GelAnalyzer does). Outside the outermost curves it
    extrapolates from the two nearest ones."""
    cs = sorted(((curve_y(c, x), c["rf"]) for c in curves), key=lambda t: t[1])
    if len(cs) < 2:
        return None
    for (ya, ra), (yb, rb) in zip(cs, cs[1:]):
        if ya <= y <= yb and yb > ya:
            return ra + (rb - ra) * (y - ya) / (yb - ya)
    (ya, ra), (yb, rb) = (cs[0], cs[1]) if y < cs[0][0] else (cs[-2], cs[-1])
    return ra + (rb - ra) * (y - ya) / (yb - ya) if yb != ya else None


# =========================================================================== #
# 4. Molecular-weight calibration
# =========================================================================== #
def calibration_interp(rf, mw):
    """Piecewise-linear interpolation of log10(MW) against Rf. Outside the
    ladder range the end segment is extrapolated (bands are flagged)."""
    pairs = sorted(zip(rf, mw))
    xs = [p[0] for p in pairs]
    ls = [math.log10(p[1]) for p in pairs]
    if any(b <= a for a, b in zip(xs, xs[1:])):
        raise RuntimeError("ladder bands with repeated or out-of-order Rf")

    def f(r):
        if r <= xs[0]:
            i = 0
        elif r >= xs[-1]:
            i = len(xs) - 2
        else:
            i = bisect_right(xs, r) - 1
        t = (r - xs[i]) / (xs[i + 1] - xs[i])
        return 10 ** (ls[i] + t * (ls[i + 1] - ls[i]))
    return f


def calibration_exp(rf, mw):
    """Global model MW = a*exp(-b*Rf) + c (a, c by least squares for each b;
    b by grid search followed by golden-section refinement)."""
    def fit_ac(b):
        e = [math.exp(-b * r) for r in rf]
        n, see, se, sy = len(rf), sum(v * v for v in e), sum(e), sum(mw)
        sey = sum(v * y for v, y in zip(e, mw))
        det = see * n - se * se
        if abs(det) < 1e-300:
            return None
        a = (sey * n - se * sy) / det
        c = (see * sy - se * sey) / det
        return a, c, sum((y - (a * v + c)) ** 2 for v, y in zip(e, mw))

    best = min((fit_ac(k * 0.05) or (0, 0, float("inf")))[2:] + (k * 0.05,)
               for k in range(1, 801))
    lo, hi = max(1e-6, best[1] - 0.05), best[1] + 0.05
    phi = (math.sqrt(5) - 1) / 2
    for _ in range(200):
        b1, b2 = hi - phi * (hi - lo), lo + phi * (hi - lo)
        if fit_ac(b1)[2] < fit_ac(b2)[2]:
            hi = b2
        else:
            lo = b1
    b = (lo + hi) / 2
    a, c, _ = fit_ac(b)
    return (lambda r: a * math.exp(-b * r) + c), dict(a=a, b=b, c=c)


# =========================================================================== #
# 5. Densitometry from the embedded image
# =========================================================================== #
def lane_profile(img, lane):
    """Mean intensity across the lane width at every pixel along the lane."""
    w, h, pix = img
    half = lane["width"] / 2
    out = []
    for pos in range(lane["length"] + 1):
        x, y = image_position(lane, pos)
        yi = int(round(y))
        if not (0 <= yi < h):
            out.append(0.0)
            continue
        x0, x1 = max(0, int(round(x - half))), min(w, int(round(x + half)) + 1)
        out.append(sum(pix[yi][x0:x1]) / max(1, x1 - x0))
    return out


def lane_background(profile, lane, window=31):
    """Baseline stored in the .gap (if any), otherwise a rolling-minimum
    morphological opening of the profile."""
    n = len(profile)
    if len(lane["baseline"]) >= 2:
        xs = [q for q, _ in lane["baseline"]]
        ys = [v for _, v in lane["baseline"]]
        return [_interp(xs, ys, i) for i in range(n)], "gelanalyzer"
    k = window // 2
    mn = [min(profile[max(0, i - k):i + k + 1]) for i in range(n)]
    mx = [max(mn[max(0, i - k):i + k + 1]) for i in range(n)]
    return mx, "rolling_minimum"


# =========================================================================== #
# 6. Automatic Rf curves (--curves auto)
#    Instead of hand-drawn curves, each "row" of bands (the same band in
#    neighbouring lanes) is tracked lane by lane, a smooth curve y(x) is fitted
#    to it, and these curves are used as Rf reference curves.
# =========================================================================== #
def _polyfit(xs, ys, degree):
    """Least squares of degree 1 or 2 without numpy. Returns [c0, c1, (c2)]
    for y = c0 + c1*x + c2*x^2, or None if singular."""
    n = degree + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for c in range(n):                                  # Gauss-Jordan
        p = max(range(c, n), key=lambda r: abs(A[r][c]))
        A[c], A[p] = A[p], A[c]
        b[c], b[p] = b[p], b[c]
        if abs(A[c][c]) < 1e-12:
            return None
        for r in range(n):
            if r != c:
                f = A[r][c] / A[c][c]
                A[r] = [a - f * q for a, q in zip(A[r], A[c])]
                b[r] -= f * b[c]
    return [b[i] / A[i][i] for i in range(n)]


class BandRow:
    """Curve y(x) fitted to one row of bands: a straight line, or a parabola
    when the curvature is real. Outside the lanes that support the row it is
    continued as a straight line with the row's mean slope, so a parabola does
    not shoot off when extrapolated to the ladder lane."""

    def __init__(self, pts):
        self.pts = sorted(pts)
        xs = [p[0] for p in self.pts]
        ys = [p[1] for p in self.pts]
        self.x0, self.x1 = xs[0], xs[-1]
        self.scale = max(1.0, self.x1 - self.x0)
        u = [(x - self.x0) / self.scale for x in xs]
        lin = _polyfit(u, ys, 1) or [sum(ys) / len(ys), 0.0]
        self.linear_slope = lin[1]
        self.c = lin + [0.0]
        self.degree = 1
        self.lanes = []
        sd_line = self._sd(u, ys)
        # parabola only with enough points and a clearly better fit than the
        # line (real smiling, not 1 px of marking noise)
        if len(xs) >= 8:
            q = _polyfit(u, ys, 2)
            if q:
                c_line, self.c = self.c, q
                if self._sd(u, ys) < sd_line - 0.3:
                    self.degree = 2
                else:
                    self.c = c_line
        self.residuals = [y - self._p(v) for v, y in zip(u, ys)]

    def _sd(self, u, ys):
        r = [y - self._p(v) for v, y in zip(u, ys)]
        return (sum(v * v for v in r) / max(1, len(r) - 1)) ** 0.5

    def _p(self, u):
        return self.c[0] + self.c[1] * u + self.c[2] * u * u

    def y(self, x):
        u = (x - self.x0) / self.scale
        if u < 0:
            return self._p(0) + self.linear_slope * u
        if u > 1:
            return self._p(1) + self.linear_slope * (u - 1)
        return self._p(u)

    def residual_sd(self):
        r = self.residuals
        return (sum(v * v for v in r) / max(1, len(r) - 1)) ** 0.5


def _track_rows(lanes, ladders, tol_px=3.0, min_frac=0.5, tilt=0.0, corrected=False):
    """Track rows of bands lane by lane (sample lanes only). A row is accepted
    if it appears in >= min_frac of the sample lanes (and at least 3).

    tilt: overall gel slope (px per px) removed before tracking; points are
    returned in original image coordinates."""
    samples = [(i, L) for i, L in enumerate(lanes, 1) if i not in ladders]
    samples.sort(key=lambda t: t[1]["x"])
    n = len(samples)
    bands = [[(b["x_img"], b["y_img"] - tilt * b["x_img"]) for b in L["bands"]]
             for _, L in samples]
    owner = [[None] * len(v) for v in bands]
    rows = []

    def original(k, j):
        return (bands[k][j][0], bands[k][j][1] + tilt * bands[k][j][0])

    # start from the lanes with most bands
    for ref in sorted(range(n), key=lambda k: -len(bands[k])):
        for bi in range(len(bands[ref])):
            if owner[ref][bi] is not None:
                continue
            track = {ref: bi}
            for step in (1, -1):
                found = [bands[ref][bi]]
                k = ref + step
                while 0 <= k < n:
                    x_k = samples[k][1]["x"]
                    # prediction: median of the last 3 positions found. Real
                    # drift between neighbouring lanes is < 1 px; extrapolating
                    # a slope would amplify marking noise (1-2 px).
                    last = sorted(q[1] for q in found[-3:])
                    pred = last[len(last) // 2]
                    dx = abs(x_k - found[-1][0])
                    tol_k = tol_px + (min(1.0, 0.004 * dx) if corrected else 0.01 * dx)
                    cand = sorted((abs(y - pred), j) for j, (x, y) in enumerate(bands[k])
                                  if owner[k][j] is None)
                    # accept the nearest band unless the runner-up is almost
                    # as close (ambiguous)
                    if cand and cand[0][0] <= tol_k and not (
                            len(cand) > 1 and cand[1][0] - cand[0][0] < 1.0):
                        j = cand[0][1]
                        track[k] = j
                        found.append(bands[k][j])
                    k += step
            if len(track) >= max(3, min_frac * n):
                for k, j in track.items():
                    owner[k][j] = len(rows)
                row = BandRow([original(k, j) for k, j in track.items()])
                row.lanes = sorted(samples[k][0] for k in track)
                rows.append(row)

    # gap filling (tilt-corrected pass only): extend each accepted row at its
    # ends, one lane at a time, with the free band nearest to the fitted curve
    # (within tol_px + 1), refit, and repeat until nothing changes
    for ri in range(len(rows) if corrected else 0):
        row = rows[ri]
        for _ in range(30):
            new = []
            row_lanes = set(row.lanes)
            ks = [k for k, (li, _) in enumerate(samples) if li in row_lanes]
            kmin, kmax = min(ks), max(ks)
            for k, (li, L) in enumerate(samples):
                if li in row_lanes or not (k == kmin - 1 or k == kmax + 1):
                    continue
                cand = [(abs((y + tilt * x) - row.y(x)), j)
                        for j, (x, y) in enumerate(bands[k]) if owner[k][j] is None]
                if cand:
                    dmin, j = min(cand)
                    if dmin <= tol_px + 1.0:
                        owner[k][j] = ri
                        new.append((k, j))
            if not new:
                break
            pts = row.pts + [original(k, j) for k, j in new]
            lanes_new = sorted(row_lanes | {samples[k][0] for k, _ in new})
            row = BandRow(pts)
            row.lanes = lanes_new
            rows[ri] = row

    # drop rows that cross, or come within 2 px of, a longer accepted row
    xs_test = [L["x"] for L in lanes]
    rows.sort(key=lambda r: -len(r.pts))
    accepted = []
    for r in rows:
        if all(min(abs(r.y(x) - g.y(x)) for x in xs_test) >= 2.0 and
               (r.y(xs_test[0]) - g.y(xs_test[0])) *
               (r.y(xs_test[-1]) - g.y(xs_test[-1])) > 0 for g in accepted):
            accepted.append(r)
    return accepted


def detect_rows(lanes, ladders, tol_px=3.0, min_frac=0.5,
                tilt_threshold_px=TILT_THRESHOLD_PX):
    """Two passes. Pass 1 finds the rows and measures the overall gel tilt
    (median of the row slopes, weighted by lanes). If the tilt across the gel
    exceeds `tilt_threshold_px`, pass 2 tracks the rows again with that tilt
    removed, so a band from another row that crosses the right height because
    the gel is tilted is not picked up by mistake.

    Returns (rows, tilt) with tilt = dict(px=end-to-end tilt, corrected=bool).
    """
    tilt = dict(px=None, corrected=False)
    rows1 = _track_rows(lanes, ladders, tol_px, min_frac, 0.0)
    if not rows1:
        return rows1, tilt
    slopes = []
    for r in rows1:
        xs = [p[0] for p in r.pts]
        ys = [p[1] for p in r.pts]
        c = _polyfit(xs, ys, 1)
        if c:
            slopes += [c[1]] * len(xs)
    if not slopes:
        return rows1, tilt
    slopes.sort()
    s = slopes[len(slopes) // 2]
    span = max(L["x"] for L in lanes) - min(L["x"] for L in lanes)
    tilt["px"] = round(s * span, 1)
    if abs(s) * span < tilt_threshold_px:        # nearly straight gel
        return rows1, tilt
    tilt["corrected"] = True
    return _track_rows(lanes, ladders, tol_px, min_frac, s, corrected=True), tilt


def automatic_curves(lanes, ladders, gap_curves, image_width,
                     tol_px=3.0, min_frac=0.5, tilt_threshold_px=TILT_THRESHOLD_PX):
    """Build Rf curves: well (Rf 0) and front (Rf 1) from the .gap if present,
    plus one curve per detected band row. The Rf of each row is its position
    on the ladder lane, measured linearly between Rf 0 and Rf 1."""
    ref = lanes[(min(ladders) if ladders else 1) - 1]
    top = [c for c in gap_curves if abs(c["rf"]) < 1e-6]
    bottom = [c for c in gap_curves if abs(c["rf"] - 1) < 1e-6]
    y_min = min(L["y0"] for L in lanes)
    y_max = max(L["y0"] + L["length"] for L in lanes)
    c_top = top[0] if top else dict(rf=0.0, pts=[(0, y_min), (image_width, y_min)])
    c_bottom = bottom[0] if bottom else dict(rf=1.0, pts=[(0, y_max), (image_width, y_max)])
    rows, tilt = detect_rows(lanes, ladders, tol_px, min_frac, tilt_threshold_px)
    xs = sorted(set([0, image_width] + [L["x"] for L in lanes]))
    xl = ref["x"]
    yt, yb = curve_y(c_top, xl), curve_y(c_bottom, xl)
    curves, info = [c_top], []
    for r in sorted(rows, key=lambda r: r.y(xl)):
        rf = (r.y(xl) - yt) / (yb - yt)
        if not (0 < rf < 1):
            continue
        curves.append(dict(rf=rf, pts=[(int(round(x)), r.y(x)) for x in xs]))
        info.append(dict(rf=round(rf, 4), n_lanes=len(r.lanes), lanes=r.lanes,
                         y_at_ladder=round(r.y(xl), 1),
                         curvature_px=round(max(r.y(x) for x in xs[1:-1]) -
                                            min(r.y(x) for x in xs[1:-1]), 1),
                         degree=r.degree, residual_sd_px=round(r.residual_sd(), 2)))
    curves.append(c_bottom)
    for L in lanes:                               # curves must never cross
        ys = [curve_y(c, L["x"]) for c in curves]
        if any(b <= a for a, b in zip(ys, ys[1:])):
            raise RuntimeError("automatic Rf curves cross each other; "
                               "try a smaller --auto-tol-px")
    return curves, info, tilt


# =========================================================================== #
# 7. Band table
# =========================================================================== #
def extract_bands(gap_path, curve_mode="auto", ladder_lane=None, mw_model="interp",
                  auto_tol_px=3.0, auto_min_frac=0.5,
                  tilt_threshold_px=TILT_THRESHOLD_PX):
    """Parse a .gap file and return (rows, meta, image, lanes, curves).

    rows   : list of dicts, one per band (the bands_<gel>.csv table)
    meta   : calibration details (the calibration_<gel>.json content)
    image  : decoded grey-scale image, or None if it could not be decoded
    curves : Rf curves actually used (empty if no correction was applied)
    """
    g = parse_gap(gap_path)
    lanes, curves = g["lanes"], g["curves"]

    if ladder_lane is None:
        ladders = {i for i, L in enumerate(lanes, 1)
                   if any(b["ladder_mw"] is not None for b in L["bands"])}
    else:
        ladders = {ladder_lane}

    for L in lanes:
        for b in L["bands"]:
            b["x_img"], b["y_img"] = image_position(L, b["peak"])

    auto_info, tilt = None, dict(px=None, corrected=False)
    if curve_mode == "auto":
        image_width = struct.unpack_from(">I", g["header"]["image"], 16)[0]  # PNG IHDR
        curves, auto_info, tilt = automatic_curves(
            lanes, ladders, curves, image_width, auto_tol_px, auto_min_frac,
            tilt_threshold_px)
    use_curves = len(curves) >= 2 and curve_mode != "none"
    if curve_mode != "none" and not use_curves:
        print("Warning: fewer than 2 Rf curves in the .gap; using "
              "Rf = peak / lane length.", file=sys.stderr)

    # Rf of each band
    for L in lanes:
        for b in L["bands"]:
            b["rf_linear"] = b["peak"] / L["length"]
            b["rf"] = (rf_from_curves(curves, b["x_img"], b["y_img"])
                       if use_curves else b["rf_linear"])

    # ladder: mean Rf of each size across ladder lanes
    by_mw = {}
    for L in lanes:
        for b in L["bands"]:
            if b["ladder_mw"]:
                by_mw.setdefault(b["ladder_mw"], []).append(b["rf"])
    mw_l = sorted(by_mw)
    rf_l = [sum(by_mw[m]) / len(by_mw[m]) for m in mw_l]
    if len(rf_l) < 3:
        raise RuntimeError(f"only {len(rf_l)} ladder band(s) have a size "
                           "assigned in this .gap; at least 3 are needed")
    if mw_model == "exp":
        cal, params = calibration_exp(rf_l, mw_l)
    else:
        cal, params = calibration_interp(rf_l, mw_l), None
    lo_l, hi_l = min(mw_l), max(mw_l)
    residuals = {str(int(m)): round((cal(r) / m - 1) * 100, 2)
                 for r, m in zip(rf_l, mw_l)}

    # densitometry (skipped if the image cannot be decoded)
    img, bg_sources = None, set()
    try:
        img = decode_png(g["header"]["image"])
    except Exception as e:                               # noqa: BLE001
        print(f"Warning: band intensity not computed ({e}).", file=sys.stderr)

    rows = []
    for li, L in enumerate(lanes, start=1):
        is_ladder = li in ladders
        sample = "LADDER" if is_ladder else f"L{li:02d}"
        dens = {}
        if img is not None and L["bands"]:
            prof = lane_profile(img, L)
            bg, source = lane_background(prof, L)
            bg_sources.add(source)
            net = [max(0.0, p - f) for p, f in zip(prof, bg)]
            vols = []
            for b in L["bands"]:
                i0, i1 = max(0, b["start"]), min(len(net) - 1, b["end"])
                vols.append(sum(net[i0:i1 + 1]) * L["width"])
            tot = sum(vols) or float("nan")
            for bi, (b, v) in enumerate(zip(L["bands"], vols), start=1):
                pk = net[min(len(net) - 1, max(0, b["peak"]))]
                dens[bi] = (round(pk, 2), round(v, 1), round(100 * v / tot, 2))
        for bi, b in enumerate(L["bands"], start=1):
            mw = cal(b["rf"])
            pk, vol, rel = dens.get(bi, (None, None, None))
            # Equivalent position in the corrected frame, as INTEGER pixels:
            # bin_gel.R treats peak_px as pixels (detection floor, locus
            # ceiling, contiguous-pixel cuts). The continuous value is kept in
            # peak_px_corrected; MW is computed from the continuous Rf.
            pc = b["rf"] * L["length"]
            shift = int(round(pc)) - b["peak"]
            rows.append(dict(
                lane=li, sample=sample, band=bi, ladder=int(is_ladder),
                peak_px=int(round(pc)),
                start_px=b["start"] + shift, end_px=b["end"] + shift,
                width_px=b["end"] - b["start"],
                lane_length_px=L["length"],
                rf=round(b["rf"], 6), mw_bp=round(mw, 3),
                ladder_mw_bp=b["ladder_mw"],
                extrapolated=int(not (lo_l * (1 - 1e-6) <= mw <= hi_l * (1 + 1e-6))),
                peak_intensity=pk, raw_volume=vol, rel_volume_pct=rel,
                peak_px_corrected=round(pc, 2),
                peak_px_original=b["peak"], rf_linear=round(b["rf_linear"], 6),
                x_img=round(b["x_img"], 1), y_img=round(b["y_img"], 1),
            ))
    if not rows:
        raise RuntimeError("no bands found in the .gap file")

    samples = [r for r in rows if not r["ladder"]]
    meta = dict(
        script="extract_gel.py",
        gap_file=Path(gap_path).name,
        signature=g["header"]["signature"], image=g["header"]["image_name"],
        n_lanes=len(lanes), n_bands=len(rows), ladder_lanes=sorted(ladders),
        n_samples=len({r["sample"] for r in samples}),
        lane_length_px=lanes[0]["length"],
        rf_curves=[dict(rf=round(c["rf"], 4), n_points=len(c["pts"]),
                        points=[(x, round(y, 2)) for x, y in c["pts"]])
                   for c in curves],
        rf_corrected_by_curves=use_curves,
        curve_source={"auto": "automatic (band rows)", "drawn": "drawn in GelAnalyzer",
                      "none": "not used"}[curve_mode],
        auto_rows=auto_info,
        gel_tilt_px=tilt["px"], tilt_corrected=tilt["corrected"],
        mw_model=("MW(Rf) = a*exp(-b*Rf) + c" if mw_model == "exp"
                  else "piecewise-linear interpolation of log10(MW)"),
        mw_model_parameters=params,
        ladder_residual_pct=residuals,
        ladder_range_bp=[lo_l, hi_l],
        sample_mw_range_bp=([min(r["mw_bp"] for r in samples),
                             max(r["mw_bp"] for r in samples)] if samples else None),
        n_extrapolated_bands=sum(r["extrapolated"] for r in samples),
        intensity=("computed from the embedded image; background: " +
                   ",".join(sorted(bg_sources))) if img else "not computed",
    )
    return rows, meta, img, lanes, (curves if use_curves else [])


def weak_lanes(rows, factor=5.0):
    """Sample lanes whose total band volume is below median / factor."""
    tot = {}
    for r in rows:
        if not r["ladder"] and r["raw_volume"] is not None:
            tot[r["lane"]] = tot.get(r["lane"], 0) + r["raw_volume"]
    if not tot:
        return []
    med = statistics.median(tot.values())
    return sorted(l for l, s in tot.items() if s < med / factor)


# =========================================================================== #
# 8. Check image (PNG writer, no dependencies)
# =========================================================================== #
_FONT = {"0": ["111", "101", "101", "101", "111"], "1": ["010", "110", "010", "010", "111"],
         "2": ["111", "001", "111", "100", "111"], "3": ["111", "001", "111", "001", "111"],
         "4": ["101", "101", "111", "001", "001"], "5": ["111", "100", "111", "001", "111"],
         "6": ["111", "100", "111", "101", "111"], "7": ["111", "001", "010", "010", "010"],
         "8": ["111", "101", "111", "101", "111"], "9": ["111", "101", "111", "001", "111"]}


def write_png(path, w, h, rgb):
    """rgb: list of rows, each a bytearray of 3*w bytes."""
    def chunk(kind, data):
        c = kind + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
    raw = b"".join(b"\x00" + bytes(r) for r in rgb)
    png = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))
    Path(path).write_bytes(png)


def check_image(path, img, lanes, curves, rows, scale=2):
    """Gel crop with Rf curves (cyan), lane axes (grey), sample bands (yellow),
    ladder bands (green) and lane numbers (red)."""
    w, h, pix = img
    ys = [r["y_img"] for r in rows]
    y0 = max(0, int(min(ys)) - 40)
    y1 = min(h - 1, int(max(ys)) + 30)
    vals = sorted(v for row in pix[y0:y1] for v in row[::7])
    top = vals[int(0.995 * (len(vals) - 1))] or 1          # contrast stretch
    W, H = w * scale, (y1 - y0) * scale
    rgb = []
    for y in range(y0, y1):
        row = bytearray()
        for v in pix[y]:
            g = int(255 * min(1.0, v / top) ** 0.7)
            row += bytes((g, g, g)) * scale
        for _ in range(scale):
            rgb.append(bytearray(row))

    def dot(x, y, colour, r=0):
        X, Y = int(round(x * scale)), int(round((y - y0) * scale))
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                if 0 <= X + dx < W and 0 <= Y + dy < H:
                    rgb[Y + dy][3 * (X + dx):3 * (X + dx) + 3] = bytes(colour)

    for c in curves:
        if 0 < c["rf"] < 1:
            for xx in range(w):
                dot(xx, curve_y(c, xx), (0, 200, 255), r=1)
    for L in lanes:
        for yy in range(y0, y1, 6):
            t = (yy - L["y0"]) / L["length"]
            dot(L["x"] + t * (L["x_bottom"] - L["x"]), yy, (70, 70, 70))
    for r in rows:
        dot(r["x_img"], r["y_img"], (80, 255, 80) if r["ladder"] else (255, 220, 0),
            r=scale + 1)
    for i, L in enumerate(lanes, 1):                       # 3x5 pixel font
        s = str(i)
        x_start = L["x"] * scale - len(s) * 4 * 3 // 2
        for k, ch in enumerate(s):
            for row_i, bits in enumerate(_FONT[ch]):
                for col_i, bit in enumerate(bits):
                    if bit == "1":
                        for sy in range(3):
                            for sx in range(3):
                                X = x_start + k * 12 + col_i * 3 + sx
                                Y = 4 + row_i * 3 + sy
                                if 0 <= X < W and 0 <= Y < H:
                                    rgb[Y][3 * X:3 * X + 3] = b"\xff\x50\x50"
    write_png(path, W, H, rgb)


# =========================================================================== #
# 9. Optional validation against a GelAnalyzer "Bands Information" export
# =========================================================================== #
def _num(v):
    v = (v or "").strip()
    if not v:
        return None
    if "," in v and "." in v:
        v = v.replace(".", "").replace(",", ".")
    else:
        v = v.replace(",", ".")
    try:
        return float(v)
    except ValueError:
        return None


def validate(rows, gelanalyzer_csv):
    """Print Rf and MW from this script next to those exported by GelAnalyzer."""
    with open(gelanalyzer_csv, newline="", encoding="utf-8-sig") as fh:
        head = fh.read(4096); fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(head, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel_tab
        ref = list(csv.DictReader(fh, dialect=dialect))
    idx = {(r["lane"], r["band"]): r for r in rows}
    print(f"\n{'lane':>4} {'band':>5} {'Rf GA':>7} {'Rf':>7} {'MW GA':>8} "
          f"{'MW':>8} {'diff %':>8}")
    diffs = []
    for r in ref:
        k = {c.strip().lower(): v for c, v in r.items() if c}

        def col(*keys):
            for c, v in k.items():
                if any(key in c for key in keys):
                    return v
            return None
        lane = _num(col("lane"))
        band = _num(col("band #", "band no", "band"))
        rf = _num(col("rf"))
        mw = _num(col("mw", "weight", "size"))
        if lane is None or band is None:
            continue
        m = idx.get((int(lane), int(band)))
        if not m:
            continue
        d = (m["mw_bp"] / mw - 1) * 100 if mw else float("nan")
        diffs.append(abs(d))
        print(f"{int(lane):>4} {int(band):>5} {rf if rf is not None else float('nan'):>7.3f} "
              f"{m['rf']:>7.3f} {mw or float('nan'):>8.1f} {m['mw_bp']:>8.1f} {d:>8.2f}")
    if diffs:
        print(f"\n{len(diffs)} bands compared; mean MW difference "
              f"{sum(diffs) / len(diffs):.2f} %, max {max(diffs):.2f} %")


# =========================================================================== #
# 10. Command line
# =========================================================================== #
def run(gap, curve_mode="auto", out_dir=None, ladder_lane=None, mw_model="interp",
        auto_tol_px=3.0, auto_min_frac=0.5, tilt_threshold_px=TILT_THRESHOLD_PX):
    gap = Path(gap)
    if not gap.exists():
        raise RuntimeError(f"file not found: {gap}")
    folder = Path(out_dir) if out_dir else gap.parent
    folder.mkdir(parents=True, exist_ok=True)
    name = gap.stem
    out_csv = folder / f"bands_{name}.csv"
    out_json = folder / f"calibration_{name}.json"
    out_png = folder / f"check_{name}.png"
    print(f"\nGel: {name}")

    rows, meta, img, lanes, curves = extract_bands(
        gap, curve_mode, ladder_lane, mw_model, auto_tol_px, auto_min_frac,
        tilt_threshold_px)

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    out_json.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    if img is not None:
        check_image(out_png, img, lanes, curves, rows)

    samples = [r for r in rows if not r["ladder"]]
    print(f"  lanes: {meta['n_lanes']} ({len(meta['ladder_lanes'])} ladder, "
          f"{len({r['lane'] for r in samples})} samples) | sample bands: {len(samples)}")
    if curve_mode == "auto":
        f = meta["auto_rows"] or []
        print(f"  Rf curves: automatic, {len(f)} band rows "
              f"(present in {', '.join(str(x['n_lanes']) for x in f)} lanes)")
        if len(f) < 2:
            print("  WARNING: few band rows found; the distortion correction is weak on this gel.")
        if meta["gel_tilt_px"] is not None:
            print(f"  gel tilt: {meta['gel_tilt_px']:+.1f} px end to end"
                  + (" -> rows re-tracked with the tilt removed"
                     if meta["tilt_corrected"] else ""))
    else:
        print(f"  Rf curves: {curve_mode}")
    print(f"  MW model: {meta['mw_model']}")
    print(f"  ladder range: {meta['ladder_range_bp'][0]:g}-{meta['ladder_range_bp'][1]:g} bp"
          f" | bands outside it: {meta['n_extrapolated_bands']}")
    weak = weak_lanes(rows)
    if weak:
        print(f"  WARNING: lanes with very weak signal (volume < 1/5 of the median): {weak}")
        print("    Missing bands in these lanes may be PCR failures. To exclude them,")
        print("    use in step 2:  --exclude " + ",".join(map(str, weak)))
    print(f"\n  saved: {out_csv.name}, {out_png.name if img else '(no image)'}, {out_json.name}")
    if img is not None:
        print(f"  CHECK {out_png.name}: each cyan curve should follow a row of bands.")
    print(f"\nNext step:\n  Rscript bin_gel.R \"{out_csv}\"")
    return out_csv, rows


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="STEP 1: extract the band table from one GelAnalyzer .gap file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  python3 extract_gel.py "Primer i5-Population2.gap"')
    ap.add_argument("gap", help="GelAnalyzer .gap file")
    ap.add_argument("--curves", choices=["auto", "drawn", "none"], default="auto",
                    help="gel distortion correction: auto (default), drawn or none")
    ap.add_argument("--out", default=None, help="output folder (default: folder of the .gap)")
    ap.add_argument("--ladder-lane", type=int, default=None,
                    help="ladder lane number (default: every lane with sizes assigned)")
    ap.add_argument("--mw-model", choices=["interp", "exp"], default="interp",
                    help="size calibration: piecewise log-linear (default) or a*exp(-b*Rf)+c")
    ap.add_argument("--auto-tol-px", type=float, default=3.0,
                    help="max vertical step (px) when tracking band rows (default 3)")
    ap.add_argument("--auto-min-frac", type=float, default=0.5,
                    help="min fraction of sample lanes a row must span (default 0.5)")
    ap.add_argument("--no-tilt-correction", action="store_true",
                    help="do not remove the overall gel tilt when tracking rows")
    ap.add_argument("--validate", default=None, metavar="CSV",
                    help="GelAnalyzer 'Bands Information' export to compare against")
    a = ap.parse_args(argv)
    try:
        _, rows = run(a.gap, a.curves, a.out, a.ladder_lane, a.mw_model,
                      a.auto_tol_px, a.auto_min_frac,
                      float("inf") if a.no_tilt_correction else TILT_THRESHOLD_PX)
    except RuntimeError as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1
    if a.validate:
        validate(rows, a.validate)
    return 0


if __name__ == "__main__":
    sys.exit(main())
