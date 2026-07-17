#!/usr/bin/env python3
"""
Plot gas density movie frames written by RAMSES movie.f90,
overplotting the sink particle (black hole) position.

Map file layout (Fortran sequential unformatted, written by movie.f90):
    record 1 : t (or aexp if tendmov<=0), delx, dely, delz    -- 4 x real(8)
    record 2 : nw_frame, nh_frame                             -- 2 x int(4)
    record 3 : real(data_frame(:,:,imap),4)                   -- nw*nh x real(4)

data_frame(ii,jj) has ii along the horizontal (x) axis and jj along the
vertical (y) axis of the projection, stored column-major.

The map header does NOT store the camera centre, so it is reconstructed with
the same logic as movie.f90: if do_center_track_movie is on, aexp exceeds
aexp_track_movie and the tracker is active, the centre is the tracker centre
of mass (read from tracker_com_*.txt); otherwise the namelist centre is used.

Usage:
    python3 plot_movie_dens.py
    python3 plot_movie_dens.py --frames 20 21 --width 2.0 --outdir frames
"""

import argparse
import glob
import os
import re

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------- constants
MSUN_G = 1.98892e33
MH_G = 1.6600000e-24
KPC_CM = 3.08567758128e21
X_H = 0.76  # hydrogen mass fraction, for the n_H conversion

# Camera settings, filled by read_namelist() from the run's &MOVIE_PARAMS block.
# The frame centre is not stored in the .map header, so it comes from the
# namelist, overridden per frame by the tracker centre (see frame_centre).
NML = {}

DEFAULT_INDIR = "./movie1"

# Where to hunt for the run namelist, in order, relative to --indir. A run can
# hold several .nml files side by side (e.g. ramses.nml next to RTMHD.nml), so
# candidates are chosen by containing a &MOVIE_PARAMS block, not by filename.
NML_SEARCH_DIRS = ("../ics_ramses", "..")



# ------------------------------------------------------------------ readers
def has_movie_params(path):
    """True if this file contains an uncommented &MOVIE_PARAMS group."""
    try:
        with open(path, errors="replace") as f:
            return re.search(r"^\s*&MOVIE_PARAMS\b", f.read(), re.M | re.I) is not None
    except OSError:
        return False


def find_namelist(indir):
    """Locate the run namelist by searching NML_SEARCH_DIRS relative to indir.

    Selection is by content (must hold a &MOVIE_PARAMS group) rather than by
    name, since e.g. ics_ramses/ may contain both ramses.nml and RTMHD.nml and
    only one of them configures the movie.
    """
    seen, candidates = set(), []
    for d in NML_SEARCH_DIRS:
        for p in sorted(glob.glob(os.path.join(indir, d, "*.nml"))):
            real = os.path.realpath(p)
            if real not in seen:
                seen.add(real)
                candidates.append(p)

    matches = [p for p in candidates if has_movie_params(p)]
    if not matches:
        looked = ", ".join(os.path.normpath(os.path.join(indir, d)) for d in NML_SEARCH_DIRS)
        found = ", ".join(os.path.basename(p) for p in candidates) or "none"
        raise SystemExit(
            f"No namelist with a &MOVIE_PARAMS block found.\n"
            f"  looked in : {looked}\n"
            f"  .nml found: {found}\n"
            f"Pass one explicitly with --nml."
        )
    if len(matches) > 1:
        print(f"  warning: several namelists define &MOVIE_PARAMS "
              f"({', '.join(matches)}); using the first. Override with --nml.")
    return matches[0]


def read_namelist(path):
    """Parse &MOVIE_PARAMS from the run's namelist into the NML dict.

    Only the entries needed to reconstruct the camera centre are pulled out.
    Fortran namelist values may be comma-separated lists (e.g. the four
    polynomial-in-aexp coefficients of xcentre_frame); movie.f90 evaluates that
    polynomial, but the higher coefficients are zero for this run, so only the
    constant term is used here.
    """
    with open(path) as f:
        txt = f.read()

    m = re.search(r"^\s*&MOVIE_PARAMS\b(.*?)^\s*/\s*$", txt, re.S | re.M | re.I)
    if not m:
        raise SystemExit(f"{path}: no &MOVIE_PARAMS block found")

    raw = {}
    for line in m.group(1).splitlines():
        line = line.split("!")[0].strip()
        if not line or "=" not in line:
            continue
        key, _, val = line.partition("=")
        raw[key.strip().lower()] = val.strip()

    def first(key, default=""):
        return raw.get(key, default).split(",")[0].strip()

    def flt(key, default="0"):
        v = first(key, default)
        return float(v.replace("D", "E").replace("d", "e")) if v else float(default)

    axis = first("proj_axis", "'z'").strip("'\" ").lower()
    if len(axis) > 1:
        # movie.f90 loops over each character, writing movie1/, movie2/, ...
        print(f"  note: proj_axis='{axis}', using '{axis[0]}' for movie1")
    axis = axis[0] if axis else "z"

    cfg = dict(
        proj_axis=axis,
        xcentre_frame=flt("xcentre_frame"),
        ycentre_frame=flt("ycentre_frame"),
        zcentre_frame=flt("zcentre_frame"),
        do_center_track_movie=first("do_center_track_movie", ".false.").lower().startswith(".t"),
        aexp_track_movie=flt("aexp_track_movie"),
    )
    NML.clear()
    NML.update(cfg)
    return cfg


def read_map(path):
    """Read one RAMSES movie .map file. Returns (time_or_aexp, del, data[nw,nh])."""
    with open(path, "rb") as f:
        raw = f.read()

    off = 0

    def record():
        """Return the payload of one Fortran sequential record."""
        nonlocal off
        n = np.frombuffer(raw, dtype=np.int32, count=1, offset=off)[0]
        payload = raw[off + 4: off + 4 + n]
        tail = np.frombuffer(raw, dtype=np.int32, count=1, offset=off + 4 + n)[0]
        if tail != n:
            raise IOError(f"{path}: corrupt record (head {n} != tail {tail})")
        off += 8 + n
        return payload

    head = np.frombuffer(record(), dtype=np.float64)
    tval, delx, dely, delz = head[0], head[1], head[2], head[3]

    dims = np.frombuffer(record(), dtype=np.int32)
    nw, nh = int(dims[0]), int(dims[1])

    data = np.frombuffer(record(), dtype=np.float32)
    if data.size != nw * nh:
        raise IOError(f"{path}: expected {nw * nh} pixels, got {data.size}")
    # Fortran column-major: data[ii, jj] with ii -> x, jj -> y
    data = data.reshape((nw, nh), order="F").astype(np.float64)

    return dict(t=tval, delx=delx, dely=dely, delz=delz, nw=nw, nh=nh, data=data)


def read_info(path):
    """Parse a movie info_*.txt into a dict of floats/ints where possible."""
    info = {}
    with open(path) as f:
        for line in f:
            if "=" not in line or line.lstrip().startswith("DOMAIN"):
                continue
            key, _, val = line.partition("=")
            key, val = key.strip(), val.strip()
            if not key:
                continue
            try:
                info[key] = float(val) if ("." in val or "E" in val.upper()) else int(val)
            except ValueError:
                info[key] = val
    return info


def read_tracker(path):
    """Read tracker_com_*.txt: returns (active, xpos) for the first tracker."""
    if not os.path.exists(path):
        return False, None
    txt = open(path).read()
    # Only look at the block for itrack = 1, before the "Tracked particle IDs" section
    txt = txt.split("Tracked particle IDs")[0]

    def scalar(name):
        m = re.search(rf"^{name}\s*=\s*([-\d.EDed+]+)", txt, re.M)
        return float(m.group(1).replace("D", "E").replace("d", "e")) if m else None

    active = scalar("active")
    xyz = [scalar("xtrack"), scalar("ytrack"), scalar("ztrack")]
    if active is None or any(v is None for v in xyz):
        return False, None
    return bool(int(active)), np.array(xyz)


def read_sinks(path):
    """Read sink_*.txt (written by output_sink_csv). Returns a list of dicts.

    Columns: id, m, x, y, z, vx, vy, vz, lx, ly, lz, age(t-tsink), dMBHoverdt
    """
    sinks = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return sinks
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 13:
                continue
            v = [float(p) for p in parts]
            sinks.append(dict(
                id=int(v[0]), m=v[1], x=v[2], y=v[3], z=v[4],
                v=np.array(v[5:8]), l=np.array(v[8:11]),
                age=v[11], dmdt=v[12],
            ))
    return sinks


# ------------------------------------------------------------------ helpers
def frame_centre(aexp, tracker_active, tracker_pos):
    """Reproduce the camera centre logic of movie.f90 (see its lines 296-320).

    The tracker centre is only used once the tracker is active AND aexp has
    passed aexp_track_movie -- before that, tracker_com_*.txt just holds the
    placeholder 0.5,0.5,0.5, so the namelist centre is the real camera centre.

    Returns (centre, source) where source is 'tracker' or 'namelist'.
    """
    if (NML["do_center_track_movie"] and aexp > NML["aexp_track_movie"]
            and tracker_active and tracker_pos is not None):
        cen, source = tracker_pos.copy(), "tracker"
    else:
        cen = np.array([NML["xcentre_frame"], NML["ycentre_frame"], NML["zcentre_frame"]])
        source = "namelist"
    cen[cen == 0.0] = 0.5  # movie.f90 falls back to boxlen/2
    return cen, source


def projected_axes(proj_axis):
    """Return (horizontal, vertical, line-of-sight) coordinate indices."""
    return {"x": (1, 2, 0), "y": (0, 2, 1), "z": (0, 1, 2)}[proj_axis]


def find_frames(indir):
    nums = []
    for p in glob.glob(os.path.join(indir, "dens_*.map")):
        m = re.search(r"dens_(\d+)\.map$", p)
        if m:
            nums.append(int(m.group(1)))
    return sorted(nums)


# ------------------------------------------------------------------- plotting
def plot_frame(num, indir, outdir, args):
    tag = f"{num:05d}"
    mp = read_map(os.path.join(indir, f"dens_{tag}.map"))
    info = read_info(os.path.join(indir, f"info_{tag}.txt"))
    tracker_active, tracker_pos = read_tracker(os.path.join(indir, f"tracker_com_{tag}.txt"))
    sinks = read_sinks(os.path.join(indir, f"sink_{tag}.txt"))

    aexp = info["aexp"]
    redshift = 1.0 / aexp - 1.0
    unit_d, unit_l = info["unit_d"], info["unit_l"]
    unit_m = unit_d * unit_l ** 3

    cen, cen_source = frame_centre(aexp, tracker_active, tracker_pos)
    ih, iv, ilos = projected_axes(NML["proj_axis"])
    hcen, vcen = cen[ih], cen[iv]

    # Physical (proper) kpc per code length unit at this epoch
    kpc_per_code = unit_l / KPC_CM

    # Map is mass-weighted mean density along the line of sight (method_frame
    # defaults to 'mean_mass'), in code units -> hydrogen number density.
    nH = mp["data"] * unit_d * X_H / MH_G

    # Rows must be the vertical axis for imshow, so transpose data[x, y].
    img = nH.T

    half_h, half_v = mp["delx"] / 2.0, mp["dely"] / 2.0
    extent = np.array([-half_h, half_h, -half_v, half_v]) * kpc_per_code

    if args.width is not None:  # optional zoom, in proper kpc
        lim = args.width / 2.0
        extent_lim = [-lim, lim, -lim, lim]
    else:
        extent_lim = [extent[0], extent[1], extent[2], extent[3]]

    pos = img[img > 0]
    if pos.size == 0:
        print(f"  frame {tag}: map is empty, skipping")
        return
    vmin = args.vmin if args.vmin is not None else np.percentile(pos, 1.0)
    vmax = args.vmax if args.vmax is not None else np.percentile(pos, 99.99)
    if vmin <= 0 or vmin >= vmax:
        vmin = max(pos.min(), vmax * 1e-6)

    fig, ax = plt.subplots(figsize=(8.2, 8.0))
    im = ax.imshow(
        np.clip(img, vmin, vmax),
        origin="lower", extent=extent, norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=args.cmap, interpolation="nearest",
    )

    # Sink particles: convert code position -> proper kpc offset from centre
    bh_label = ""
    for s in sinks:
        sh = (s[["x", "y", "z"][ih]] - hcen) * kpc_per_code
        sv = (s[["x", "y", "z"][iv]] - vcen) * kpc_per_code
        mass_msun = s["m"] * unit_m / MSUN_G
        bh_label += f"   BH : {mass_msun:.2e} " + r"M$_\odot$"
        ax.plot(sh, sv, marker="o", ms=6, mfc="black", mec="orange", mew=1.8, zorder=5)
        #ax.plot(sh, sv, marker="+", ms=8, color="white", mew=1.4, zorder=6)
        #ax.annotate(
        #    f"BH : {mass_msun:.2e} " + r"M$_\odot$",
        #    xy=(sh, sv), xytext=(10, 10), textcoords="offset points",
        #    color="white", fontsize=9, zorder=6,
        #)

    ax.set_xlim(extent_lim[0], extent_lim[1])
    ax.set_ylim(extent_lim[2], extent_lim[3])
    hl, vl = "xyz"[ih], "xyz"[iv]
    ax.set_xlabel(f"{hl} [proper kpc]")
    ax.set_ylabel(f"{vl} [proper kpc]")
    ax.set_title(f"frame {tag}   z = {redshift:.2f}   a = {aexp:.4f}" + bh_label)

    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(r"mass-weighted $n_{\rm H}$ [cm$^{-3}$]")

    out = os.path.join(outdir, f"dens_{tag}.png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  frame {tag}: z={redshift:5.2f}  centre={cen_source:8s}  "
          f"sinks={len(sinks)}  -> {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indir", default=DEFAULT_INDIR, help="movie directory")
    p.add_argument("--nml", default=None,
                   help="run namelist holding &MOVIE_PARAMS (default: search "
                        + " and ".join(NML_SEARCH_DIRS) + " relative to --indir)")
    p.add_argument("--outdir", default="dens_images", help="output directory for PNGs")
    p.add_argument("--frames", type=int, nargs="+", default=None,
                   help="frame numbers to plot (default: all)")
    p.add_argument("--width", type=float, default=None,
                   help="zoom width in proper kpc (default: full frame)")
    p.add_argument("--vmin", type=float, default=None, help="colour scale min [cm^-3]")
    p.add_argument("--vmax", type=float, default=None, help="colour scale max [cm^-3]")
    p.add_argument("--cmap", default="inferno")
    p.add_argument("--dpi", type=int, default=140)
    args = p.parse_args()

    if args.nml:
        nml_path = args.nml
        if not os.path.exists(nml_path):
            raise SystemExit(f"Namelist not found: {nml_path}")
    else:
        nml_path = find_namelist(args.indir)
    cfg = read_namelist(nml_path)
    print(f"Camera from {nml_path}:")
    print(f"  proj_axis={cfg['proj_axis']}  centre=("
          f"{cfg['xcentre_frame']:.9f}, {cfg['ycentre_frame']:.9f}, {cfg['zcentre_frame']:.9f})")
    print(f"  do_center_track_movie={cfg['do_center_track_movie']}  "
          f"aexp_track_movie={cfg['aexp_track_movie']:.6g}")

    os.makedirs(args.outdir, exist_ok=True)
    frames = args.frames if args.frames is not None else find_frames(args.indir)
    if not frames:
        raise SystemExit(f"No dens_*.map files found in {args.indir}")

    print(f"Plotting {len(frames)} frame(s) from {args.indir}")
    for num in frames:
        plot_frame(num, args.indir, args.outdir, args)


if __name__ == "__main__":
    main()
