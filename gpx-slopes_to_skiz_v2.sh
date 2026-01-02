#!/usr/bin/env bash
# gpx_to_skiz (release)
# Convert GPX (.gpx) and Slopes (.slopes) to Ski Tracks .skiz.
#
# Why this exists:
# - Ski Tracks can import .skiz (a zip with specific files).
# - GPX alone doesn't contain Ski Tracks' proprietary sensor streams, but we can
#   generate a structurally compatible SKIZ that survives "Refresh segments".
#
# Outputs:
#   ./gpx_to_skiz/<base>.skiz
# Logs:
#   ./gpx_to_skiz/logs/<base>_*.log
#
# Supported inputs:
#   - .gpx     (GPX track)
#   - .slopes  (Slopes app export: contains GPS.csv, RawGPS.csv, Metadata.xml)
#
# Exit codes:
#   0 OK, 1 partial/failed conversions
set -euo pipefail
IFS=$'\n\t'

SCRIPT_NAME="$(basename "$0")"
OUT_DIR="gpx_to_skiz"
LOG_DIR="$OUT_DIR/logs"
VERBOSE=1

# Default tz written in Track.xml for GPX inputs. Slopes inputs default to their Metadata timeZoneOffset.
TZ_OFFSET="+01:00"
MODE="v7.1"  # v7 (conservative) or v7.1 (sensitive) for GPX segmentation only

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [OPTIONS] [FILES...]

Converts GPX (.gpx) and Slopes (.slopes) to Ski Tracks .skiz.

Options:
  -o, --out-dir DIR        Output directory (default: $OUT_DIR)
  --tz +HH:MM|-HH:MM       Timezone offset for GPX Track.xml (default: $TZ_OFFSET)
  --mode v7|v7.1           GPX segmentation mode (default: $MODE)
  -q, --quiet              Less console output
  -h, --help               Show help

If no FILES are provided, converts all *.gpx and *.slopes in the current directory.

Examples:
  ./$SCRIPT_NAME
  ./$SCRIPT_NAME activity.gpx
  ./$SCRIPT_NAME "1 janvier 2026 - Grand Massif.slopes"
  ./$SCRIPT_NAME --mode v7 --tz +01:00 *.gpx *.slopes
EOF
}

mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/${SCRIPT_NAME%.sh}_$(date -u '+%Y%m%dT%H%M%SZ').log"
: >"$MAIN_LOG"

log() {
  local level="$1"; shift
  local ts
  ts="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  if [[ "$VERBOSE" -eq 1 ]]; then
    printf '[%s] %-5s %s\n' "$ts" "$level" "$*" | tee -a "$MAIN_LOG" >&2
  else
    printf '[%s] %-5s %s\n' "$ts" "$level" "$*" >>"$MAIN_LOG"
  fi
}

die() {
  log ERROR "$*"
  exit 1
}

need_cmd() { command -v "$1" >/dev/null 2>&1; }

as_root() {
  if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
    "$@"
  elif need_cmd sudo; then
    sudo "$@"
  else
    die "Need root privileges to install dependencies (run as root or install sudo)."
  fi
}

install_deps() {
  local missing=()
  for c in python3 zip unzip; do
    need_cmd "$c" || missing+=("$c")
  done

  if [[ "${#missing[@]}" -eq 0 ]]; then
    log INFO "Dependencies already installed."
    return 0
  fi

  log INFO "Installing missing dependencies: ${missing[*]}"

  if need_cmd apk; then
    log INFO "Detected Alpine (apk)."
    as_root apk add --no-cache python3 zip unzip ca-certificates >/dev/null
  elif need_cmd apt-get; then
    log INFO "Detected Debian/Ubuntu (apt)."
    as_root apt-get update -y >/dev/null
    as_root apt-get install -y python3 zip unzip ca-certificates >/dev/null
  else
    die "Unsupported distro (no apk/apt-get). Please install python3, zip, unzip."
  fi

  for c in python3 zip unzip; do
    need_cmd "$c" || die "Failed to install dependency: $c"
  done
}

# -----------------
# Argument parsing
# -----------------
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out-dir)
      shift || true
      [[ $# -gt 0 ]] || die "Missing value for --out-dir"
      OUT_DIR="$1"
      LOG_DIR="$OUT_DIR/logs"
      mkdir -p "$LOG_DIR"
      ;;
    --tz)
      shift || true
      [[ $# -gt 0 ]] || die "Missing value for --tz"
      TZ_OFFSET="$1"
      ;;
    --mode)
      shift || true
      [[ $# -gt 0 ]] || die "Missing value for --mode"
      MODE="$1"
      ;;
    -q|--quiet)
      VERBOSE=0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift || true
      while [[ $# -gt 0 ]]; do ARGS+=("$1"); shift || true; done
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      ARGS+=("$1")
      ;;
  esac
  shift || true
done

case "$MODE" in
  v7|v7.1) ;;
  *) die "Invalid --mode: $MODE (expected v7 or v7.1)" ;;
esac

log INFO "Starting $SCRIPT_NAME"
log INFO "Output directory: $OUT_DIR"
log INFO "Mode: $MODE"
log INFO "TZ (GPX): $TZ_OFFSET"
log INFO "Main log: $MAIN_LOG"

install_deps

# Collect inputs
FILES=()
if [[ "${#ARGS[@]}" -eq 0 ]]; then
  while IFS= read -r f; do FILES+=("$f"); done < <(find . -maxdepth 1 -type f \( -name '*.gpx' -o -name '*.slopes' \) -print | sort)
else
  FILES=("${ARGS[@]}")
fi

[[ "${#FILES[@]}" -gt 0 ]] || die "No .gpx or .slopes files found."

convert_one() {
  local inpath="$1"
  [[ -f "$inpath" ]] || { log ERROR "File not found: $inpath"; return 1; }

  local base per_log out_skiz
  base="$(basename "$inpath")"
  base="${base%.*}"
  out_skiz="$OUT_DIR/$base.skiz"
  per_log="$LOG_DIR/${base}_$(date -u '+%Y%m%dT%H%M%SZ').log"
  : >"$per_log"
  mkdir -p "$OUT_DIR"

  log INFO "Converting: $inpath -> $out_skiz"

  if ! python3 - "$inpath" "$out_skiz" "$TZ_OFFSET" "$MODE" >>"$per_log" 2>&1 <<'PY'
import sys, os, math, uuid, random, re, zipfile, tempfile, shutil, datetime
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

INPATH   = sys.argv[1]
OUT_SKIZ = sys.argv[2]
TZ_GPX   = sys.argv[3]   # used for GPX inputs
MODE     = sys.argv[4]   # v7 or v7.1 (GPX segmentation thresholds only)

def iso_z(ts: float) -> str:
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc)
    return dt.isoformat(timespec="milliseconds").replace("+00:00", "Z")

def tz_minutes(s: str) -> int:
    m = re.match(r"^([+-])(\d\d):(\d\d)$", s)
    if not m:
        return 0
    sign = 1 if m.group(1) == "+" else -1
    return sign * (int(m.group(2)) * 60 + int(m.group(3)))

def iso_with_tz(ts: float, tz_str: str) -> str:
    tzinfo = datetime.timezone(datetime.timedelta(minutes=tz_minutes(tz_str)))
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc).astimezone(tzinfo)
    return dt.isoformat(timespec="milliseconds")

def parse_tz_slopes(s: str) -> str:
    # Slopes: "+0100" / "-0500"
    m = re.match(r"^([+-])(\d\d)(\d\d)$", s.strip())
    if not m:
        return TZ_GPX
    return f"{m.group(1)}{m.group(2)}:{m.group(3)}"

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000.0
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1); dl = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dl/2)**2
    return 2 * R * math.asin(math.sqrt(a))

def bearing(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360

def pressure_kpa_from_alt_m(alt_m, p0=101.325):
    # ISA troposphere approximation
    T0=288.15; L=0.0065; g=9.80665; M=0.0289644; R=8.3144598
    h = max(-500.0, float(alt_m))
    factor = 1.0 - (L*h)/T0
    if factor <= 0:
        factor = 1e-6
    exp = (g*M)/(R*L)
    return float(p0 * (factor ** exp))

# ---------- GPX parsing ----------
def parse_gpx(path):
    raw = open(path, "rb").read()
    root = ET.fromstring(raw)
    ns = {}
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}")[0][1:]
        ns = {"g": ns_uri}

    trk = root.find("g:trk", ns) if ns else root.find("trk")
    if trk is None:
        raise SystemExit("No <trk> found in GPX")

    name = (trk.findtext("g:name", default="", namespaces=ns) if ns else trk.findtext("name", default="")) or os.path.splitext(os.path.basename(path))[0]
    desc = (trk.findtext("g:desc", default="", namespaces=ns) if ns else trk.findtext("desc", default="")) or ""

    pts = []
    trksegs = trk.findall("g:trkseg", ns) if ns else trk.findall("trkseg")
    for seg in trksegs:
        trkpts = seg.findall("g:trkpt", ns) if ns else seg.findall("trkpt")
        for pt in trkpts:
            lat = float(pt.attrib.get("lat"))
            lon = float(pt.attrib.get("lon"))
            ele_txt = (pt.findtext("g:ele", default=None, namespaces=ns) if ns else pt.findtext("ele", default=None))
            time_txt = (pt.findtext("g:time", default=None, namespaces=ns) if ns else pt.findtext("time", default=None))
            if time_txt is None:
                continue
            dt = datetime.datetime.fromisoformat(time_txt.replace("Z", "+00:00"))
            ts = dt.timestamp()
            ele = float(ele_txt) if ele_txt is not None else float("nan")
            pts.append((ts, lat, lon, ele))

    pts.sort(key=lambda x: x[0])
    if len(pts) < 2:
        raise SystemExit(f"Not enough timed points in GPX (found {len(pts)})")

    # Fill missing elevation by interpolation
    times = [p[0] for p in pts]
    eles = [p[3] for p in pts]
    if any(math.isnan(e) for e in eles):
        ok = [(t, e) for t, e in zip(times, eles) if not math.isnan(e)]
        if len(ok) >= 2:
            ok_t = [t for t, _ in ok]
            ok_e = [e for _, e in ok]
            def interp(t):
                if t <= ok_t[0]: return ok_e[0]
                if t >= ok_t[-1]: return ok_e[-1]
                import bisect
                j = bisect.bisect_right(ok_t, t) - 1
                t0, t1 = ok_t[j], ok_t[j+1]
                e0, e1 = ok_e[j], ok_e[j+1]
                w = (t - t0) / (t1 - t0)
                return e0*(1-w) + e1*w
            pts = [(t, la, lo, (interp(t) if math.isnan(e) else e)) for (t, la, lo, e) in pts]
        else:
            pts = [(t, la, lo, (0.0 if math.isnan(e) else e)) for (t, la, lo, e) in pts]

    return name, desc, pts

# ---------- Slopes parsing ----------
def read_slopes(zip_path):
    with zipfile.ZipFile(zip_path, "r") as z:
        gps = z.read("GPS.csv").decode("utf-8", errors="replace").splitlines()
        raw = z.read("RawGPS.csv").decode("utf-8", errors="replace").splitlines()
        meta = ET.fromstring(z.read("Metadata.xml"))

    # GPS.csv columns (8):
    # time,lat,lon,alt,course,speed,x,y (x,y are accuracies in Slopes)
    pts = []
    for ln in gps:
        if not ln.strip():
            continue
        c = ln.split(",")
        if len(c) < 6:
            continue
        ts = float(c[0])
        pts.append((ts, float(c[1]), float(c[2]), float(c[3])))

    if len(pts) < 2:
        raise SystemExit("Slopes GPS.csv has too few points")

    # Use Metadata for name/desc/tz + actions for segments
    name = meta.attrib.get("locationName") or os.path.splitext(os.path.basename(zip_path))[0]
    desc = f"{meta.attrib.get('locationName','')}".strip()
    tz = parse_tz_slopes(meta.attrib.get("timeZoneOffset", "").strip() or "+0000")

    actions_el = meta.find("actions")
    actions = []
    if actions_el is not None:
        for a in actions_el.findall("Action"):
            at = a.attrib
            typ = at.get("type", "")
            start = at.get("start", "")
            end = at.get("end", "")
            # "2026-01-01 11:31:18 +0100"
            def parse_dt(s):
                m = re.match(r"^(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2}) ([+-]\d{4})$", s.strip())
                if not m:
                    return None
                # parse with offset
                off = parse_tz_slopes(m.group(3))
                tzinfo = datetime.timezone(datetime.timedelta(minutes=tz_minutes(off)))
                dt = datetime.datetime.strptime(f"{m.group(1)} {m.group(2)}", "%Y-%m-%d %H:%M:%S")
                dt = dt.replace(tzinfo=tzinfo)
                return dt.timestamp()
            st = parse_dt(start)
            en = parse_dt(end)
            if st is None or en is None:
                continue
            actions.append((st, en, typ))

    return name, desc, tz, pts, gps, raw, actions

# ---------- Shared helpers ----------
def build_nodes_repo_style(pts, course, speed):
    # Repo nodes.style expects:
    # time(epoch with .000), lat(6), lon(6), alt(2), heading(2), speed(6), hdop, 0.25
    lines=[]
    for i,(ts,lat,lon,alt) in enumerate(pts):
        hdop = 0.0
        lines.append(f"{ts:.3f},{lat:.6f},{lon:.6f},{alt:.2f},{course[i]:.2f},{speed[i]:.6f},{hdop:.1f},0.25")
    return lines

def build_rawlocations_9cols(pts, course, speed, hacc=46.0, vacc=7.58):
    # 9 columns we used for stability with Ski Tracks:
    # ts,lat,lon,alt,course,speed,hacc,vacc,ts2
    lines=[]
    for i,(ts,lat,lon,alt) in enumerate(pts):
        h = 41.642 if i==0 else hacc
        v = 7.577 if i==0 else vacc
        lines.append(f"{ts:.3f},{lat:.9f},{lon:.9f},{alt:.3f},{course[i]:.2f},{speed[i]:.6f},{h:.3f},{v:.3f},{(ts+0.009):.3f}")
    return lines

def step_metrics(pts):
    n=len(pts)
    dist=[0.0]*n
    dt=[0.0]*n
    speed=[0.0]*n
    course=[0.0]*n
    grade=[0.0]*n
    for i in range(1,n):
        t0,la0,lo0,e0=pts[i-1]
        t1,la1,lo1,e1=pts[i]
        d=haversine(la0,lo0,la1,lo1)
        delta=max(0.001, t1-t0)
        dist[i]=d
        dt[i]=delta
        speed[i]=d/delta if delta>0 else 0.0
        course[i]=bearing(la0,lo0,la1,lo1) if d>0 else course[i-1]
        de=e1-e0
        grade[i]=(abs(de)/d*100.0) if d>0 else 0.0
    # sanitize first point: copy second (avoid -1 sentinels)
    if n >= 2:
        speed[0]=speed[1]
        course[0]=course[1]
    return dist,dt,speed,course,grade

def alt_at_linear(pts, t):
    # pts sorted by ts
    if t <= pts[0][0]: return pts[0][3]
    if t >= pts[-1][0]: return pts[-1][3]
    import bisect
    ts=[p[0] for p in pts]
    j=bisect.bisect_right(ts, t)-1
    t0,t1=ts[j],ts[j+1]
    a0,a1=pts[j][3],pts[j+1][3]
    if t1==t0: return a0
    w=(t-t0)/(t1-t0)
    return a0*(1-w)+a1*w

def build_relative_alt_sensor(pts, interval=1.315, lag=1.650):
    start_ts=pts[0][0]; end_ts=pts[-1][0]
    times=[]
    t=start_ts+1.3
    while t <= end_ts:
        times.append(t); t += interval
    if not times:
        times=[start_ts]
    alt0=alt_at_linear(pts, times[0])
    lines=[]
    for t in times:
        alt=alt_at_linear(pts, t)
        rel=alt-alt0
        p=pressure_kpa_from_alt_m(alt)+random.uniform(-0.02,0.02)
        lines.append(f"{iso_z(t)},{iso_z(max(start_ts,t-lag))},{p:.4f},{rel:.4f}")
    return lines

def segmentize_gpx(pts, dt, mode):
    # classify by vertical rate with hysteresis-ish memory via "flat inherits previous"
    if mode == "v7":
        vr_thr=0.15
        min_seg_s=180.0
    else:
        vr_thr=0.12
        min_seg_s=90.0

    kinds=["flat"]*len(pts)
    for i in range(1,len(pts)):
        vr=(pts[i][3]-pts[i-1][3])/dt[i] if dt[i]>0 else 0.0
        if vr > vr_thr:
            kinds[i]="ascent"
        elif vr < -vr_thr:
            kinds[i]="descent"
        else:
            kinds[i]=kinds[i-1] if i>1 else "flat"

    segs=[]
    start=0
    cur=kinds[1] if len(pts)>1 else "flat"
    for i in range(2,len(pts)):
        if kinds[i] != cur:
            segs.append((start,i-1,cur))
            start=i-1
            cur=kinds[i]
    segs.append((start,len(pts)-1,cur))

    # merge too-short
    merged=[]
    for s in segs:
        if not merged:
            merged.append(s); continue
        a0,a1,ak=merged[-1]
        b0,b1,bk=s
        dur=pts[b1][0]-pts[b0][0]
        if dur < min_seg_s:
            merged[-1]=(a0,b1,ak)
        else:
            merged.append(s)

    return merged if merged else [(0,len(pts)-1,"descent")]

def reduce_to_target(pts, segs, target=24):
    segs=list(segs)
    def dur(seg):
        i0,i1,_=seg
        return pts[i1][0]-pts[i0][0]
    while len(segs) > target and len(segs) > 1:
        idx=min(range(len(segs)), key=lambda i: dur(segs[i]))
        i0,i1,_=segs[idx]
        if idx > 0:
            p0,p1,pk=segs[idx-1]
            segs[idx-1]=(p0,i1,pk)
            segs.pop(idx)
        else:
            n0,n1,nk=segs[1]
            segs[1]=(i0,n1,nk)
            segs.pop(0)
    return segs

def seg_stats(pts, dist, dt, speed, grade, i0, i1):
    start=pts[i0][0]; end=pts[i1][0]
    dur=max(0.0,end-start)
    seg_dist=sum(dist[i0+1:i1+1])
    seg_max_speed=max(speed[i0+1:i1+1] or [0.0])
    elev_change=pts[i1][3]-pts[i0][3]
    avg_speed=(seg_dist/dur) if dur>0 else 0.0
    avg_grade=(abs(elev_change)/seg_dist*100.0) if seg_dist>0 else 0.0
    max_grade=max(grade[i0+1:i1+1] or [0.0])
    alts=[p[3] for p in pts[i0:i1+1]]
    return dict(
        start=start,end=end,dur=dur,dist=seg_dist,avg_speed=avg_speed,max_speed=seg_max_speed,
        elev_change=elev_change,avg_grade=avg_grade,max_grade=max_grade,
        start_alt=pts[i0][3],end_alt=pts[i1][3],max_alt=max(alts),min_alt=min(alts),
    )

def build_segment_csv(pts, dist, dt, speed, grade, segs):
    lines=["4"]
    asc_id=0; des_id=0
    for (i0,i1,kind) in segs:
        st=seg_stats(pts, dist, dt, speed, grade, i0, i1)
        if kind=="ascent":
            seg_type="10"; subtype="ski_lift"
            asc_id += 1; sid=asc_id
            label=f"Remonte-pente {asc_id}"
        else:
            seg_type="8"; subtype="ski_run"
            des_id += 1; sid=des_id
            label=f"Descente à ski {des_id}"
        row=[
            f"{st['start']:.6f}",
            f"{st['end']:.6f}",
            seg_type,
            "0",
            str(sid),
            label,
            "",
            subtype,
            "",
            "",
            "",
            f"{st['dur']:.3f}",
            f"{st['avg_speed']:.6f}",
            f"{st['dist']:.6f}",
            f"{st['elev_change']:.6f}",
            f"{st['max_speed']:.6f}",
            f"{st['avg_grade']:.6f}",
            f"{st['max_grade']:.6f}",
            f"{st['start_alt']:.2f}",
            f"{st['max_alt']:.2f}",
            f"{st['min_alt']:.2f}",
            f"{st['end_alt']:.2f}",
        ]
        if len(row) != 22:
            raise SystemExit(f"Segment row has {len(row)} cols")
        lines.append(",".join(row))
    return lines, asc_id, des_id

def build_events_xml(start_ts, end_ts):
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
        '<events>\n'
        f'\t<event start="{iso_z(start_ts)}" end="{iso_z(start_ts)}" type="start"/>\n'
        f'\t<event start="{iso_z(end_ts)}" end="{iso_z(end_ts)}" type="stop"/>\n'
        '</events>\n'
    )

def build_track_xml(name, desc, tz, pts, dist, dt, speed, grade, segs, asc_id, des_id):
    start_ts=pts[0][0]; end_ts=pts[-1][0]
    total_dist=sum(dist)
    total_dur=end_ts-start_ts

    ascent=0.0; descent=0.0
    for i in range(1,len(pts)):
        de=pts[i][3]-pts[i-1][3]
        if de>0: ascent += de
        else: descent += -de

    max_sp=max(speed or [0.0])
    min_alt=min(p[3] for p in pts)
    max_alt=max(p[3] for p in pts)

    # Segment split totals
    ascent_dist=0.0; descent_dist=0.0
    ascent_time=0.0; descent_time=0.0
    max_asc_sp=0.0; max_des_sp=0.0
    max_asc_grade=0.0; max_des_grade=0.0

    for (i0,i1,kind) in segs:
        st=seg_stats(pts, dist, dt, speed, grade, i0, i1)
        if kind=="ascent":
            ascent_dist += st["dist"]; ascent_time += st["dur"]
            max_asc_sp = max(max_asc_sp, st["max_speed"])
            max_asc_grade = max(max_asc_grade, st["max_grade"])
        else:
            descent_dist += st["dist"]; descent_time += st["dur"]
            max_des_sp = max(max_des_sp, st["max_speed"])
            max_des_grade = max(max_des_grade, st["max_grade"])

    avg_speed=(total_dist/total_dur) if total_dur>0 else 0.0
    avg_asc_sp=(ascent_dist/ascent_time) if ascent_time>0 else 0.0
    avg_des_sp=(descent_dist/descent_time) if descent_time>0 else 0.0

    moving_time=0.0
    for i in range(1,len(pts)):
        if speed[i] > 0.5:
            moving_time += dt[i]
    moving_avg=(total_dist/moving_time) if moving_time>0 else 0.0

    # Metrics order as Ski Tracks export
    order=[
        "maxspeed","maxdescentspeed","maxascentspeed","maxdescentsteepness","maxascentsteepness",
        "totalascent","totaldescent","maxaltitude","minaltitude","distance","descentdistance","ascentdistance",
        "averagespeed","averagedescentspeed","averageascentspeed","movingaveragespeed","duration",
        "startaltitude","finishaltitude","ascents","descents"
    ]
    vals={
        "maxspeed":max_sp,
        "maxdescentspeed":max_des_sp,
        "maxascentspeed":max_asc_sp,
        "maxdescentsteepness":max_des_grade,
        "maxascentsteepness":max_asc_grade,
        "totalascent":ascent,
        "totaldescent":descent,
        "maxaltitude":max_alt,
        "minaltitude":min_alt,
        "distance":total_dist,
        "descentdistance":descent_dist,
        "ascentdistance":ascent_dist,
        "averagespeed":avg_speed,
        "averagedescentspeed":avg_des_sp,
        "averageascentspeed":avg_asc_sp,
        "movingaveragespeed":moving_avg,
        "duration":moving_time,
        "startaltitude":pts[0][3],
        "finishaltitude":pts[-1][3],
        "ascents":asc_id,
        "descents":des_id,
    }

    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
        '<track\n'
        '\tversion="1.0"\n'
        f'\tname="{escape(name)}"\n'
        f'\tdescription="{escape(desc)}"\n'
        '\tactivity="skiing"\n'
        f'\tstart="{iso_with_tz(start_ts, tz)}"\n'
        f'\tfinish="{iso_with_tz(end_ts, tz)}"\n'
        f'\ttz="{tz}"\n'
        f'\tduration="{(end_ts-start_ts):.3f}"\n'
        '\tplatform="gpx_to_skiz_v7.1_fixed2"\n'
        f'\tsyncIdentifier="{uuid.uuid4()}"\n'
        '\tsyncVersion="1">\n'
        '\t<extensions>\n'
        '\t\t<location>\n'
        f'\t\t\t<lat>{pts[0][1]:.14f}</lat>\n'
        f'\t\t\t<lon>{pts[0][2]:.14f}</lon>\n'
        '\t\t</location>\n'
        '\t</extensions>\n'
        '\t<metrics>\n'
    )
    for k in order:
        v=vals[k]
        if k in ("ascents","descents"):
            xml += f"\t\t<{k}>{int(v)}</{k}>\n"
        elif k == "duration":
            xml += f"\t\t<{k}>{float(v):.3f}</{k}>\n"
        elif k in ("totalascent","totaldescent","maxaltitude","minaltitude","distance","descentdistance","ascentdistance","startaltitude","finishaltitude"):
            xml += f"\t\t<{k}>{float(v):.1f}</{k}>\n"
        elif k in ("maxdescentsteepness","maxascentsteepness"):
            xml += f"\t\t<{k}>{float(v):.1f}</{k}>\n"
        else:
            xml += f"\t\t<{k}>{float(v):.2f}</{k}>\n"
    xml += "\t</metrics>\n</track>\n"
    return xml

def write_skiz(out_skiz, files_dict):
    tmp = tempfile.mkdtemp()
    try:
        for fname, content in files_dict.items():
            with open(os.path.join(tmp, fname), "w", encoding="utf-8", newline="\n") as f:
                f.write(content)
        with zipfile.ZipFile(out_skiz, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for fname in ["Segment.csv","Events.xml","RelativeAltitudeSensor.csv","Samples.xml","Photos.csv","Track.xml","RawLocations.csv","Nodes.csv"]:
                z.write(os.path.join(tmp, fname), arcname=fname)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

# ---------- Build depending on input ----------
ext = os.path.splitext(INPATH)[1].lower()

if ext == ".slopes":
    name, desc, tz, pts, gps_lines, raw_lines, actions = read_slopes(INPATH)
    pts.sort(key=lambda x: x[0])
    dist, dt, spd, crs, grd = step_metrics(pts)

    # For slopes, segments come from Metadata.xml actions when available
    segs=[]
    if actions:
        # Map action start/end times to nearest point indices
        ts_list=[p[0] for p in pts]
        import bisect
        for st,en,typ in actions:
            i0=max(0, bisect.bisect_left(ts_list, st) - 1)
            i1=min(len(pts)-1, bisect.bisect_right(ts_list, en) - 1)
            kind = "ascent" if typ.lower()=="lift" else "descent"
            if i1 > i0:
                segs.append((i0,i1,kind))
    if not segs:
        segs = reduce_to_target(pts, segmentize_gpx(pts, dt, mode="v7.1"), target=24)

    seg_csv, asc_id, des_id = build_segment_csv(pts, dist, dt, spd, grd, segs)

    # Nodes: reformat Slopes GPS.csv into repo-style (stable under refresh)
    nodes_lines = build_nodes_repo_style(pts, crs, spd)

    # RawLocations: use last 2 cols from Slopes GPS as "accuracies" if present; add ts2
    # Slopes GPS.csv has 8 cols: ts,lat,lon,alt,course,speed,a,b -> we map a->hacc, b->vacc
    hacc=46.0; vacc=7.58
    if gps_lines:
        parts=gps_lines[0].split(",")
        if len(parts) >= 8:
            try:
                hacc=float(parts[6]); vacc=float(parts[7])
            except Exception:
                pass
    rawloc_lines = build_rawlocations_9cols(pts, crs, spd, hacc=hacc, vacc=vacc)

    ras_lines = build_relative_alt_sensor(pts)
    events_xml = build_events_xml(pts[0][0], pts[-1][0])
    track_xml = build_track_xml(name, desc, tz, pts, dist, dt, spd, grd, segs, asc_id, des_id)

else:
    # GPX (default)
    name, desc, pts = parse_gpx(INPATH)
    pts.sort(key=lambda x: x[0])
    tz = TZ_GPX
    dist, dt, spd, crs, grd = step_metrics(pts)

    segs = segmentize_gpx(pts, dt, mode=MODE)
    segs = reduce_to_target(pts, segs, target=24)

    seg_csv, asc_id, des_id = build_segment_csv(pts, dist, dt, spd, grd, segs)
    nodes_lines = build_nodes_repo_style(pts, crs, spd)
    rawloc_lines = build_rawlocations_9cols(pts, crs, spd)
    ras_lines = build_relative_alt_sensor(pts)
    events_xml = build_events_xml(pts[0][0], pts[-1][0])
    track_xml = build_track_xml(name, desc, tz, pts, dist, dt, spd, grd, segs, asc_id, des_id)

samples_xml = '<?xml version="1.0" encoding="UTF-8"?>\n<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n<samples/>\n'

files = {
    "Segment.csv": "\n".join(seg_csv) + "\n",
    "Events.xml": events_xml,
    "RelativeAltitudeSensor.csv": "\n".join(ras_lines) + "\n",
    "Samples.xml": samples_xml,
    "Photos.csv": "",
    "Track.xml": track_xml,
    "RawLocations.csv": "\n".join(rawloc_lines) + "\n",
    "Nodes.csv": "\n".join(nodes_lines) + "\n",
}

# Validate Track.xml is parseable
ET.fromstring(files["Track.xml"].encode("utf-8"))

write_skiz(OUT_SKIZ, files)

PY
  then
    echo "Python conversion failed; see log." >&2
    exit 1
  fi

  if unzip -tq "$out_skiz" >>"$per_log" 2>&1; then
    log INFO "OK: $out_skiz"
  else
    log ERROR "Zip integrity failed: $out_skiz (see $per_log)"
    return 1
  fi

  return 0
}

fail=0
for f in "${FILES[@]}"; do
  if ! convert_one "$f"; then
    fail=1
  fi
done

if [[ "$fail" -ne 0 ]]; then
  die "Some conversions failed. Check logs in: $LOG_DIR"
fi

log INFO "Done. Output files are in: $OUT_DIR"
