\
#!/usr/bin/env bash
# Convert GPX tracks to Ski Tracks .skiz (zip) files
# Focus: make "Refresh segments" NOT wipe the main stats.
# Strategy: output a Nodes.csv schema aligned with alexbelgium/gpx_to_skiz (gpsbabel style),
# because Ski Tracks appears to recompute totals from Nodes.csv during refresh.
#
# Works on Alpine (apk) and Debian/Ubuntu (apt).
set -euo pipefail
IFS=$'\n\t'

SCRIPT_NAME="$(basename "$0")"
OUT_DIR="gpx_to_skiz"
LOG_DIR="$OUT_DIR/logs"
VERBOSE=1
TZ_OFFSET="+01:00"

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [OPTIONS] [FILE1.gpx FILE2.gpx ...]

Converts GPX files to .skiz (Ski Tracks import).

Options:
  -o, --out-dir DIR      Output directory (default: $OUT_DIR)
  --tz +HH:MM|-HH:MM     Timezone offset used in Track.xml (default: $TZ_OFFSET)
  -q, --quiet            Less console output
  -h, --help             Show help

Examples:
  $SCRIPT_NAME
  $SCRIPT_NAME activity.gpx
  $SCRIPT_NAME --tz +01:00 -o out *.gpx
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

die() { log ERROR "$*"; exit 1; }
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
    die "Unsupported distribution (no apk/apt-get). Please install python3, zip, unzip."
  fi

  for c in python3 zip unzip; do
    need_cmd "$c" || die "Failed to install dependency: $c"
  done
}

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
    -q|--quiet) VERBOSE=0 ;;
    -h|--help) usage; exit 0 ;;
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

log INFO "Starting $SCRIPT_NAME"
log INFO "Output directory: $OUT_DIR"
log INFO "Timezone offset: $TZ_OFFSET"
log INFO "Main log: $MAIN_LOG"

install_deps

GPX_FILES=()
if [[ "${#ARGS[@]}" -eq 0 ]]; then
  while IFS= read -r f; do GPX_FILES+=("$f"); done < <(find . -maxdepth 1 -type f -name '*.gpx' -print | sort)
else
  GPX_FILES=("${ARGS[@]}")
fi
[[ "${#GPX_FILES[@]}" -gt 0 ]] || die "No GPX files found."

convert_one() {
  local gpx="$1"
  [[ -f "$gpx" ]] || { log ERROR "File not found: $gpx"; return 1; }

  local base per_log out_skiz
  base="$(basename "$gpx")"
  base="${base%.gpx}"
  out_skiz="$OUT_DIR/$base.skiz"
  per_log="$LOG_DIR/${base}_$(date -u '+%Y%m%dT%H%M%SZ').log"
  : >"$per_log"
  mkdir -p "$OUT_DIR"

  log INFO "Converting: $gpx -> $out_skiz"

  if ! python3 - "$gpx" "$out_skiz" "$TZ_OFFSET" >>"$per_log" 2>&1 <<'PY'
import sys, os, math, uuid, datetime, tempfile, shutil, zipfile, random, re
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

gpx_path = sys.argv[1]
out_skiz  = sys.argv[2]
tz_str    = sys.argv[3]  # +HH:MM or -HH:MM

def tz_minutes(s):
    m = re.match(r'^([+-])(\d\d):(\d\d)$', s)
    if not m:
        return 0
    sign = 1 if m.group(1) == '+' else -1
    return sign * (int(m.group(2))*60 + int(m.group(3)))

tzinfo = datetime.timezone(datetime.timedelta(minutes=tz_minutes(tz_str)))

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000.0
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dphi = math.radians(lat2-lat1); dl = math.radians(lon2-lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(a))

def bearing(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dl = math.radians(lon2-lon1)
    y = math.sin(dl)*math.cos(phi2)
    x = math.cos(phi1)*math.sin(phi2) - math.sin(phi1)*math.cos(phi2)*math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360

def iso_z(ts):
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc)
    return dt.isoformat(timespec='milliseconds').replace('+00:00','Z')

def iso_with_tz(ts):
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc).astimezone(tzinfo)
    return dt.isoformat(timespec='milliseconds')

def parse_gpx(path):
    raw = open(path,'rb').read()
    root = ET.fromstring(raw)
    ns={}
    if root.tag.startswith('{'):
        ns_uri = root.tag.split('}')[0][1:]
        ns={'g':ns_uri}
    trk = root.find('g:trk', ns) if ns else root.find('trk')
    if trk is None:
        raise SystemExit("No <trk> in GPX")
    name = (trk.findtext('g:name', default='', namespaces=ns) if ns else trk.findtext('name', default='')) or os.path.splitext(os.path.basename(path))[0]
    desc = (trk.findtext('g:desc', default='', namespaces=ns) if ns else trk.findtext('desc', default='')) or ''
    pts=[]
    for seg in (trk.findall('g:trkseg', ns) if ns else trk.findall('trkseg')):
        for pt in (seg.findall('g:trkpt', ns) if ns else seg.findall('trkpt')):
            lat = float(pt.attrib.get('lat'))
            lon = float(pt.attrib.get('lon'))
            ele_txt  = (pt.findtext('g:ele', default=None, namespaces=ns) if ns else pt.findtext('ele', default=None))
            time_txt = (pt.findtext('g:time', default=None, namespaces=ns) if ns else pt.findtext('time', default=None))
            if not time_txt:
                continue
            dt = datetime.datetime.fromisoformat(time_txt.replace('Z', '+00:00'))
            ts = dt.timestamp()
            ele = float(ele_txt) if ele_txt is not None else 0.0
            pts.append((ts,lat,lon,ele))
    pts.sort(key=lambda x:x[0])
    if len(pts) < 2:
        raise SystemExit(f"Not enough timed points in GPX: {len(pts)}")
    return name, desc, pts

def downsample_time(pts, step=3.0):
    # keep first, then keep point whenever >= step seconds elapsed
    out=[pts[0]]
    last_t=pts[0][0]
    for p in pts[1:]:
        if p[0] - last_t >= step - 1e-6:
            out.append(p)
            last_t=p[0]
    if out[-1][0] != pts[-1][0]:
        out.append(pts[-1])
    return out

def compute_course_speed(pts):
    # returns course,speed for each point, using forward difference (like gpsbabel track,speed,course tends to do)
    n=len(pts)
    course=[-1.0]*n
    speed=[-1.0]*n
    for i in range(n-1):
        t0,la0,lo0,_=pts[i]
        t1,la1,lo1,_=pts[i+1]
        dt=max(0.0, t1-t0)
        d=haversine(la0,lo0,la1,lo1)
        course[i]=bearing(la0,lo0,la1,lo1) if d>0 else (course[i-1] if i>0 else -1.0)
        speed[i]=(d/dt) if dt>0 else 0.0
    course[-1]=course[-2] if n>1 else -1.0
    speed[-1]=speed[-2] if n>1 else -1.0
    return course,speed

def pressure_kpa_from_alt_m(alt_m, p0=101.325):
    T0=288.15; L=0.0065; g=9.80665; M=0.0289644; R=8.3144598
    h=max(-500.0, alt_m)
    factor = 1.0 - (L*h)/T0
    if factor <= 0:
        factor = 1e-6
    exp = (g*M)/(R*L)
    return p0*(factor**exp)

def build_relative_alt_sensor(pts_raw):
    start_ts=pts_raw[0][0]; end_ts=pts_raw[-1][0]
    interval=1.0
    lag=1.65
    # simple linear interpolation on altitude
    t_arr=[p[0] for p in pts_raw]
    e_arr=[p[3] for p in pts_raw]
    def alt_at(t):
        if t <= t_arr[0]: return e_arr[0]
        if t >= t_arr[-1]: return e_arr[-1]
        import bisect
        j=bisect.bisect_right(t_arr, t)-1
        t0,t1=t_arr[j],t_arr[j+1]
        e0,e1=e_arr[j],e_arr[j+1]
        w=(t-t0)/(t1-t0) if t1!=t0 else 0.0
        return e0*(1-w)+e1*w

    times=[]
    t=start_ts+3.5
    while t <= end_ts + 1e-6:
        times.append(t)
        t += interval
    if not times:
        times=[start_ts]
    alt0=alt_at(times[0])
    lines=[]
    for t in times:
        alt=alt_at(t)
        rel=alt-alt0
        p=pressure_kpa_from_alt_m(alt) + random.uniform(-0.02,0.02)
        lines.append(f"{iso_z(t)},{iso_z(max(start_ts, t-lag))},{p:.4f},{rel:.4f}")
    return lines

def segmentize_simple(pts_nodes, course, speed):
    # minimal segments: one "ski_run" spanning all (Ski Tracks can refresh anyway)
    start=pts_nodes[0][0]; end=pts_nodes[-1][0]
    # columns based on observed 22-col Segment.csv rows
    row=[
        f"{start:.6f}", f"{end:.6f}",
        "8","0","1","Ski segment 1","",
        "ski_run","","","",
        f"{(end-start):.3f}",
        "0.000000","0.000000","0.000000","0.000000",
        "0.000000","0.000000",
        f"{pts_nodes[0][3]:.2f}",
        f"{max(p[3] for p in pts_nodes):.2f}",
        f"{min(p[3] for p in pts_nodes):.2f}",
        f"{pts_nodes[-1][3]:.2f}",
    ]
    return ["4", ",".join(row)]

name, desc, pts_raw = parse_gpx(gpx_path)

# Build Nodes in "gpsbabel style" (repo nodes.style):
#  time as integer seconds with .000, lat/lon 6 decimals, altitude 2 decimals,
#  course float, speed 6 decimals, GPS precision (vdop) -> 0 if unknown, last field constant 0.25
pts_nodes = downsample_time(pts_raw, step=3.0)
course, speed = compute_course_speed(pts_nodes)

nodes_lines=[]
for i,(ts,la,lo,el) in enumerate(pts_nodes):
    t_int=int(round(ts))
    nodes_lines.append(
        f"{t_int:d}.000,{la:.6f},{lo:.6f},{el:.2f},{course[i]:.2f},{speed[i]:.6f},0.000000,0.250000"
    )

# RawLocations: keep original GPX points; add plausible accuracy + "system time" column
course_r, speed_r = compute_course_speed(pts_raw)
raw_lines=[]
for i,(ts,la,lo,el) in enumerate(pts_raw):
    # accuracy: if GPX has none, use a plausible constant (Ski Tracks doesn't seem strict here)
    hacc=25.0
    vacc=7.5
    sys_ts = ts + 0.009
    raw_lines.append(
        f"{ts:.3f},{la:.9f},{lo:.9f},{el:.3f},{course_r[i]:.2f},{speed_r[i]:.2f},{hacc:.3f},{vacc:.3f},{sys_ts:.3f}"
    )

ras_lines = build_relative_alt_sensor(pts_raw)
seg_lines = segmentize_simple(pts_nodes, course, speed)

# Events.xml: start/stop
start_ts=pts_raw[0][0]; end_ts=pts_raw[-1][0]
events_xml = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
    '<events>\n'
    f'\t<event start="{iso_z(start_ts)}" end="{iso_z(start_ts)}" type="start"/>\n'
    f'\t<event start="{iso_z(end_ts)}" end="{iso_z(end_ts)}" type="stop"/>\n'
    '</events>\n'
)
samples_xml = '<?xml version="1.0" encoding="UTF-8"?>\n<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n<samples/>\n'
photos_csv=""

# Track.xml: keep a full metrics set so import looks correct,
# but make sure the underlying Nodes schema is compatible with refresh.
# Compute basic totals from raw points.
dist=0.0
ascent=0.0; descent=0.0
max_speed=0.0
for i in range(1,len(pts_raw)):
    d=haversine(pts_raw[i-1][1],pts_raw[i-1][2],pts_raw[i][1],pts_raw[i][2])
    dist += d
    dt=max(0.0, pts_raw[i][0]-pts_raw[i-1][0])
    if dt>0:
        max_speed=max(max_speed, d/dt)
    de=pts_raw[i][3]-pts_raw[i-1][3]
    if de>0: ascent += de
    else: descent += -de

duration_total=end_ts-start_ts
avg_speed=(dist/duration_total) if duration_total>0 else 0.0
min_alt=min(p[3] for p in pts_raw); max_alt=max(p[3] for p in pts_raw)

track_xml = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
    '<track\n'
    '\tversion="1.0"\n'
    f'\tname="{escape(name)}"\n'
    f'\tdescription="{escape(desc)}"\n'
    '\tactivity="skiing"\n'
    f'\tstart="{iso_with_tz(start_ts)}"\n'
    f'\tfinish="{iso_with_tz(end_ts)}"\n'
    f'\ttz="{tz_str}"\n'
    f'\tduration="{duration_total:.3f}"\n'
    '\tplatform="iPhone/ios/SkiTracks"\n'
    f'\tsyncIdentifier="{str(uuid.uuid4()).upper()}"\n'
    '\tsyncVersion="1">\n'
    '\t<extensions>\n'
    '\t\t<location>\n'
    f'\t\t\t<lat>{pts_raw[0][1]:.14f}</lat>\n'
    f'\t\t\t<lon>{pts_raw[0][2]:.14f}</lon>\n'
    '\t\t</location>\n'
    '\t</extensions>\n'
    '\t<metrics>\n'
    f'\t\t<maxspeed>{max_speed*3.6/3.6:.2f}</maxspeed>\n'
    f'\t\t<maxdescentspeed>{max_speed:.2f}</maxdescentspeed>\n'
    f'\t\t<maxascentspeed>{max_speed:.2f}</maxascentspeed>\n'
    '\t\t<maxdescentsteepness>0.0</maxdescentsteepness>\n'
    '\t\t<maxascentsteepness>0.0</maxascentsteepness>\n'
    f'\t\t<totalascent>{ascent:.1f}</totalascent>\n'
    f'\t\t<totaldescent>{descent:.1f}</totaldescent>\n'
    f'\t\t<maxaltitude>{max_alt:.1f}</maxaltitude>\n'
    f'\t\t<minaltitude>{min_alt:.1f}</minaltitude>\n'
    f'\t\t<distance>{dist:.1f}</distance>\n'
    '\t\t<descentdistance>0.0</descentdistance>\n'
    '\t\t<ascentdistance>0.0</ascentdistance>\n'
    f'\t\t<averagespeed>{avg_speed:.2f}</averagespeed>\n'
    '\t\t<averagedescentspeed>0.00</averagedescentspeed>\n'
    '\t\t<averageascentspeed>0.00</averageascentspeed>\n'
    f'\t\t<movingaveragespeed>{avg_speed:.2f}</movingaveragespeed>\n'
    f'\t\t<duration>{duration_total:.3f}</duration>\n'
    f'\t\t<startaltitude>{pts_raw[0][3]:.1f}</startaltitude>\n'
    f'\t\t<finishaltitude>{pts_raw[-1][3]:.1f}</finishaltitude>\n'
    '\t\t<ascents>0</ascents>\n'
    '\t\t<descents>0</descents>\n'
    '\t</metrics>\n'
    '</track>\n'
)

tmpdir=tempfile.mkdtemp()
try:
    def w(fname, content):
        with open(os.path.join(tmpdir, fname), 'w', encoding='utf-8', newline='\n') as f:
            f.write(content)

    w("Nodes.csv", "\n".join(nodes_lines) + "\n")
    w("RawLocations.csv", "\n".join(raw_lines) + "\n")
    w("RelativeAltitudeSensor.csv", "\n".join(ras_lines) + "\n")
    w("Segment.csv", "\n".join(seg_lines) + "\n")
    w("Events.xml", events_xml)
    w("Samples.xml", samples_xml)
    w("Photos.csv", photos_csv)
    w("Track.xml", track_xml)

    with zipfile.ZipFile(out_skiz, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for fname in ["Segment.csv","Events.xml","RelativeAltitudeSensor.csv","Samples.xml","Photos.csv","Track.xml","RawLocations.csv","Nodes.csv"]:
            z.write(os.path.join(tmpdir, fname), arcname=fname)
finally:
    shutil.rmtree(tmpdir, ignore_errors=True)
PY
  then
    log ERROR "Conversion failed for $gpx. See: $per_log"
    return 1
  fi

  if unzip -tq "$out_skiz" >>"$per_log" 2>&1; then
    log INFO "Created: $out_skiz (OK)"
  else
    log ERROR "Zip integrity test failed: $out_skiz (see $per_log)"
    return 1
  fi
}

fail=0
for f in "${GPX_FILES[@]}"; do
  if ! convert_one "$f"; then
    fail=1
  fi
done

if [[ "$fail" -ne 0 ]]; then
  die "Some conversions failed. Check logs in: $LOG_DIR"
fi

log INFO "Done. Output files are in: $OUT_DIR"
