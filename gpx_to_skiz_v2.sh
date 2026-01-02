#!/usr/bin/env bash
# GPX/.slopes -> Ski Tracks .skiz
# - GPX: keep original trackpoints (no resample), never output negative course/speed,
#        Nodes.csv matches repo-style (epoch seconds .000, vdop=0.0, last=0.25).
# - .slopes: reads GPS.csv/RawGPS.csv/Metadata.xml (actions -> segments)
# Works on Alpine (apk) and Debian/Ubuntu (apt)
set -euo pipefail
IFS=$'\n\t'

SCRIPT_NAME="$(basename "$0")"
OUT_DIR="gpx_to_skiz"
LOG_DIR=""
VERBOSE=1
TZ_OFFSET="+01:00"

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [OPTIONS] [FILES...]

Converts GPX (.gpx) and Slopes exports (.slopes) into Ski Tracks .skiz.

Options:
  -o, --out-dir DIR      Output directory (default: $OUT_DIR)
  --tz +HH:MM|-HH:MM     Timezone offset for GPX Track.xml (default: $TZ_OFFSET)
  -q, --quiet            Less console output
  -h, --help             Show this help

If no input files are provided, converts all *.gpx and *.slopes in current directory.

Outputs:
  <out-dir>/<name>.skiz
Logs:
  <out-dir>/logs/*.log
EOF
}

need_cmd() { command -v "$1" >/dev/null 2>&1; }

as_root() {
  if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
    "$@"
  elif need_cmd sudo; then
    sudo "$@"
  else
    echo "ERROR: Need root privileges to install dependencies (run as root or install sudo)." >&2
    exit 1
  fi
}

install_deps() {
  local miss=()
  for c in python3 unzip; do
    need_cmd "$c" || miss+=("$c")
  done
  # zip is optional (python uses zipfile), but unzip helps self-check
  if [[ "${#miss[@]}" -eq 0 ]]; then
    return 0
  fi

  if need_cmd apk; then
    as_root apk add --no-cache python3 unzip ca-certificates >/dev/null
  elif need_cmd apt-get; then
    as_root apt-get update -y >/dev/null
    as_root apt-get install -y python3 unzip ca-certificates >/dev/null
  else
    echo "ERROR: Unsupported distro (no apk/apt-get). Install python3+unzip manually." >&2
    exit 1
  fi
}

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out-dir)
      shift
      OUT_DIR="$1"
      ;;
    --tz)
      shift
      TZ_OFFSET="$1"
      ;;
    -q|--quiet)
      VERBOSE=0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do ARGS+=("$1"); shift; done
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      exit 1
      ;;
    *)
      ARGS+=("$1")
      ;;
  esac
  shift || true
done

LOG_DIR="$OUT_DIR/logs"
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

install_deps

# Collect files
FILES=()
if [[ "${#ARGS[@]}" -eq 0 ]]; then
  while IFS= read -r f; do FILES+=("$f"); done < <(find . -maxdepth 1 -type f \( -name '*.gpx' -o -name '*.slopes' \) -print | sort)
else
  FILES=("${ARGS[@]}")
fi

if [[ "${#FILES[@]}" -eq 0 ]]; then
  log ERROR "No .gpx or .slopes found."
  exit 1
fi

convert_one() {
  local inpath="$1"
  [[ -f "$inpath" ]] || { log ERROR "Missing file: $inpath"; return 1; }

  local base ext out_skiz per_log
  base="$(basename "$inpath")"
  ext="${base##*.}"
  base="${base%.*}"
  out_skiz="$OUT_DIR/$base.skiz"
  per_log="$LOG_DIR/${base}_$(date -u '+%Y%m%dT%H%M%SZ').log"
  : >"$per_log"

  mkdir -p "$OUT_DIR"

  log INFO "Converting: $inpath -> $out_skiz"

  if ! python3 - "$inpath" "$out_skiz" "$TZ_OFFSET" >>"$per_log" 2>&1 <<'PY'
import sys, os, math, re, uuid, zipfile, datetime, tempfile, shutil
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

INPATH = sys.argv[1]
OUTSKIZ = sys.argv[2]
DEFAULT_TZ = sys.argv[3]  # used for GPX only

def tzinfo_from_offset(tz_str: str) -> datetime.tzinfo:
    m = re.match(r'^([+-])(\d\d):(\d\d)$', tz_str.strip())
    if not m:
        return datetime.timezone.utc
    sign = 1 if m.group(1) == '+' else -1
    minutes = sign * (int(m.group(2)) * 60 + int(m.group(3)))
    return datetime.timezone(datetime.timedelta(minutes=minutes))

def iso_z(ts: float) -> str:
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc)
    return dt.isoformat(timespec='milliseconds').replace('+00:00', 'Z')

def iso_with_tz(ts: float, tzinfo: datetime.tzinfo) -> str:
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc).astimezone(tzinfo)
    return dt.isoformat(timespec='milliseconds')

def haversine(lat1, lon1, lat2, lon2) -> float:
    R = 6371000.0
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(a))

def bearing(lat1, lon1, lat2, lon2) -> float:
    phi1 = math.radians(lat1); phi2 = math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl)*math.cos(phi2)
    x = math.cos(phi1)*math.sin(phi2) - math.sin(phi1)*math.cos(phi2)*math.cos(dl)
    br = (math.degrees(math.atan2(y, x)) + 360) % 360
    return br

def pressure_kpa_from_alt_m(alt_m: float, p0=101.325) -> float:
    # ISA troposphere approximation
    T0=288.15; L=0.0065; g=9.80665; M=0.0289644; R=8.3144598
    h=max(-500.0, float(alt_m))
    factor = 1.0 - (L*h)/T0
    if factor <= 0:
        factor = 1e-6
    exp = (g*M)/(R*L)
    return float(p0*(factor**exp))

def safe_epoch_int(ts: float, prev: int | None) -> int:
    t = int(round(ts))
    if prev is not None and t <= prev:
        t = prev + 1
    return t

# --------------------
# GPX parsing
# --------------------
def parse_gpx(path: str):
    raw = open(path, 'rb').read()
    root = ET.fromstring(raw)
    ns = {}
    if root.tag.startswith('{'):
        ns_uri = root.tag.split('}')[0][1:]
        ns = {'g': ns_uri}

    trk = root.find('g:trk', ns) if ns else root.find('trk')
    if trk is None:
        raise SystemExit('No <trk> found')

    name = (trk.findtext('g:name', default='', namespaces=ns) if ns else trk.findtext('name', default='')) or os.path.splitext(os.path.basename(path))[0]
    desc = (trk.findtext('g:desc', default='', namespaces=ns) if ns else trk.findtext('desc', default='')) or ''

    pts = []
    trksegs = trk.findall('g:trkseg', ns) if ns else trk.findall('trkseg')
    for seg in trksegs:
        trkpts = seg.findall('g:trkpt', ns) if ns else seg.findall('trkpt')
        for pt in trkpts:
            lat = float(pt.attrib.get('lat'))
            lon = float(pt.attrib.get('lon'))
            ele_txt = (pt.findtext('g:ele', default=None, namespaces=ns) if ns else pt.findtext('ele', default=None))
            time_txt = (pt.findtext('g:time', default=None, namespaces=ns) if ns else pt.findtext('time', default=None))
            if time_txt is None:
                continue
            dt = datetime.datetime.fromisoformat(time_txt.replace('Z', '+00:00'))
            ts = dt.timestamp()
            ele = float(ele_txt) if ele_txt is not None else float('nan')
            pts.append((ts, lat, lon, ele))

    pts.sort(key=lambda x: x[0])
    if len(pts) < 2:
        raise SystemExit(f'Not enough points: {len(pts)}')

    # Fill NaN elevations by linear interpolation over time
    times = [p[0] for p in pts]
    eles  = [p[3] for p in pts]
    if any(math.isnan(e) for e in eles):
        ok = [(t,e) for t,e in zip(times, eles) if not math.isnan(e)]
        if len(ok) >= 2:
            ok_t = [t for t,_ in ok]
            ok_e = [e for _,e in ok]
            import bisect
            def interp(t):
                if t <= ok_t[0]: return ok_e[0]
                if t >= ok_t[-1]: return ok_e[-1]
                j = bisect.bisect_right(ok_t, t) - 1
                t0,t1 = ok_t[j], ok_t[j+1]
                e0,e1 = ok_e[j], ok_e[j+1]
                w = (t - t0) / (t1 - t0) if t1 != t0 else 0.0
                return e0*(1-w) + e1*w
            pts = [(t,la,lo,(interp(t) if math.isnan(e) else e)) for (t,la,lo,e) in pts]
        else:
            pts = [(t,la,lo,(0.0 if math.isnan(e) else e)) for (t,la,lo,e) in pts]

    return name, desc, pts

def compute_course_speed(pts):
    # returns course_deg, speed_m_s lists aligned to pts
    n = len(pts)
    course = [0.0]*n
    speed  = [0.0]*n
    for i in range(1, n):
        t0,la0,lo0,_ = pts[i-1]
        t1,la1,lo1,_ = pts[i]
        d = haversine(la0,lo0,la1,lo1)
        dt = max(0.0, t1 - t0)
        sp = (d/dt) if dt > 0 else 0.0
        br = bearing(la0,lo0,la1,lo1) if d > 0 else course[i-1]
        course[i] = br
        speed[i] = sp
    # First point: copy from second point (never negative)
    course[0] = course[1]
    speed[0] = speed[1] if speed[1] >= 0 else 0.0
    return course, speed

def build_segments_from_alt(pts, min_seg_s=90.0, vr_thr=0.12, target=24):
    # heuristic: ascent vs descent based on vertical rate
    n = len(pts)
    kinds = ['flat'] * n
    for i in range(1, n):
        t0,_,_,e0 = pts[i-1]
        t1,_,_,e1 = pts[i]
        dt = max(0.0, t1 - t0)
        if dt <= 0:
            kinds[i] = kinds[i-1]
            continue
        vr = (e1 - e0) / dt
        if vr > vr_thr:
            kinds[i] = 'ascent'
        elif vr < -vr_thr:
            kinds[i] = 'descent'
        else:
            kinds[i] = kinds[i-1] if i > 1 else 'flat'

    segs=[]
    start=0
    cur=kinds[1] if n > 1 else 'descent'
    for i in range(2,n):
        if kinds[i] != cur:
            segs.append((start, i-1, cur))
            start = i-1
            cur = kinds[i]
    segs.append((start, n-1, cur))

    # merge short
    merged=[]
    for s in segs:
        if not merged:
            merged.append(s); continue
        a0,a1,ak = merged[-1]
        b0,b1,bk = s
        dur = pts[b1][0] - pts[b0][0]
        if dur < min_seg_s:
            merged[-1] = (a0, b1, ak)
        else:
            merged.append(s)

    segs = merged if merged else [(0,n-1,'descent')]

    # reduce to target by merging shortest
    def dur(seg): 
        i0,i1,_ = seg
        return pts[i1][0] - pts[i0][0]
    while len(segs) > target and len(segs) > 1:
        idx = min(range(len(segs)), key=lambda i: dur(segs[i]))
        i0,i1,_ = segs[idx]
        if idx > 0:
            p0,p1,pk = segs[idx-1]
            segs[idx-1] = (p0, i1, pk)
            segs.pop(idx)
        else:
            n0,n1,nk = segs[1]
            segs[1] = (i0, n1, nk)
            segs.pop(0)

    return segs

def seg_stats(pts, course, speed, i0, i1):
    # compute dist, max speed, alt stats
    start = pts[i0][0]
    end = pts[i1][0]
    dur = max(0.0, end - start)
    dist = 0.0
    max_speed = 0.0
    for i in range(i0+1, i1+1):
        dist += haversine(pts[i-1][1], pts[i-1][2], pts[i][1], pts[i][2])
        max_speed = max(max_speed, speed[i])
    elev_change = pts[i1][3] - pts[i0][3]
    avg_speed = (dist/dur) if dur > 0 else 0.0
    avg_grade = (abs(elev_change)/dist*100.0) if dist > 0 else 0.0
    # rough max grade using consecutive points
    max_grade = 0.0
    for i in range(i0+1, i1+1):
        d = haversine(pts[i-1][1], pts[i-1][2], pts[i][1], pts[i][2])
        if d > 0:
            de = pts[i][3] - pts[i-1][3]
            max_grade = max(max_grade, abs(de)/d*100.0)
    alts = [p[3] for p in pts[i0:i1+1]]
    return dict(start=start,end=end,dur=dur,dist=dist,avg_speed=avg_speed,max_speed=max_speed,
                elev_change=elev_change,avg_grade=avg_grade,max_grade=max_grade,
                start_alt=pts[i0][3],end_alt=pts[i1][3],max_alt=max(alts),min_alt=min(alts))

def build_ras(pts):
    # 4 cols like Ski Tracks: t1_iso,t2_iso,pressure_kpa,relative_alt
    interval = 1.315
    lag = 1.650
    times = [p[0] for p in pts]
    alts  = [p[3] for p in pts]
    import bisect
    def alt_at(t):
        if t <= times[0]: return alts[0]
        if t >= times[-1]: return alts[-1]
        j = bisect.bisect_right(times, t) - 1
        t0,t1 = times[j], times[j+1]
        a0,a1 = alts[j], alts[j+1]
        w = (t - t0) / (t1 - t0) if t1 != t0 else 0.0
        return a0*(1-w) + a1*w

    start_ts = times[0]
    end_ts = times[-1]
    t = start_ts + 1.3
    if t > end_ts:
        t = start_ts
    base_alt = alt_at(t)
    lines=[]
    while t <= end_ts + 1e-6:
        a = alt_at(t)
        rel = a - base_alt
        p = pressure_kpa_from_alt_m(a)
        # deterministic tiny jitter to mimic sensor (no random)
        p += ( (int(t*1000) % 7) - 3 ) * 0.002
        lines.append(f"{iso_z(t)},{iso_z(max(start_ts, t-lag))},{p:.4f},{rel:.4f}")
        t += interval
    return lines

# --------------------
# .slopes parsing
# --------------------
def parse_slopes(path: str):
    with zipfile.ZipFile(path, 'r') as z:
        gps = z.read('GPS.csv').decode('utf-8', errors='replace').splitlines()
        rawgps = z.read('RawGPS.csv').decode('utf-8', errors='replace').splitlines()
        meta = z.read('Metadata.xml').decode('utf-8', errors='replace')
    # GPS.csv columns: ts,lat,lon,alt,course,speed,hacc,vacc
    gps_rows = []
    for ln in gps:
        if not ln.strip(): 
            continue
        cols = ln.split(',')
        if len(cols) < 8:
            continue
        gps_rows.append([float(cols[0]), float(cols[1]), float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7])])
    raw_rows = []
    for ln in rawgps:
        if not ln.strip():
            continue
        cols = ln.split(',')
        if len(cols) < 8:
            continue
        raw_rows.append([float(cols[0]), float(cols[1]), float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7])])

    mroot = ET.fromstring(meta)
    # Activity attributes
    attrs = mroot.attrib
    loc_name = attrs.get('locationName', 'Slopes')
    activity_name = loc_name
    desc = f"Imported from .slopes ({loc_name})"
    tz_off_hours = attrs.get('timeZoneOffset', '1')
    # build +HH:MM from hours
    try:
        h = int(float(tz_off_hours))
    except Exception:
        h = 0
    tz_str = f"{h:+03d}:00".replace('+0', '+00').replace('-0', '-00') if abs(h) < 10 else f"{h:+03d}:00"
    # normalize to +HH:MM
    tz_str = f"{h:+03d}:00"
    tz_str = tz_str.replace('+', '+').replace('-', '-')
    if len(tz_str) == 5:
        tz_str = tz_str[0:3] + ':' + tz_str[3:5]

    tzinfo = tzinfo_from_offset(tz_str)

    # Parse Activity start/end (format: '2025-12-31 10:52:02 +0100')
    def parse_meta_dt(s: str) -> float:
        dt = datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S %z')
        return dt.timestamp()

    start_ts = parse_meta_dt(attrs['start'])
    end_ts = parse_meta_dt(attrs['end'])

    # Actions -> segments
    actions = mroot.find('actions')
    segs = []
    if actions is not None:
        for a in actions.findall('Action'):
            st = parse_meta_dt(a.attrib['start'])
            en = parse_meta_dt(a.attrib['end'])
            typ = a.attrib.get('type','Run')
            segs.append((st, en, typ, a.attrib))
    return activity_name, desc, tz_str, tzinfo, start_ts, end_ts, gps_rows, raw_rows, segs, attrs

# --------------------
# Build SKIZ writers
# --------------------
def write_skiz(files: dict[str, str], outpath: str):
    tmp = tempfile.mkdtemp()
    try:
        for name, content in files.items():
            with open(os.path.join(tmp, name), 'w', encoding='utf-8', newline='\n') as f:
                f.write(content)
        with zipfile.ZipFile(outpath, 'w', compression=zipfile.ZIP_DEFLATED) as z:
            for name in ['Segment.csv','Events.xml','RelativeAltitudeSensor.csv','Samples.xml','Photos.csv','Track.xml','RawLocations.csv','Nodes.csv']:
                z.write(os.path.join(tmp, name), arcname=name)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

def make_events_xml(start_ts: float, end_ts: float) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
        '<events>\n'
        f'\t<event start="{iso_z(start_ts)}" end="{iso_z(start_ts)}" type="start"/>\n'
        f'\t<event start="{iso_z(end_ts)}" end="{iso_z(end_ts)}" type="stop"/>\n'
        '</events>\n'
    )

SAMPLES_XML = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
    '<samples/>\n'
)

def make_track_xml(name, desc, tz_str, tzinfo, start_ts, end_ts, metrics: dict) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- Created using CCXML v1.1.0 encoder library, © Core Coders Ltd. -->\n'
        '<track\n'
        '\tversion="1.0"\n'
        f'\tname="{escape(name)}"\n'
        f'\tdescription="{escape(desc)}"\n'
        '\tactivity="skiing"\n'
        f'\tstart="{iso_with_tz(start_ts, tzinfo)}"\n'
        f'\tfinish="{iso_with_tz(end_ts, tzinfo)}"\n'
        f'\ttz="{tz_str}"\n'
        f'\tduration="{(end_ts-start_ts):.3f}"\n'
        '\tplatform="gpx_to_skiz_v7_1_fixed2"\n'
        f'\tsyncIdentifier="{uuid.uuid4()}"\n'
        '\tsyncVersion="1">\n'
        '\t<extensions>\n'
        '\t\t<location>\n'
        f'\t\t\t<lat>{metrics["start_lat"]:.14f}</lat>\n'
        f'\t\t\t<lon>{metrics["start_lon"]:.14f}</lon>\n'
        '\t\t</location>\n'
        '\t</extensions>\n'
        '\t<metrics>\n'
        f'\t\t<maxspeed>{metrics["maxspeed"]:.2f}</maxspeed>\n'
        f'\t\t<maxdescentspeed>{metrics["maxdescentspeed"]:.2f}</maxdescentspeed>\n'
        f'\t\t<maxascentspeed>{metrics["maxascentspeed"]:.2f}</maxascentspeed>\n'
        f'\t\t<maxdescentsteepness>{metrics["maxdescentsteepness"]:.1f}</maxdescentsteepness>\n'
        f'\t\t<maxascentsteepness>{metrics["maxascentsteepness"]:.1f}</maxascentsteepness>\n'
        f'\t\t<totalascent>{metrics["totalascent"]:.1f}</totalascent>\n'
        f'\t\t<totaldescent>{metrics["totaldescent"]:.1f}</totaldescent>\n'
        f'\t\t<maxaltitude>{metrics["maxaltitude"]:.1f}</maxaltitude>\n'
        f'\t\t<minaltitude>{metrics["minaltitude"]:.1f}</minaltitude>\n'
        f'\t\t<distance>{metrics["distance"]:.1f}</distance>\n'
        f'\t\t<descentdistance>{metrics["descentdistance"]:.1f}</descentdistance>\n'
        f'\t\t<ascentdistance>{metrics["ascentdistance"]:.1f}</ascentdistance>\n'
        f'\t\t<averagespeed>{metrics["averagespeed"]:.2f}</averagespeed>\n'
        f'\t\t<averagedescentspeed>{metrics["averagedescentspeed"]:.2f}</averagedescentspeed>\n'
        f'\t\t<averageascentspeed>{metrics["averageascentspeed"]:.2f}</averageascentspeed>\n'
        f'\t\t<movingaveragespeed>{metrics["movingaveragespeed"]:.2f}</movingaveragespeed>\n'
        f'\t\t<duration>{metrics["duration"]:.3f}</duration>\n'
        f'\t\t<startaltitude>{metrics["startaltitude"]:.1f}</startaltitude>\n'
        f'\t\t<finishaltitude>{metrics["finishaltitude"]:.1f}</finishaltitude>\n'
        f'\t\t<ascents>{int(metrics["ascents"])}</ascents>\n'
        f'\t\t<descents>{int(metrics["descents"])}</descents>\n'
        '\t</metrics>\n'
        '</track>\n'
    )

def build_from_gpx(path: str, outpath: str, tz_str: str):
    tzinfo = tzinfo_from_offset(tz_str)
    name, desc, pts = parse_gpx(path)
    course, speed = compute_course_speed(pts)

    # RawLocations.csv (9 cols): timestamp,lat,lon,alt,course,speed,hacc,vacc,timestamp2
    raw_lines=[]
    for i,(ts,la,lo,el) in enumerate(pts):
        # sane accuracies (close to Slopes)
        hacc=5.9
        vacc=14.0
        raw_lines.append(f"{ts:.6f},{la:.6f},{lo:.6f},{el:.2f},{course[i]:.1f},{speed[i]:.6f},{hacc:.1f},{vacc:.1f},{ts:.6f}")
    raw_csv="\n".join(raw_lines)+"\n"

    # Nodes.csv (8 cols): epoch(.000),lat,lon,alt,course,speed,vdop,0.25
    node_lines=[]
    prev=None
    for i,(ts,la,lo,el) in enumerate(pts):
        t_int = safe_epoch_int(ts, prev)
        prev=t_int
        node_lines.append(f"{t_int}.000,{la:.6f},{lo:.6f},{el:.2f},{course[i]:.1f},{speed[i]:.6f},0.0,0.25")
    nodes_csv="\n".join(node_lines)+"\n"

    # RAS
    ras_csv="\n".join(build_ras(pts))+"\n"

    # Segments heuristic (for initial display; refresh will rebuild)
    segs = build_segments_from_alt(pts, min_seg_s=90.0, vr_thr=0.12, target=24)

    seg_lines=['4']
    asc_id=0; des_id=0
    max_asc_sp=0.0; max_des_sp=0.0
    max_asc_grade=0.0; max_des_grade=0.0
    asc_dist=0.0; des_dist=0.0
    asc_time=0.0; des_time=0.0

    for (i0,i1,kind) in segs:
        st = seg_stats(pts, course, speed, i0, i1)
        if kind == 'ascent':
            seg_type='10'; subtype='ski_lift'
            asc_id += 1; sid=asc_id
            label=f"Remonte-pente {asc_id}"
            asc_dist += st['dist']; asc_time += st['dur']
            max_asc_sp = max(max_asc_sp, st['max_speed'])
            max_asc_grade = max(max_asc_grade, st['max_grade'])
        else:
            seg_type='8'; subtype='ski_run'
            des_id += 1; sid=des_id
            label=f"Descente à ski {des_id}"
            des_dist += st['dist']; des_time += st['dur']
            max_des_sp = max(max_des_sp, st['max_speed'])
            max_des_grade = max(max_des_grade, st['max_grade'])

        row = [
            f"{pts[i0][0]:.6f}",
            f"{pts[i1][0]:.6f}",
            seg_type,
            '0',
            str(sid),
            label,
            '',
            subtype,
            '',
            '',
            '',
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
        seg_lines.append(','.join(row))
    seg_csv="\n".join(seg_lines)+"\n"

    # Metrics for Track.xml
    start_ts=pts[0][0]; end_ts=pts[-1][0]
    total_dist=0.0
    total_ascent=0.0; total_descent=0.0
    min_alt=min(p[3] for p in pts); max_alt=max(p[3] for p in pts)
    for i in range(1,len(pts)):
        d=haversine(pts[i-1][1], pts[i-1][2], pts[i][1], pts[i][2])
        total_dist += d
        de=pts[i][3]-pts[i-1][3]
        if de>0: total_ascent += de
        else: total_descent += -de

    total_dur = max(0.0, end_ts - start_ts)
    avg_speed = (total_dist/total_dur) if total_dur > 0 else 0.0
    moving_time = 0.0
    for i in range(1,len(pts)):
        dt = max(0.0, pts[i][0]-pts[i-1][0])
        if speed[i] > 0.5:
            moving_time += dt
    moving_avg = (total_dist/moving_time) if moving_time > 0 else 0.0

    avg_asc = (asc_dist/asc_time) if asc_time>0 else 0.0
    avg_des = (des_dist/des_time) if des_time>0 else 0.0
    max_sp = max(speed) if speed else 0.0

    metrics = dict(
        start_lat=pts[0][1], start_lon=pts[0][2],
        maxspeed=max_sp,
        maxdescentspeed=max_des_sp,
        maxascentspeed=max_asc_sp,
        maxdescentsteepness=max_des_grade,
        maxascentsteepness=max_asc_grade,
        totalascent=total_ascent,
        totaldescent=total_descent,
        maxaltitude=max_alt,
        minaltitude=min_alt,
        distance=total_dist,
        descentdistance=des_dist,
        ascentdistance=asc_dist,
        averagespeed=avg_speed,
        averagedescentspeed=avg_des,
        averageascentspeed=avg_asc,
        movingaveragespeed=moving_avg,
        duration=moving_time,
        startaltitude=pts[0][3],
        finishaltitude=pts[-1][3],
        ascents=asc_id,
        descents=des_id,
    )

    track_xml = make_track_xml(name, desc, tz_str, tzinfo, start_ts, end_ts, metrics)
    events_xml = make_events_xml(start_ts, end_ts)

    files = {
        'Nodes.csv': nodes_csv,
        'RawLocations.csv': raw_csv,
        'RelativeAltitudeSensor.csv': ras_csv,
        'Segment.csv': seg_csv,
        'Events.xml': events_xml,
        'Samples.xml': SAMPLES_XML,
        'Photos.csv': '',
        'Track.xml': track_xml,
    }
    write_skiz(files, outpath)

def build_from_slopes(path: str, outpath: str):
    name, desc, tz_str, tzinfo, start_ts, end_ts, gps_rows, raw_rows, segs, attrs = parse_slopes(path)

    # Nodes.csv (repo style)
    node_lines=[]
    prev=None
    for r in gps_rows:
        ts, la, lo, alt, course, spd, _, _ = r
        t_int = safe_epoch_int(ts, prev)
        prev=t_int
        node_lines.append(f"{t_int}.000,{la:.6f},{lo:.6f},{alt:.2f},{course:.1f},{spd:.6f},0.0,0.25")
    nodes_csv="\n".join(node_lines)+"\n"

    # RawLocations.csv: 9 cols (add ts2)
    raw_lines=[]
    for r in raw_rows:
        ts, la, lo, alt, course, spd, hacc, vacc = r
        raw_lines.append(f"{ts:.6f},{la:.6f},{lo:.6f},{alt:.3f},{course:.2f},{spd:.6f},{hacc:.3f},{vacc:.3f},{ts:.6f}")
    raw_csv="\n".join(raw_lines)+"\n"

    # RAS from altitude in raw_rows
    pts=[(r[0], r[1], r[2], r[3]) for r in raw_rows]
    ras_csv="\n".join(build_ras(pts))+"\n"

    # Segment.csv from actions
    seg_lines=['4']
    asc_id=0; des_id=0
    max_asc_sp=0.0; max_des_sp=0.0
    max_asc_grade=0.0; max_des_grade=0.0
    asc_dist=0.0; des_dist=0.0
    asc_time=0.0; des_time=0.0

    for st,en,typ,a in segs:
        if typ.lower() == 'lift':
            seg_type='10'; subtype='ski_lift'
            asc_id += 1; sid=asc_id
            label=f"Remonte-pente {asc_id}"
        else:
            seg_type='8'; subtype='ski_run'
            des_id += 1; sid=des_id
            label=f"Descente à ski {des_id}"

        dur=float(a.get('duration','0'))
        dist=float(a.get('distance','0'))
        avg_sp=float(a.get('avgSpeed','0'))
        top_sp=float(a.get('topSpeed','0'))
        vertical=float(a.get('vertical','0'))
        max_alt=float(a.get('maxAlt','0'))
        min_alt=float(a.get('minAlt','0'))

        # steepness approximations
        avg_grade = (abs(vertical)/dist*100.0) if dist>0 else 0.0
        max_grade = avg_grade

        if seg_type=='10':
            asc_dist += dist; asc_time += dur
            max_asc_sp=max(max_asc_sp, top_sp)
            max_asc_grade=max(max_asc_grade, max_grade)
        else:
            des_dist += dist; des_time += dur
            max_des_sp=max(max_des_sp, top_sp)
            max_des_grade=max(max_des_grade, max_grade)

        row = [
            f"{st:.6f}",
            f"{en:.6f}",
            seg_type,
            '0',
            str(sid),
            label,
            '',
            subtype,
            '',
            '',
            '',
            f"{dur:.3f}",
            f"{avg_sp:.6f}",
            f"{dist:.6f}",
            f"{vertical:.6f}",
            f"{top_sp:.6f}",
            f"{avg_grade:.6f}",
            f"{max_grade:.6f}",
            f"{min_alt:.2f}",
            f"{max_alt:.2f}",
            f"{min_alt:.2f}",
            f"{max_alt:.2f}",
        ]
        seg_lines.append(','.join(row))
    seg_csv="\n".join(seg_lines)+"\n"

    # Metrics largely from attrs + fallbacks from raw
    try:
        total_dist=float(attrs.get('distance','0'))
    except Exception:
        total_dist=0.0
    try:
        total_dur=float(attrs.get('duration','0'))
    except Exception:
        total_dur=max(0.0, end_ts-start_ts)

    # altitude stats from pts
    min_alt=min(p[3] for p in pts) if pts else 0.0
    max_alt=max(p[3] for p in pts) if pts else 0.0

    # vertical totals from attrs if present
    total_descent=float(attrs.get('vertical','0')) if attrs.get('vertical') else 0.0
    # crude ascent from actions if missing
    total_ascent=0.0
    for _,_,typ,a in segs:
        if typ.lower()=='lift':
            try: total_ascent += float(a.get('vertical','0'))
            except Exception: pass

    # speed max from attrs
    max_speed=float(attrs.get('topSpeed','0')) if attrs.get('topSpeed') else 0.0

    avg_speed=(total_dist/total_dur) if total_dur>0 else 0.0
    avg_asc=(asc_dist/asc_time) if asc_time>0 else 0.0
    avg_des=(des_dist/des_time) if des_time>0 else 0.0

    moving_avg=avg_speed
    moving_time=total_dur

    metrics = dict(
        start_lat=float(attrs.get('centerLat','0')) if attrs.get('centerLat') else (pts[0][1] if pts else 0.0),
        start_lon=float(attrs.get('centerLong','0')) if attrs.get('centerLong') else (pts[0][2] if pts else 0.0),
        maxspeed=max_speed,
        maxdescentspeed=max_des_sp,
        maxascentspeed=max_asc_sp,
        maxdescentsteepness=max_des_grade,
        maxascentsteepness=max_asc_grade,
        totalascent=total_ascent,
        totaldescent=total_descent,
        maxaltitude=max_alt,
        minaltitude=min_alt,
        distance=total_dist,
        descentdistance=des_dist,
        ascentdistance=asc_dist,
        averagespeed=avg_speed,
        averagedescentspeed=avg_des,
        averageascentspeed=avg_asc,
        movingaveragespeed=moving_avg,
        duration=moving_time,
        startaltitude=min_alt,
        finishaltitude=min_alt,
        ascents=asc_id,
        descents=des_id,
    )

    track_xml = make_track_xml(name, desc, tz_str, tzinfo, start_ts, end_ts, metrics)
    events_xml = make_events_xml(start_ts, end_ts)

    files = {
        'Nodes.csv': nodes_csv,
        'RawLocations.csv': raw_csv,
        'RelativeAltitudeSensor.csv': ras_csv,
        'Segment.csv': seg_csv,
        'Events.xml': events_xml,
        'Samples.xml': SAMPLES_XML,
        'Photos.csv': '',
        'Track.xml': track_xml,
    }
    write_skiz(files, outpath)

# --------------------
# Dispatch by extension
# --------------------
lower = INPATH.lower()
if lower.endswith('.slopes'):
    build_from_slopes(INPATH, OUTSKIZ)
elif lower.endswith('.gpx'):
    build_from_gpx(INPATH, OUTSKIZ, DEFAULT_TZ)
else:
    raise SystemExit('Unsupported input extension (need .gpx or .slopes)')
PY
  then
    log ERROR "Failed: $inpath (see $per_log)"
    return 1
  fi

  if unzip -tq "$out_skiz" >>"$per_log" 2>&1; then
    log INFO "OK: $out_skiz"
    return 0
  fi

  log ERROR "Zip integrity check failed: $out_skiz (see $per_log)"
  return 1
}

fail=0
for f in "${FILES[@]}"; do
  convert_one "$f" || fail=1
done

if [[ "$fail" -ne 0 ]]; then
  log ERROR "Some conversions failed. Check logs in $LOG_DIR"
  exit 1
fi

log INFO "Done. Output in $OUT_DIR"
