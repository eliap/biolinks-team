#!/usr/bin/env python3
"""
Build biolinks_data.js for the Upper Wimmera Landcare groups.

Python port of biolinks_data_builder.html (same thresholds, same output schema) that
works over the whole Landcare area instead of a small test square.

  .venv\\Scripts\\python build_biolinks_data.py                    # whole Upper Wimmera (catchment + all groups)
  .venv\\Scripts\\python build_biolinks_data.py --groups-only       # only land inside the Landcare groups
  .venv\\Scripts\\python build_biolinks_data.py --group "Stawell Urban"
  .venv\\Scripts\\python build_biolinks_data.py --square=-37.214978,142.797641,15  # same area as the browser builder
  .venv\\Scripts\\python build_biolinks_data.py --refresh-inat      # re-download iNaturalist only
  .venv\\Scripts\\python build_biolinks_data.py --refresh           # re-download everything

Downloads are cached in cache/ (one file per page), so an interrupted run carries on
where it stopped and re-runs only re-process.
"""
import argparse
import csv
import gzip
import json
import math
import os
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
import shapefile
import shapely
from pyproj import Transformer
from shapely import STRtree
from shapely.geometry import MultiPolygon, Polygon, box, shape
from shapely.ops import substring, unary_union

HERE = Path(__file__).resolve().parent
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

# ---------------- Settings (match biolinks_data_builder.html) ----------------
WFS = 'https://opendata.maps.vic.gov.au/geoserver/wfs'
WFS_PAGE = 5000
# layer -> (sort key needed for paging, fields to download)
WFS_LAYERS = {
    'v_property_mp': ('prop_pfi', ['prop_pfi', 'prop_propnum', 'geom']),
    'v_parcel_mp': ('parcel_pfi', ['parcel_pfi', 'parcel_spi', 'parcel_road', 'geom']),
    'address': ('pfi', ['pfi', 'property_pfi', 'ezi_address']),
    'nv2005_evcbcs': (None, ['evc', 'x_evcname', 'evc_bcs_desc', 'geom']),
    'road_casement_polygon': ('pfi', ['pfi', 'geom']),
    'tr_road': ('ufi', ['ufi', 'ezi_road_name_label', 'geom']),
}
INAT = 'https://api.inaturalist.org/v1/observations'
USER_AGENT = 'biolinks-data-builder/1.0 (Project Platypus / Upper Wimmera Landcare habitat map; python-requests)'
MAX_ACCURACY_M = 200
ROADSIDE_SECTION_M = 500
ROAD_BUFFER_M = 40
MIN_SECTION_M2 = 200          # pieces smaller than this are dropped
LEFTOVER_IGNORE_M2 = 30       # leftover casement slivers smaller than this are ignored
LEFTOVER_MERGE_M2 = 3000      # leftovers smaller than this merge into a touching section
LEFTOVER_KEEP_M2 = 500        # larger leftovers become "unformed" sections
UNFORMED = 'Road reserve with no mapped road (unformed)'
VEG_SIMPLIFY_M = 1.0          # NV2005 patches are blocky raster shapes
OTHER_SIMPLIFY_M = 0.5        # removes redundant points only (output is rounded to ~1 m anyway)
GROUP_SIMPLIFY_M = 5.0

LIFEFORMS = {"TREE": "Canopy tree", "UTLS": "Understorey tree or large shrub", "MS": "Medium shrub", "SS": "Small or prostrate shrub",
             "LH": "Large herb", "MH": "Medium herb", "SPH": "Small or prostrate herb", "LTG": "Large tufted graminoid",
             "MTG": "Medium to tiny tufted graminoid", "NTG": "Non-tufted graminoid", "GF": "Ground fern",
             "SCC": "Scrambler or climber", "BRY": "Bryophyte or liverwort", "AQ": "Aquatic herb", "MIS": "Mistletoe"}
EXCLUDED = ["TREE", "MIS"]

# EPSG:3111 = VicGrid94 (metres) for all geometry work
LL_TO_M = Transformer.from_crs(4326, 3111, always_xy=True)
M_TO_LL = Transformer.from_crs(3111, 4326, always_xy=True)
MGA_TO_M = Transformer.from_crs(28354, 3111, always_xy=True)

# ---------------- Logging ----------------
T0 = time.time()
_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


class Stage:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.time()
        log(f'== {self.name}')
        return self

    def __exit__(self, *a):
        log(f'   {self.name} done in {time.time() - self.t:.1f} s')


# ---------------- Geometry helpers ----------------
def reproject(geoms, transformer):
    def f(c):
        x, y = transformer.transform(c[:, 0], c[:, 1])
        return np.column_stack([x, y])
    return shapely.transform(geoms, f)


def polyonly(g):
    """Keep only the polygon parts of a geometry (intersections can return lines/points)."""
    if g is None or g.is_empty:
        return None
    t = g.geom_type
    if t in ('Polygon', 'MultiPolygon'):
        return g
    if t == 'GeometryCollection':
        ps = []
        for p in shapely.get_parts(g):
            q = polyonly(p)
            if q is not None:
                ps.extend(shapely.get_parts(q))
        return MultiPolygon(ps) if ps else None
    return None


def clean_polys(arr):
    arr = np.asarray(arr, dtype=object)
    bad = ~shapely.is_valid(arr)
    if bad.any():
        arr[bad] = shapely.make_valid(arr[bad])
        for i in np.nonzero(bad)[0]:
            arr[i] = polyonly(arr[i])
    return arr


def safe(op, a, b, fallback=None):
    try:
        return op(a, b)
    except shapely.errors.GEOSException:
        try:
            return op(shapely.make_valid(a), shapely.make_valid(b))
        except shapely.errors.GEOSException:
            return fallback


def bbox_overlap(a, b):
    return not (a[0] > b[2] or a[2] < b[0] or a[1] > b[3] or a[3] < b[1])


def encode_geom(g):
    """Same compact encoding as the browser builder: per ring, round(lon*1e5), round(lat*1e5),
    first pair absolute then deltas, consecutive duplicates skipped, rings < 4 points dropped,
    polygons whose outer ring was dropped are dropped."""
    out = []
    if g is None or g.is_empty:
        return out
    for poly in shapely.get_parts(g):
        if poly.geom_type != 'Polygon':
            continue
        rings = []
        for k, ring in enumerate([poly.exterior, *poly.interiors]):
            c = shapely.get_coordinates(ring)
            ix = np.floor(c[:, 0] * 1e5 + 0.5).astype(np.int64)   # JS Math.round
            iy = np.floor(c[:, 1] * 1e5 + 0.5).astype(np.int64)
            keep = np.ones(len(ix), bool)
            keep[1:] = (ix[1:] != ix[:-1]) | (iy[1:] != iy[:-1])
            ix, iy = ix[keep], iy[keep]
            if len(ix) < 4:
                if k == 0:
                    break          # outer ring dropped -> drop the polygon
                continue
            d = np.empty(2 * len(ix), np.int64)
            d[0::2] = np.diff(ix, prepend=0)
            d[1::2] = np.diff(iy, prepend=0)
            rings.append(d.tolist())
        else:
            out.append(rings)
    return out


def encode_m(geoms_m, tol):
    """Simplify (metres), reproject to lon/lat and encode a list of geometries."""
    arr = np.asarray(geoms_m, dtype=object)
    if tol:
        arr = shapely.simplify(arr, tol, preserve_topology=True)
    arr = reproject(arr, M_TO_LL)
    return [encode_geom(g) for g in arr]


# ---------------- HTTP ----------------
session = requests.Session()
session.headers['User-Agent'] = USER_AGENT


def get_json(url, params, tries=8, label=''):
    for i in range(tries):
        try:
            r = session.get(url, params=params, timeout=(30, 90 if "inaturalist" in url else 600))
            if r.status_code == 429 or r.status_code >= 500:
                wait = float(r.headers.get('Retry-After') or 0)
                raise RuntimeError(f'HTTP {r.status_code}', wait)
            if not r.ok:
                raise RuntimeError(f'HTTP {r.status_code}: {r.text[:200]}')
            try:
                return r.json()
            except ValueError:
                txt = ' '.join(r.text.replace('<', ' <').split())
                raise RuntimeError('server error: ' + txt[:240])
        except (requests.RequestException, RuntimeError) as e:
            if i == tries - 1:
                raise
            wait = max(min(5 * 2 ** i, 120), e.args[1] if isinstance(e, RuntimeError) and len(e.args) > 1 else 0)
            log(f'  {label} retrying in {wait:.0f} s after error ({e.args[0] if e.args else e})')
            time.sleep(wait)


def bbox_key(b):
    return '_'.join(f'{v:.4f}' for v in b)


def covers(outer, inner):
    return outer[0] <= inner[0] + 1e-9 and outer[1] <= inner[1] + 1e-9 and outer[2] >= inner[2] - 1e-9 and outer[3] >= inner[3] - 1e-9


def write_gz(path, obj):
    tmp = path.with_suffix(path.suffix + '.tmp')
    with gzip.open(tmp, 'wt', encoding='utf-8') as f:
        json.dump(obj, f, separators=(',', ':'))
    os.replace(tmp, path)


def read_gz(path):
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        return json.load(f)


def find_cache(kind_dir, prefix, bbox):
    """A complete cached download whose bbox covers the one we need (e.g. a group inside the full area)."""
    if not kind_dir.exists():
        return None
    best = None
    for d in kind_dir.glob(prefix + '_*'):
        done = d / 'complete.json'
        if done.exists():
            meta = json.loads(done.read_text())
            if covers(meta['bbox'], bbox):
                area = (meta['bbox'][2] - meta['bbox'][0]) * (meta['bbox'][3] - meta['bbox'][1])
                if best is None or area < best[0]:
                    best = (area, d, meta)
    return best


# ---------------- Downloads ----------------
def download_wfs(cache, layer, bbox, refresh):
    sort_key, fields = WFS_LAYERS[layer]
    wdir = cache / 'wfs'
    if not refresh:
        hit = find_cache(wdir, layer, bbox)
        if hit:
            log(f'{layer}: using cache ({hit[2]["count"]:,} features, downloaded {hit[2]["downloaded"][:10]})')
            return load_pages(hit[1])
    d = wdir / f'{layer}_{bbox_key(bbox)}'
    if refresh and d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)
    params = {'service': 'WFS', 'version': '2.0.0', 'request': 'GetFeature', 'typeNames': 'open-data-platform:' + layer,
              'outputFormat': 'application/json', 'srsName': 'EPSG:4326', 'bbox': ','.join(map(str, bbox)) + ',EPSG:4326',
              'count': WFS_PAGE, 'propertyName': ','.join(fields)}
    if sort_key:
        params['sortBy'] = sort_key
    start, total = 0, 0
    while True:
        page = d / f'page_{start:07d}.json.gz'
        if page.exists():
            n = len(read_gz(page)['features'])
        else:
            j = get_json(WFS, dict(params, startIndex=start), label=layer)
            feats = [{'id': f.get('id'), 'p': f.get('properties') or {}, 'g': f.get('geometry')} for f in j.get('features', [])]
            write_gz(page, {'features': feats})
            n = len(feats)
            log(f'  {layer}: {start + n:,} downloaded…')
        total += n
        if n < WFS_PAGE:
            break
        start += WFS_PAGE
    (d / 'complete.json').write_text(json.dumps({'bbox': bbox, 'count': total, 'layer': layer,
                                                 'downloaded': datetime.now(timezone.utc).isoformat()}))
    log(f'{layer}: {total:,} features downloaded')
    return load_pages(d)


def load_pages(d):
    seen, out = set(), []
    for page in sorted(d.glob('page_*.json.gz')):
        for f in read_gz(page)['features']:
            if f['id'] in seen:
                continue
            seen.add(f['id'])
            out.append(f)
    return out


def trim_obs(o):
    t = o.get('taxon')
    g = o.get('geojson')
    return {'id': o['id'], 'c': g.get('coordinates') if g else None, 'acc': o.get('positional_accuracy'),
            'd': o.get('observed_on') or '', 'q': o.get('quality_grade'),
            't': {k: t.get(k) for k in ('id', 'name', 'rank_level', 'min_species_taxon_id', 'preferred_common_name', 'introduced')} if t else None}


def download_inat(cache, bbox, refresh):
    idir = cache / 'inat'
    if not refresh:
        hit = find_cache(idir, 'inat', bbox)
        if hit:
            log(f'iNaturalist: using cache ({hit[2]["count"]:,} records, downloaded {hit[2]["downloaded"][:10]}; --refresh-inat for new records)')
            return load_inat(hit[1])
    d = idir / f'inat_{bbox_key(bbox)}'
    if refresh and d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)
    W, S, E, N = bbox
    pages = sorted(d.glob('page_*.json.gz'))
    id_above, n = 0, 0
    if pages:   # resume an interrupted download
        for p in pages:
            recs = read_gz(p)['results']
            n += len(recs)
            if recs:
                id_above = max(id_above, recs[-1]['id'])
        log(f'iNaturalist: resuming after {n:,} cached records')
    page_no = len(pages)
    last = 0.0
    while True:
        params = {'iconic_taxa': 'Plantae', 'verifiable': 'true', 'geoprivacy': 'open', 'taxon_geoprivacy': 'open',
                  'swlat': S, 'swlng': W, 'nelat': N, 'nelng': E, 'preferred_place_id': 7830, 'locale': 'en',
                  'per_page': 200, 'order_by': 'id', 'order': 'asc', 'id_above': id_above}
        wait = 1.0 - (time.time() - last)        # about one request per second
        if wait > 0:
            time.sleep(wait)
        last = time.time()
        j = get_json(INAT, params, label='iNaturalist')
        res = [trim_obs(o) for o in j.get('results', [])]
        if page_no == 0 or page_no % 20 == 0:
            log(f'  iNaturalist: {n:,} of ~{n + j.get("total_results", 0):,} records…')
        if res:
            write_gz(d / f'page_{page_no:05d}.json.gz', {'results': res})
            page_no += 1
            n += len(res)
            id_above = res[-1]['id']
        if len(res) < 200:
            break
    (d / 'complete.json').write_text(json.dumps({'bbox': bbox, 'count': n, 'downloaded': datetime.now(timezone.utc).isoformat()}))
    log(f'iNaturalist: {n:,} plant records downloaded')
    return load_inat(d)


def load_inat(d):
    seen, out = set(), []
    for p in sorted(d.glob('page_*.json.gz')):
        for o in read_gz(p)['results']:
            if o['id'] not in seen:
                seen.add(o['id'])
                out.append(o)
    return out


# ---------------- Lookup ----------------
def origin_code(text):
    t = (text or '').lower().strip()
    if t.startswith('introduced') or t == 'i':
        return 'I'
    if 'not local' in t or 'planted' in t or t == 'p':
        return 'P'
    return 'N'


def load_lookup(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
    head = [h.strip().lower() for h in rows[0]]
    need = ['taxon_id', 'scientific_name', 'origin', 'lifeform_code', 'counts_as_native_understorey']
    missing = [n for n in need if n not in head]
    if missing:
        sys.exit(f'{path.name}: missing column(s): {", ".join(missing)}')
    col = {h: i for i, h in enumerate(head)}
    taxa = {}
    for r in rows[1:]:
        r = r + [''] * (len(head) - len(r))
        try:
            tid = int(r[col['taxon_id']].strip())
        except ValueError:
            continue
        lf = r[col['lifeform_code']].strip().upper()
        if lf and lf not in LIFEFORMS:
            sys.exit(f'{path.name}: unknown lifeform code "{lf}" for {r[col["scientific_name"]]}')
        taxa[tid] = [r[col['scientific_name']].strip(), lf, origin_code(r[col['origin']]),
                     'Y' if r[col['counts_as_native_understorey']].strip().upper() == 'Y' else 'N']
    seen = {}
    for sci, lf, org, _ in taxa.values():
        seen.setdefault(sci.split(' ')[0], set()).add((lf, org))
    genus = {g: list(next(iter(s))) for g, s in seen.items() if len(s) == 1}
    return taxa, genus


def classify(lookup, tid, name, inat_introduced):
    taxa, genus = lookup
    L = taxa.get(tid)
    if L:
        return L[1] or '', L[2], L[3] == 'Y' and L[1] not in EXCLUDED, 'lookup'
    G = genus.get((name or '').split(' ')[0])
    if G:
        return G[0], G[1], G[1] == 'N' and G[0] not in EXCLUDED, 'genus'
    origin = 'I' if inat_introduced else 'N'
    return '', origin, origin == 'N', 'unclassified'


# ---------------- Area ----------------
def load_groups(path):
    sf = shapefile.Reader(str(path))
    fields = [f[0] for f in sf.fields[1:]]
    order, geoms, abbr = [], {}, {}
    for sr in sf.iterShapeRecords():
        rec = dict(zip(fields, sr.record))
        name = (rec.get('GROUP_NAME') or rec.get('Group') or '').strip()
        g = shapely.make_valid(shape(sr.shape.__geo_interface__))
        if name not in geoms:
            order.append(name)
            geoms[name] = []
            abbr[name] = (rec.get('Group_Abbr') or '').strip()
        geoms[name].append(g)
    out = []
    for name in order:   # dissolve by group name (Crowlands has two polygons)
        g = polyonly(unary_union(reproject(np.array(geoms[name], dtype=object), MGA_TO_M)))
        out.append({'name': name, 'abbr': abbr[name], 'geom': g})
    return out


def load_catchment(path, region):
    """Upper Wimmera outline: the 'Upper Catchment' polygon of the 2007 Wimmera Landcare network boundaries."""
    sf = shapefile.Reader(str(path))
    fields = [f[0] for f in sf.fields[1:]]
    parts = []
    for sr in sf.iterShapeRecords():
        rec = dict(zip(fields, sr.record))
        if region.lower() in (str(rec.get('LCRegion', '')).lower(), str(rec.get('Network', '')).lower()):
            parts.append(shapely.make_valid(shape(sr.shape.__geo_interface__)))
    if not parts:
        sys.exit(f'{Path(path).name}: no polygon with LCRegion/Network "{region}"')
    return polyonly(unary_union(reproject(np.array(parts, dtype=object), MGA_TO_M)))


def square_bbox(lat, lng, size):
    half = size / 2
    dlat = half / 110.574
    dlon = half / (111.320 * math.cos(lat * math.pi / 180))
    return [lng - dlon, lat - dlat, lng + dlon, lat + dlat]


# ---------------- Roadside sections (port of buildRoadsides) ----------------
class Grid:
    def __init__(self, size=1000.0):
        self.size, self.cells = size, {}

    def _keys(self, b):
        s = self.size
        for x in range(int(math.floor(b[0] / s)), int(math.floor(b[2] / s)) + 1):
            for y in range(int(math.floor(b[1] / s)), int(math.floor(b[3] / s)) + 1):
                yield (x, y)

    def add(self, i, b):
        for k in self._keys(b):
            lst = self.cells.setdefault(k, [])
            if not lst or lst[-1] != i:
                lst.append(i)

    def near(self, b):
        s = set()
        for k in self._keys(b):
            s.update(self.cells.get(k, ()))
        return sorted(s)


def build_roadsides(cas, roads):
    """cas: list of casement polygons (m); roads: list of (name, line geometry m)."""
    cells = []   # dicts: geom, bounds, road, lengthM
    grid = Grid()
    cas_tree = STRtree(cas)

    def add_piece(piece, road, length_m):
        pb = piece.bounds
        prev = [cells[i]['geom'] for i in grid.near(pb) if bbox_overlap(pb, cells[i]['bounds'])]
        if prev:
            prev = [p for p in prev if p.intersects(piece)]
        for q in prev:                       # sections never overlap, so subtract them one by one
            piece = polyonly(safe(shapely.difference, piece, q, piece))
            if piece is None:
                return
        if piece.area > MIN_SECTION_M2:
            c = {'geom': piece, 'bounds': piece.bounds, 'road': road, 'lengthM': length_m}
            cells.append(c)
            grid.add(len(cells) - 1, c['bounds'])

    t = time.time()
    nchunks = 0
    for ri, (name, geom) in enumerate(roads):
        if ri and ri % 2000 == 0:
            log(f'  roads {ri:,}/{len(roads):,} · {len(cells):,} sections · {time.time() - t:.0f} s')
        for ln in shapely.get_parts(geom):
            L = ln.length
            if L <= 0:
                continue
            s = 0.0
            while s < L - 1e-9:
                chunk = substring(ln, s, min(s + ROADSIDE_SECTION_M, L)) if L > ROADSIDE_SECTION_M else ln
                s += ROADSIDE_SECTION_M
                nchunks += 1
                length_m = chunk.length
                buf = chunk.buffer(ROAD_BUFFER_M)
                for ci in sorted(cas_tree.query(buf, predicate='intersects')):
                    piece = polyonly(safe(shapely.intersection, buf, cas[ci]))
                    if piece is not None:
                        add_piece(piece, name, length_m)
    log(f'  {nchunks:,} centreline chunks → {len(cells):,} sections; now sorting out leftover casement')

    merged = unformed = 0
    for c in cas:
        cb = c.bounds
        prev = [cells[i]['geom'] for i in grid.near(cb) if bbox_overlap(cb, cells[i]['bounds'])]
        prev = [p for p in prev if p.intersects(c)]
        rest = c
        for q in prev:
            rest = polyonly(safe(shapely.difference, rest, q, rest))
            if rest is None:
                break
        if rest is None:
            continue
        for part in shapely.get_parts(rest):
            a = part.area
            if a < LEFTOVER_IGNORE_M2:
                continue
            touching = None
            if a < LEFTOVER_MERGE_M2:
                grown = part.buffer(3)
                gb = grown.bounds
                for i in grid.near(gb):
                    if bbox_overlap(gb, cells[i]['bounds']) and grown.intersects(cells[i]['geom']):
                        touching = i
                        break
            if touching is not None:
                u = polyonly(safe(shapely.union, cells[touching]['geom'], part))
                if u is not None:
                    cells[touching]['geom'] = u
                    cells[touching]['bounds'] = u.bounds
                    grid.add(touching, u.bounds)
                    merged += 1
            elif a > LEFTOVER_KEEP_M2:
                cells.append({'geom': part, 'bounds': part.bounds, 'road': UNFORMED, 'lengthM': None})
                grid.add(len(cells) - 1, part.bounds)
                unformed += 1
    log(f'  leftovers: {merged:,} slivers merged into sections, {unformed:,} unformed road reserve sections')
    return cells


# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--group', action='append', help='only build this Landcare group (name or part of it); repeatable')
    ap.add_argument('--square', help='--square=lat,lng,km (note the =): build a square like the browser builder instead of the group area')
    ap.add_argument('--refresh', action='store_true', help='re-download everything')
    ap.add_argument('--refresh-inat', action='store_true', help='re-download iNaturalist records only')
    ap.add_argument('--groups-only', action='store_true', help='only build land inside the Landcare group polygons')
    ap.add_argument('--boundaries', default=str(HERE / 'landcare' / 'All Upper Wimmera Landcares.shp'))
    ap.add_argument('--catchment', default=str(HERE / 'landcare' / 'network data' / 'wca_LandcareNetworks_2007_mga54.shp'),
                    help='shapefile (MGA54) holding the Upper Wimmera outline')
    ap.add_argument('--catchment-region', default='Upper Catchment', help='LCRegion/Network value of the outline in --catchment')
    ap.add_argument('--lookup', default=str(HERE / 'understorey_lookup.csv'))
    ap.add_argument('--cache', default=str(HERE / 'cache'))
    ap.add_argument('--out', default=str(HERE / 'biolinks_data.js'))
    ap.add_argument('--unclassified', default=str(HERE / 'unclassified_plants.csv'))
    ap.add_argument('--check-parcel', action='append', default=[], help=r'print plant counts for a parcel SPI, e.g. S6\PP3207')
    args = ap.parse_args()
    cache = Path(args.cache)

    # ---- 0. area & lookup ----
    with Stage('Area and lookup'):
        lookup_path = Path(args.lookup)
        lookup = load_lookup(lookup_path)
        lookup_text = f'{lookup_path.name} ({len(lookup[0])} taxa)'
        log(f'Lookup: {lookup_text}, {len(lookup[1])} genera with a single lifeform/origin')
        groups_all = load_groups(Path(args.boundaries))
        log(f'Landcare groups: {len(groups_all)} after dissolving by name')
        if args.square:
            lat, lng, size = map(float, args.square.split(','))
            bbox = square_bbox(lat, lng, size)
            area = reproject(shapely.segmentize(box(*bbox), 0.001), LL_TO_M)
            catchment = None
            groups = [g for g in groups_all if g['geom'].intersects(area)]
            mode = 'square'
            log(f'Area: {size:g} × {size:g} km square around {lat}, {lng} (overlaps {len(groups)} groups)')
        else:
            groups = groups_all
            if args.group:
                sel = []
                for q in args.group:
                    ql = q.lower()
                    m = [g for g in groups_all if g['name'].lower() == ql or g['abbr'].lower() == ql]
                    m = m or [g for g in groups_all if ql in g['name'].lower() or ql in g['abbr'].lower()]
                    if not m:
                        sys.exit(f'No group matches "{q}". Groups: ' + '; '.join(g['name'] for g in groups_all))
                    sel += [g for g in m if g not in sel]
                groups = sel
            area = unary_union([g['geom'] for g in groups])
            catchment = None
            if not args.group and not args.groups_only:
                catchment = unary_union([load_catchment(args.catchment, args.catchment_region), area])
                log(f'Upper Wimmera outline ({Path(args.catchment).name}, {args.catchment_region}) + groups: '
                    f'{catchment.area / 1e6:,.0f} km², of which {catchment.difference(area).area / 1e6:,.0f} km² is not in any group')
                area = catchment
            mode = 'catchment' if catchment is not None else 'groups'
            ll = reproject(area, M_TO_LL).bounds
            pad = 0.002
            bbox = [round(ll[0] - pad, 4), round(ll[1] - pad, 4), round(ll[2] + pad, 4), round(ll[3] + pad, 4)]
            log(f'Area: {"Upper Wimmera" if catchment is not None else ", ".join(g["name"] for g in groups)} · {area.area / 1e6:,.0f} km²')
        log(f'Bounding box (W,S,E,N): {bbox}')
        shapely.prepare(area)

    # ---- 1. downloads ----
    with Stage('Downloads'):
        with ThreadPoolExecutor(max_workers=4) as ex:
            inat_f = ex.submit(download_inat, cache, bbox, args.refresh or args.refresh_inat)
            futs = {layer: ex.submit(download_wfs, cache, layer, bbox, args.refresh) for layer in
                    ['nv2005_evcbcs', 'v_parcel_mp', 'v_property_mp', 'tr_road', 'road_casement_polygon', 'address']}
            raw = {k: f.result() for k, f in futs.items()}
            obs_raw = inat_f.result()

    # ---- 2. geometry into metres, keep what touches the area ----
    def to_m(feats, polygonal=True):
        feats = [f for f in feats if f['g']]
        geoms = np.array([shape(f['g']) for f in feats], dtype=object)
        geoms = reproject(geoms, LL_TO_M)
        if polygonal:
            geoms = clean_polys(geoms)
        keep = np.array([g is not None and not g.is_empty for g in geoms], bool)
        keep[keep] = shapely.intersects(area, geoms[keep])
        return [f for f, k in zip(feats, keep) if k], list(geoms[keep])

    with Stage('Clip to area'):
        parcel_f, parcel_g = to_m(raw.pop('v_parcel_mp'))
        prop_f, prop_g = to_m(raw.pop('v_property_mp'))
        veg_f, veg_g = to_m(raw.pop('nv2005_evcbcs'))
        cas_f, cas_g = to_m(raw.pop('road_casement_polygon'))
        road_feats = [f for f in raw.pop('tr_road') if f['g']]
        road_g = reproject(np.array([shape(f['g']) for f in road_feats], dtype=object), LL_TO_M)
        if cas_g:   # only centrelines that could reach a kept road reserve
            idx = np.unique(STRtree(cas_g).query(road_g, predicate='dwithin', distance=ROAD_BUFFER_M)[0])
        else:
            idx = np.array([], int)
        roads = [((road_feats[i]['p'].get('ezi_road_name_label') or 'Unnamed road'), road_g[i]) for i in idx]
        addr_f = raw.pop('address')
        log(f'Parcels {len(parcel_f):,} · properties {len(prop_f):,} · NV2005 patches {len(veg_f):,} · '
            f'road reserves {len(cas_f):,} · road centrelines {len(roads):,} · addresses {len(addr_f):,}')

    # ---- 3. roadside sections ----
    with Stage(f'Roadside sections ({ROADSIDE_SECTION_M} m)'):
        roadsides = build_roadsides(cas_g, roads)
        n_all = len(roadsides)
        roadsides = [c for c in roadsides if area.intersects(c['geom'])]   # road reserves run on past the boundary
        if n_all > len(roadsides):
            log(f'  {n_all - len(roadsides):,} sections outside the area dropped')
        road_g_m = [c['geom'] for c in roadsides]
        log(f'Roadside sections: {len(roadsides):,}')

    # ---- 4. place iNat records ----
    with Stage('Place iNaturalist records'):
        taxa = {}
        cand = []
        coarse = 0
        for o in obs_raw:
            if not o['c'] or not o['t']:
                continue
            if o['acc'] and o['acc'] > MAX_ACCURACY_M:
                coarse += 1
                continue
            cand.append(o)
        xy = np.array([o['c'][:2] for o in cand], float).reshape(-1, 2)
        pts_m = shapely.points(np.column_stack(LL_TO_M.transform(xy[:, 0], xy[:, 1]))) if len(cand) else np.array([], object)
        n = len(cand)
        kind = np.zeros(n, int)
        unit = np.full(n, -1, int)
        if parcel_g and n:
            pi, ui = STRtree(parcel_g).query(pts_m, predicate='intersects')
            best = np.full(n, np.iinfo(np.int64).max)
            np.minimum.at(best, pi, ui)          # first parcel in list order wins (as the builder)
            hit = best < np.iinfo(np.int64).max
            kind[hit], unit[hit] = 1, best[hit]
        if road_g_m and n:
            pi, ui = STRtree(road_g_m).query(pts_m, predicate='intersects')
            best = np.full(n, -1)
            np.maximum.at(best, pi, ui)          # builder keeps the last roadside hit
            hit = (best >= 0) & (kind == 0)
            kind[hit], unit[hit] = 2, best[hit]
        inside = shapely.intersects(area, pts_m) if n else np.array([], bool)
        keep = (kind > 0) | inside
        outside = int((~keep).sum())
        obs_out, obs_pts = [], []
        for i in np.nonzero(keep)[0]:
            o = cand[i]
            t = o['t']
            rl = t.get('rank_level') or 0
            tid = t['min_species_taxon_id'] if rl < 10 and t.get('min_species_taxon_id') else t['id']
            if tid not in taxa:
                name = ' '.join(t['name'].split(' ')[:2]) if rl < 10 else t['name']
                lf, origin, und, src = classify(lookup, tid, name, bool(t.get('introduced')))
                taxa[tid] = [name, t.get('preferred_common_name') or '', lf, origin, 1 if und else 0, src, 1 if rl <= 10 else 0]
            x, y = o['c'][:2]
            obs_out.append([round(x, 5), round(y, 5), tid, rl, o['d'], 1 if o['q'] == 'research' else 0, o['id'], int(kind[i]), int(unit[i])])
            obs_pts.append(pts_m[i])
        placed_p = sum(1 for o in obs_out if o[7] == 1)
        placed_r = sum(1 for o in obs_out if o[7] == 2)
        unplaced = sum(1 for o in obs_out if o[7] == 0)
        log(f'Records in parcels: {placed_p:,} · on roadsides: {placed_r:,} · elsewhere: {unplaced:,} · '
            f'left out (GPS > {MAX_ACCURACY_M} m): {coarse:,} · outside the area: {outside:,}')

    # ---- 5. native veg clipped to parcels with records ----
    with Stage('Native veg in parcels with records'):
        recorded = sorted({o[8] for o in obs_out if o[7] == 1})
        veg_tree = STRtree(veg_g)
        parcel_veg = {}
        for k, pi in enumerate(recorded):
            p = parcel_g[pi]
            c = veg_tree.query(p, predicate='intersects')
            pieces = []
            for v in c:
                g = polyonly(safe(shapely.intersection, p, veg_g[v]))
                if g is not None and g.area > 0:
                    pieces.append(g)
            parcel_veg[pi] = pieces
        log(f'{len(recorded):,} parcels with records, {sum(len(v) for v in parcel_veg.values()):,} veg pieces')

    # ---- 6. Landcare group of every parcel, roadside section and record ----
    with Stage('Assign Landcare groups'):
        g_geoms = [g['geom'] for g in groups]
        g_area = np.array([g.area for g in g_geoms])
        g_tree = STRtree(g_geoms)

        def assign(geoms):
            n = len(geoms)
            out = np.full(n, -1, int)
            if not n:
                return out
            rep = shapely.point_on_surface(np.asarray(geoms, dtype=object))
            ii, gi = g_tree.query(rep, predicate='intersects')
            best = np.full(n, np.inf)
            for i, g in zip(ii, gi):          # overlapping groups: smallest group wins
                if g_area[g] < best[i]:
                    best[i], out[i] = g_area[g], g
            miss = np.nonzero(out < 0)[0]
            if len(miss):                     # straddlers: group with the largest overlap
                ii, gi = g_tree.query(np.asarray(geoms, dtype=object)[miss], predicate='intersects')
                ov = {}
                for i, g in zip(ii, gi):
                    a = safe(shapely.intersection, geoms[miss[i]], g_geoms[g])
                    a = a.area if a is not None else 0
                    if a > ov.get(i, (-1, -1))[0]:
                        ov[i] = (a, g)
                for i, (a, g) in ov.items():
                    out[miss[i]] = g
            return out

        parcel_group = assign(parcel_g)
        roadside_group = assign(road_g_m)
        obs_group = np.full(len(obs_out), -1, int)
        if obs_out:
            ii, gi = g_tree.query(np.array(obs_pts, dtype=object), predicate='intersects')
            best = np.full(len(obs_out), np.inf)
            for i, g in zip(ii, gi):
                if g_area[g] < best[i]:
                    best[i], obs_group[i] = g_area[g], g
            for i, o in enumerate(obs_out):
                if obs_group[i] < 0 and o[7]:
                    obs_group[i] = parcel_group[o[8]] if o[7] == 1 else roadside_group[o[8]]
        # group area & native veg area
        veg_arr = np.array(veg_g, dtype=object)
        def veg_ha(poly):
            shapely.prepare(poly)
            c = veg_tree.query(poly, predicate='intersects')
            if not len(c):
                return 0.0
            within = shapely.within(veg_arr[c], poly)
            ha = float(shapely.area(veg_arr[c][within]).sum())
            for v in c[~within]:
                x = safe(shapely.intersection, veg_g[v], poly)
                ha += x.area if x is not None else 0
            return ha / 1e4

        group_stats = [[round(gg.area / 1e4, 1), round(veg_ha(gg), 1)] for gg in g_geoms]
        area_veg_ha = round(veg_ha(area), 1)
        outside_stats = None
        if catchment is not None:
            rest = polyonly(area.difference(unary_union(g_geoms)))
            outside_stats = [round(rest.area / 1e4, 1), round(veg_ha(rest), 1)] if rest is not None else [0, 0]
            log(f'  not in a Landcare group: {outside_stats[0]:,.0f} ha, native veg (2005) {outside_stats[1]:,.0f} ha')
        log(f'  whole area: {area.area / 1e4:,.0f} ha, native veg (2005) {area_veg_ha:,.0f} ha')
        for g, s in zip(groups, group_stats):
            log(f'  {g["name"]}: {s[0]:,.0f} ha, native veg (2005) {s[1]:,.0f} ha')

    # ---- 7. package ----
    with Stage('Encode and write'):
        evc_idx, evcs, veg_out = {}, [], []
        veg_enc = encode_m(veg_g, VEG_SIMPLIFY_M)
        for f, enc in zip(veg_f, veg_enc):
            p = f['p']
            evc = p.get('evc')
            evc = int(evc) if isinstance(evc, (int, float)) and float(evc).is_integer() else evc
            k = (p.get('x_evcname'), evc)
            if k not in evc_idx:
                evc_idx[k] = len(evcs)
                evcs.append([p.get('x_evcname'), evc, []])
            e = evcs[evc_idx[k]]
            if p.get('evc_bcs_desc') and p['evc_bcs_desc'] not in e[2]:
                e[2].append(p['evc_bcs_desc'])
            veg_out.append([evc_idx[k], enc])
        parcel_enc = encode_m(parcel_g, OTHER_SIMPLIFY_M)
        parcels_out = [[f['p'].get('parcel_spi') or '', f['p'].get('parcel_pfi') or '', 1 if f['p'].get('parcel_road') == 'Y' else 0,
                        round(g.area / 1e4, 2), enc] for f, g, enc in zip(parcel_f, parcel_g, parcel_enc)]
        parcel_veg_out = {}
        for pi in recorded:
            pieces = parcel_veg[pi]
            parcel_veg_out[str(pi)] = [round(sum(g.area for g in pieces) / 1e4, 2), [e for e in encode_m(pieces, OTHER_SIMPLIFY_M) if e]]
        addr_by_prop = {}
        for a in addr_f:
            k = a['p'].get('property_pfi')
            if k and a['p'].get('ezi_address'):
                lst = addr_by_prop.setdefault(k, [])
                if a['p']['ezi_address'] not in lst:
                    lst.append(a['p']['ezi_address'])
        prop_enc = encode_m(prop_g, OTHER_SIMPLIFY_M)
        props_out = [[f['p'].get('prop_pfi'), f['p'].get('prop_propnum') or '', addr_by_prop.get(f['p'].get('prop_pfi'), []), enc]
                     for f, enc in zip(prop_f, prop_enc)]
        road_enc = encode_m(road_g_m, OTHER_SIMPLIFY_M)
        roadsides_out = [[c['road'], round(c['lengthM']) if c['lengthM'] else 0, round(c['geom'].area / 1e4, 3), enc]
                         for c, enc in zip(roadsides, road_enc)]
        group_enc = encode_m(g_geoms, GROUP_SIMPLIFY_M)
        groups_out = [[g['name'], g['abbr'], enc] for g, enc in zip(groups, group_enc)]

        ll = bbox
        cx, cy = (ll[0] + ll[2]) / 2, (ll[1] + ll[3]) / 2
        w_km = (ll[2] - ll[0]) * 111.320 * math.cos(math.radians(cy))
        h_km = (ll[3] - ll[1]) * 110.574
        output = {
            'meta': {
                'built': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
                'centre': [round(cy, 6), round(cx, 6)], 'sizeKm': round(max(w_km, h_km), 1), 'bbox': bbox,
                'counts': {'parcels': len(parcels_out), 'properties': len(props_out), 'veg': len(veg_out), 'roadsides': len(roadsides_out),
                           'obs': len(obs_out), 'coarse': coarse, 'unplaced': unplaced, 'outsideArea': outside},
                'lookup': lookup_text, 'roadsideSectionM': ROADSIDE_SECTION_M,
                'mode': mode, 'areaKm2': round(area.area / 1e6, 1), 'areaHa': round(area.area / 1e4, 1), 'areaNativeVegHa': area_veg_ha, 'builder': 'build_biolinks_data.py',
            },
            'lifeforms': LIFEFORMS,
            'excludedLifeforms': EXCLUDED,
            'taxa': {str(k): v for k, v in taxa.items()},
            'evcs': evcs, 'veg': veg_out,
            'parcels': parcels_out, 'parcelVeg': parcel_veg_out, 'props': props_out,
            'roadsides': roadsides_out, 'obs': obs_out,
            # additions for the group view (existing fields unchanged)
            'groups': groups_out, 'groupStats': group_stats,
            'catchment': ['Upper Wimmera', encode_m([catchment], GROUP_SIMPLIFY_M)[0], outside_stats] if catchment is not None else None,
            'parcelGroup': parcel_group.tolist(), 'roadsideGroup': roadside_group.tolist(), 'obsGroup': obs_group.tolist(),
        }
        text = ('// Biolinks map data – built ' + output['meta']['built'] + ' by build_biolinks_data.py\n'
                'window.BIOLINKS_DATA = ' + json.dumps(output, separators=(',', ':')) + ';\n')
        out = Path(args.out)
        tmp = out.with_suffix('.tmp')
        tmp.write_text(text, encoding='utf-8')
        os.replace(tmp, out)
        log(f'Wrote {out.name}: {len(text) / 1048576:.1f} MB')

        new_taxa = [(k, t) for k, t in taxa.items() if t[5] != 'lookup' and t[6]]
        with open(args.unclassified, 'w', encoding='utf-8-sig', newline='') as f:
            w = csv.writer(f, quoting=csv.QUOTE_ALL)
            w.writerow(['taxon_id', 'scientific_name', 'common_name', 'origin', 'lifeform_code', 'lifeform', 'counts_as_native_understorey', 'note'])
            for k, t in sorted(new_taxa, key=lambda x: x[1][0]):
                w.writerow([k, t[0], t[1], 'introduced' if t[3] == 'I' else 'native to Australia but not local (likely planted)' if t[3] == 'P' else 'native',
                            t[2], LIFEFORMS.get(t[2], '') if t[2] else '', 'Y' if t[4] else 'N',
                            'guessed from genus – please check' if t[5] == 'genus' else 'not yet classified'])
        log(f'Wrote {Path(args.unclassified).name}: {len(new_taxa):,} species-level taxa not in the lookup')

    # ---- summary ----
    def stats(units_kind, gfilter=None):
        st = {}
        for gi, o in enumerate(obs_out):
            if o[7] != units_kind:
                continue
            s = st.setdefault(o[8], {'n': 0, 'all': set(), 'und': set(), 'intro': set(), 'lfs': set()})
            t = taxa[o[2]]
            s['n'] += 1
            if t[4] and t[2] and t[2] not in EXCLUDED:
                s['lfs'].add(t[2])
            if not t[6]:
                continue
            s['all'].add(o[2])
            if t[4]:
                s['und'].add(o[2])
            if t[3] == 'I':
                s['intro'].add(o[2])
        return st

    pst, rst = stats(1), stats(2)
    print()
    print('Summary')
    print(f'  parcels {len(parcels_out):,} · properties {len(props_out):,} · NV2005 patches {len(veg_out):,} · roadside sections {len(roadsides_out):,}')
    print(f'  records used {len(obs_out):,}: parcels {placed_p:,} · roadsides {placed_r:,} · unplaced {unplaced:,}')
    print(f'  dropped: GPS accuracy > {MAX_ACCURACY_M} m {coarse:,} · outside the area {outside:,}')
    print(f'  plant taxa {len(taxa):,} · native understorey {sum(1 for t in taxa.values() if t[4]):,} · species not in lookup {len(new_taxa):,}')
    print(f'  {"group":38} {"parcels w/ rec":>15} {"roadsides w/ rec":>17} {"records":>8} {"und spp":>8} {"plant grps":>10}')
    rows = list(enumerate(groups)) + ([(-1, {'name': '(not in a Landcare group)'})] if catchment is not None else [])
    for gi, g in rows:
        pt = int((parcel_group == gi).sum())
        pr = sum(1 for k in pst if parcel_group[k] == gi)
        rt = int((roadside_group == gi).sum())
        rr = sum(1 for k in rst if roadside_group[k] == gi)
        und, lfs, nrec = set(), set(), 0
        for o, og in zip(obs_out, obs_group):
            if og != gi:
                continue
            nrec += 1
            t = taxa[o[2]]
            if t[4] and t[2] and t[2] not in EXCLUDED:
                lfs.add(t[2])
            if t[4] and t[6]:
                und.add(o[2])
        print(f'  {g["name"]:38} {f"{pr:,}/{pt:,}":>15} {f"{rr:,}/{rt:,}":>17} {nrec:>8,} {len(und):>8,} {len(lfs):>10}')
    for spi in args.check_parcel:
        hits = [i for i, p in enumerate(parcels_out) if p[0] == spi]
        if not hits:
            print(f'  parcel {spi}: not in the output')
        for i in hits:
            s = pst.get(i)
            if s:
                print(f'  parcel {spi} ({parcels_out[i][3]} ha): {s["n"]} records · {len(s["lfs"])} understorey plant groups '
                      f'({", ".join(sorted(s["lfs"]))}) · {len(s["und"])} native understorey species · {len(s["all"])} species · {len(s["intro"])} introduced')
            else:
                print(f'  parcel {spi} ({parcels_out[i][3]} ha): no records')
    print(f'  output {Path(args.out).name}: {len(text) / 1048576:.1f} MB · total time {time.time() - T0:.0f} s')


if __name__ == '__main__':
    main()
