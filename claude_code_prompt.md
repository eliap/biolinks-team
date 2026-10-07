I'm a Landcare facilitator (Project Platypus, Upper Wimmera, Victoria). In this folder there's a working pilot of a habitat map. I want you to scale it up to the whole Upper Wimmera Landcare area. Do the heavy processing with a Python script on this computer, not in the browser.

## What exists now (read these first)
- `pentland_creek_property_map.html`: the map viewer (Leaflet + Turf from cdnjs). It loads `biolinks_data.js` (`window.BIOLINKS_DATA = {...}`) from the same folder through a `<script src>` tag, so it works when opened straight from disk (file://). Keep it working that way, with no local server and no fetch() of local files.
- `biolinks_data_builder.html`: the browser-based builder that currently makes `biolinks_data.js` for a 15 × 15 km test square. **Use it as the reference implementation**: its logic, thresholds and output format are what I want reproduced in Python. Read it carefully.
- `understorey_lookup.csv`: plant classification table (taxon_id, scientific_name, common_name, origin, lifeform_code, lifeform, counts_as_native_understorey, note). It is the source of truth for native/introduced/planted status and Habitat Hectares-style lifeform groups. I may edit it by hand.
- The Landcare group boundaries are in `All Upper Wimmera Landcares.shp` (+ .dbf/.shx/.prj, GDA94 / MGA zone 54, EPSG:28354) or the `.gpkg`. Ask me where they are if you can't find them. There are 13 polygons. Use the `Group` / `GROUP_NAME` and `Group_Abbr` fields. Crowlands appears twice (two polygons), so dissolve by group name. The total area is about 2,540 km².

## Goal
1. Write `build_biolinks_data.py`. It downloads everything for the area covered by the Landcare group polygons, processes it, and writes `biolinks_data.js` in the **same schema the current viewer reads** (plus the additions below).
2. Update the viewer for the bigger area: show group boundaries, add a group dropdown and per-group summary, and keep it fast.
3. Keep the browser builder working for small test areas. Don't break it, but the Python script becomes the main way to build.

## Data sources (these all work from a normal internet connection)
**Vicmap / DEECA WFS**: `https://opendata.maps.vic.gov.au/geoserver/wfs`, `service=WFS&version=2.0.0&request=GetFeature&outputFormat=application/json&srsName=EPSG:4326&bbox=W,S,E,N,EPSG:4326` (lon/lat order as written worked for us).
- Layers (prefix `open-data-platform:`): `v_property_mp` (properties), `v_parcel_mp` (parcels), `address`, `nv2005_evcbcs` (native vegetation 2005 EVCs), `road_casement_polygon` (road reserves), `tr_road` (road centrelines).
- **Gotcha:** paging with `startIndex` fails on some layers ("Cannot do natural order without a primary key") unless you also pass `sortBy`. Keys that work: v_parcel_mp→`parcel_pfi`, v_property_mp→`prop_pfi`, address→`pfi`, road_casement_polygon→`pfi`, tr_road→`ufi`; nv2005_evcbcs pages fine without one. Use `count=5000` pages.
- Expect roughly 32,000 parcels and 85,000 NV2005 patches in the bounding box (about half that inside the group polygons). Fetch by bbox in tiles if single requests are slow, de-duplicate by PFI/ID, and keep only features that intersect the dissolved group area.
- Cache raw downloads to a `cache/` subfolder so re-runs don't re-download unless I pass `--refresh`.

**iNaturalist API**: `https://api.inaturalist.org/v1/observations` with `iconic_taxa=Plantae&verifiable=true&geoprivacy=open&taxon_geoprivacy=open&preferred_place_id=7830&locale=en&per_page=200&order_by=id&order=asc&id_above=<last id>` plus swlat/swlng/nelat/nelng.
- Stay at about **1 request per second**, retry with backoff on 429/5xx, and send a descriptive User-Agent. Tile the bbox if needed and de-duplicate by observation id.
- Obscured records are deliberately excluded because their locations are randomised. Drop records with `positional_accuracy` > 200 m and report the count.

## Processing (match the browser builder)
- **Work in metres:** use EPSG:3111 (VicGrid94) or MGA54 for all geometry operations, and use shapely STRtree / geopandas spatial joins for speed.
- **Roadside sections:** cut road reserve polygons into ~500 m sections along the road centrelines.
  - Chunk each centreline into 500 m pieces, buffer each piece 40 m, intersect with the casement, and subtract any area already assigned to earlier sections so sections don't overlap.
  - Leftover casement slivers under 3,000 m² that touch a section merge into it. Larger leftovers become "Road reserve with no mapped road (unformed)". Drop pieces under 200 m².
  - Keep road name (`ezi_road_name_label`), section length and area.
- **Placing records:** assign each iNat record to the parcel containing it. If none, assign it to the roadside section containing it; otherwise it's unplaced. Subspecies roll up to species (`min_species_taxon_id` when `rank_level < 10`).
- **Classifying taxa** by `understorey_lookup.csv`:
  1. Use the exact taxon_id if it's in the lookup.
  2. Otherwise use the genus, if every lookup species in that genus shares one lifeform and origin.
  3. Otherwise mark it "unclassified": native unless iNat says introduced, with no lifeform.
  - TREE and MIS lifeforms never count as understorey.
  - Write `unclassified_plants.csv` listing species-level taxa not in the lookup (same columns as the lookup) so I can classify them and re-run.
- **Native veg clipped to parcels:** only for parcels that have records. The viewer draws all other veg as a "no records" base layer.
- **Simplify NV2005 patches** at ~1 m tolerance (they're blocky raster-derived shapes) before encoding.

## Output: `biolinks_data.js` schema (the current viewer depends on this)
`window.BIOLINKS_DATA = { meta, lifeforms, excludedLifeforms, taxa, evcs, veg, parcels, parcelVeg, props, roadsides, obs }`
- Geometry encoding: every polygon is stored as a MultiPolygon → array of polygons → array of rings → flat integer array `[x0, y0, dx1, dy1, ...]`. x,y are `round(lon*1e5)`, `round(lat*1e5)`; the first pair is absolute and the rest are deltas. Skip consecutive duplicate points, drop rings with fewer than 4 points, and drop polygons whose outer ring was dropped.
- `meta`: `{built (ISO), centre [lat,lng], sizeKm, bbox [W,S,E,N], counts {parcels, properties, veg, roadsides, obs, coarse, unplaced}, lookup (text), roadsideSectionM}`
- `lifeforms`: `{code: name}`; `excludedLifeforms`: `["TREE","MIS"]`
- `taxa`: `{id: [name, commonName, lifeformCode or "", origin "N"/"I"/"P", understorey 1/0, source "lookup"/"genus"/"unclassified", isSpeciesLevel 1/0]}`
- `evcs`: `[[evcName, evcNumber, [bioregional conservation status descriptions]]]`; `veg`: `[[evcIndex, geom]]`
- `parcels`: `[[parcel_spi, parcel_pfi, isRoad 1/0, hectares, geom]]`; `parcelVeg`: `{parcelIndex: [nativeVegHa, [geom, ...]]}`
- `props`: `[[prop_pfi, prop_propnum, [addresses], geom]]`
- `roadsides`: `[[roadName, lengthM, hectares, geom]]`
- `obs`: `[[lon, lat (5 dp), taxonId, rankLevel, observedOn, researchGrade 1/0, observationId, kind (0 none, 1 parcel, 2 roadside), unitIndex]]`
- **New:** `groups: [[groupName, abbreviation, geom]]`, and give every parcel, roadside and observation its group index. Either add an extra trailing field to each record, or include separate arrays (e.g. `parcelGroup`, `roadsideGroup`); keep the existing fields in place.

## Viewer changes
- **Group boundaries:** draw the group boundaries as a labelled outline layer (on by default) instead of the "Study area" square. Fit the map to the groups on load.
- **Group dropdown:** "All groups" or one group. Choosing a group zooms to it and shows a summary in the panel:
  - area and native veg hectares
  - parcels and roadside sections with records / total
  - native understorey species and plant groups recorded across the group
  - top 5 parcels and top 5 roadside sections by the current colour mode
- **Speed:** it must stay responsive with ~50k veg patches and ~20k roadside sections. It already uses canvas renderers; if that isn't enough, only draw veg/parcels/roadsides at zoom ≥ 12 and show the group-level view below that.
- **Unchanged behaviour:**
  - colour modes: understorey plant groups (default), native understorey species, all species, introduced species
  - "native veg in parcel" vs "whole parcel" switch
  - greyscale faded satellite base with Esri topo / CARTO street alternatives (do NOT use tile.openstreetmap.org; it blocks file:// pages)
  - yellow-green → deep turquoise-blue colour bins
  - sand "no records" fill
  - default layers: plants per parcel + plants per roadside on; records, properties, parcels and EVCs off
  - click popups for parcel, roadside and property

## How to work
- Before writing code, tell me your plan and an estimate of run time and output file size. Ask me before installing anything big. Use a virtual environment (there's a `.venv` here already; check what's in it).
- Build and verify on **one small group first** (Stawell Urban or Navarre) with a `--group` option. Then run the full catchment.
- Verify by:
  - opening the map in a browser and checking the console for errors
  - spot-checking that parcel S6\PP3207 (near Rhymney, 135 ha) shows 8 understorey plant groups and ~39 native understorey species, matching the browser builder's result for the same area
  - reporting counts (parcels, veg patches, roadside sections, records placed / unplaced / dropped for accuracy) and the final file size
- Print progress and timings for each stage. The full run should be restartable from the cache.
- Don't delete or overwrite my other files. Keep the old viewer as `pentland_creek_property_map_pilot.html` before changing it.
