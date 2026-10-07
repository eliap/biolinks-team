# Biolinks habitat map – Upper Wimmera

Native understorey and habitat by parcel and roadside section for the Upper Wimmera (Project Platypus Biolinks Team), from iNaturalist plant records, Vicmap property and road data and DEECA native vegetation 2005 (EVC) mapping.

**Open the map:** https://eliap.github.io/biolinks-team/inaturalist-map/

## Files
- `inaturalist-map/index.html` – the map viewer (Leaflet): where has nobody recorded plants on iNaturalist yet? Loads `biolinks_data.js` from the same folder; also works opened straight from disk.
- `inaturalist-map/biolinks_data.js` – pre-processed map data.
- `pentland_creek_property_map.html` – old link, redirects to the map.
- `build_biolinks_data.py` – builds `biolinks_data.js` for the whole Upper Wimmera (or one group / a test square).
- `biolinks_data_builder.html` – browser builder for small test squares (reference implementation).
- `understorey_lookup.csv` – native / introduced / planted status and Habitat Hectares lifeform groups. Edit by hand, then rebuild.
- `unclassified_plants.csv` – species recorded but not yet in the lookup (same columns), to classify and copy into the lookup.
- `landcare/` – Landcare group boundaries (MGA zone 54) and the 2007 Wimmera Landcare network boundaries (Upper Catchment outline).
- `pentland_creek_property_map_pilot.html` – the original 15 km pilot viewer.

## Rebuilding the data
```
python -m venv .venv
.venv\Scripts\python -m pip install shapely pyproj pyshp requests
.venv\Scripts\python build_biolinks_data.py                 # whole Upper Wimmera (~35 min first time, ~3 min from cache)
.venv\Scripts\python build_biolinks_data.py --refresh-inat  # pick up new iNaturalist records
.venv\Scripts\python build_biolinks_data.py --group "Navarre"
.venv\Scripts\python build_biolinks_data.py --groups-only   # only land inside the Landcare groups
```
Downloads are cached in `cache/` (not in the repo). The script writes `inaturalist-map/biolinks_data.js`; commit and push it to update the online map.

## Notes
- iNaturalist records with obscured locations are not included, and records with GPS accuracy worse than 200 m are left out.
- "No records" means not surveyed, not that there are no plants.
- Native vegetation 2005 is mapped at 1:100,000 – check on the ground.
