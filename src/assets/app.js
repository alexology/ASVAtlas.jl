"use strict";

const $ = (id) => document.getElementById(id);
let dataset = null;
let current = null;   // last search response
let selected = null;  // query shown alone on the map, or null for all
let map = null;
let layer = null;

function el(tag, attrs, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") node.className = v;
    else node.setAttribute(k, v);
  }
  for (const c of children) {
    if (c === null || c === undefined) continue;
    node.append(c instanceof Node ? c : String(c));
  }
  return node;
}

const pct = (x) => (Math.floor(x * 100) / 100).toFixed(2) + "%";
const bucket = (id) => (id >= 100 - 1e-9 ? "b100" : id >= 99 ? "b99" : id >= 97 ? "b97" : "blow");
const colorOf = (b) => getComputedStyle(document.documentElement).getPropertyValue("--" + b).trim();

const fmt = (n) => Number(n).toLocaleString("en");

function setStatus(text, error) {
  $("status").textContent = text;
  $("status").className = error ? "error" : "";
}

// Taxon names are italic only at genus level and below.
function taxonNode(name, rankIndex) {
  if (!name) return el("span", { class: "lineage" }, "unknown");
  const ranks = (current && current.ranks) || (dataset && dataset.ranks) || [];
  const r = ranks[rankIndex] || "";
  const italic = ["genus", "subgenus", "species", "subspecies"].includes(r);
  return el("span", { class: italic ? "taxon" : "" }, name);
}

function deepest(lineage) {
  for (let i = lineage.length - 1; i >= 0; i--) if (lineage[i]) return [lineage[i], i];
  return ["", -1];
}

async function loadDataset() {
  try {
    const res = await fetch("api/dataset");
    dataset = await res.json();
  } catch (e) {
    setStatus("The server cannot be reached.", true);
    return;
  }
  $("stat-name").textContent = dataset.name;
  $("stat-name").title = dataset.ranks.length ? `Taxonomy: ${dataset.ranks.join(", ")}` : "No taxonomy";
  // Marker and primers of the reference, when known: another marker or primer
  // pair is another dataset.
  for (const k of ["marker", "primers"]) {
    $(`stat-${k}`).textContent = dataset[k] || "—";
    $(`stat-${k}-box`).hidden = !dataset[k];
  }
  $("stat-references").textContent = fmt(dataset.references);
  $("stat-georeferenced").textContent = fmt(dataset.georeferenced);
  $("stat-occurrences").textContent = fmt(dataset.occurrences);
  $("start-name").textContent = dataset.name +
    (dataset.marker || dataset.primers ? ` (${[dataset.marker, dataset.primers && `primers ${dataset.primers}`]
      .filter((x) => x).join(", ")})` : "");
  $("example").hidden = !dataset.example;
  // Similarity levels, from identical records (100%) down to the server's
  // minimum; changing the level repeats the search at once.
  const group = $("identity-group");
  for (const v of dataset.identity_levels) {
    const input = el("input", { type: "radio", name: "identity", value: v });
    input.checked = v === 100;
    // A new level filters the results already in the page (no new search).
    input.addEventListener("change", () => { if (full) showLevel(v); });
    group.append(el("label", { title: v === 100 ? "Identical records" : `At least ${v}% similar` },
      input, el("span", {}, v + "%")));
  }
  $("max-hits").max = dataset.max_hits;
}

function chosenIdentity() {
  const checked = document.querySelector('input[name="identity"]:checked');
  return checked ? Number(checked.value) : 100;
}

function setIdentity(v) {
  const input = document.querySelector(`input[name="identity"][value="${v}"]`);
  if (input) input.checked = true;
}

// The next, wider similarity level, or null at the server's minimum.
function widerLevel(v) {
  const levels = (dataset && dataset.identity_levels) || [];
  return levels.find((x) => x < v) ?? null;
}

const levelText = (v) => (v === 100 ? "identical" : `≥ ${v}% similar`);

// Number of sequences in the text area: FASTA headers, or non-empty lines.
function countSequences() {
  const lines = $("sequences").value.split(/\r?\n/).map((l) => l.trim()).filter((l) => l);
  const n = lines.length && lines[0].startsWith(">") ? lines.filter((l) => l.startsWith(">")).length
                                                      : lines.length;
  $("seq-count").textContent = n ? `${fmt(n)} sequence${n === 1 ? "" : "s"}` : "No sequences";
  return n;
}

function readFile(f) {
  f.text().then((t) => {
    $("sequences").value = t;
    countSequences();
    updateClear();
    setStatus(`Loaded ${f.name}.`);
  });
}

// Two base maps: the bundled world map (Natural Earth countries, public
// domain; works offline) and OpenStreetMap tiles (online, more detailed).
let worldData = null;   // the world map, fetched once for all maps

function initMap() {
  if (map) return;
  if (typeof L === "undefined") {
    showNote("The map library could not be loaded.");
    return;
  }
  ({ map, layer } = makeMap("map", showNote));
}

// Country of a place, from the polygons of the world map (null at sea, on a
// coast left out by the 1:50m outlines, or before the map is loaded).
let countryIndex = null;
const countryCache = new Map();

function buildCountryIndex(data) {
  countryIndex = data.features.map((f) => {
    const g = f.geometry;
    const polys = g.type === "Polygon" ? [g.coordinates] : g.type === "MultiPolygon" ? g.coordinates : [];
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const poly of polys) for (const [x, y] of poly[0]) {
      x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y);
    }
    return { name: f.properties.name, polys, box: [x0, y0, x1, y1] };
  });
}

function inRing(x, y, ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

function countryAt(lat, lon) {
  if (!countryIndex) return null;
  const key = lat.toFixed(4) + "," + lon.toFixed(4);
  if (countryCache.has(key)) return countryCache.get(key);
  let name = null;
  for (const c of countryIndex) {
    const [x0, y0, x1, y1] = c.box;
    if (lon < x0 || lon > x1 || lat < y0 || lat > y1) continue;
    if (c.polys.some((p) => inRing(lon, lat, p[0]) && !p.slice(1).some((h) => inRing(lon, lat, h)))) {
      name = c.name;
      break;
    }
  }
  countryCache.set(key, name);
  return name;
}

// "Norway", "Norway and Italy", "Italy, Norway and 3 more countries".
function countryList(names) {
  const list = [...names].filter((x) => x).sort();
  if (!list.length) return "";
  if (list.length === 1) return list[0];
  if (list.length === 2) return `${list[0]} and ${list[1]}`;
  if (list.length === 3) return `${list[0]}, ${list[1]} and ${list[2]}`;
  return `${list[0]}, ${list[1]} and ${list.length - 2} more countries`;
}

// A map with the two base maps and a layer for the markers.
function makeMap(id, note, layersPosition = "topright") {
  const map = L.map(id, {
    minZoom: 2, maxZoom: 18,
    maxBounds: [[-90, -240], [90, 240]], maxBoundsViscosity: 1,
  }).setView([30, 10], 2);
  map.getContainer().classList.add("world");
  const world = L.geoJSON(null, {
    renderer: L.canvas({ padding: 0.5 }),
    interactive: false,
    style: () => ({
      color: colorOf("map-border"), weight: 0.6,
      fillColor: colorOf("map-land"), fillOpacity: 1,
    }),
    attribution: 'World map: <a href="https://www.naturalearthdata.com/">Natural Earth</a>',
  }).addTo(map);
  worldData = worldData || loadWorld();
  worldData.then((data) => world.addData(data))
    .catch(() => note("The world map could not be loaded."));
  const osm = L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
    maxZoom: 18,
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
  });
  L.control.layers({ "World map (offline)": world, "OpenStreetMap (online)": osm },
    null, { position: layersPosition }).addTo(map);
  map.on("baselayerchange", (e) =>
    map.getContainer().classList.toggle("world", e.layer === world));
  return { map, layer: L.layerGroup().addTo(map) };
}

function showNote(text) {
  $("map-note").textContent = text;
  $("map-note").hidden = !text;
}

// One marker per location; its popup lists the reference ASVs found there
// and the queries matching them. Returns the number of places.
// ---- Sampling locations of the user's ASVs ---------------------------------------
// Optional coordinates from a CSV/TSV table (id, latitude, longitude). They
// stay in the page: only the sequences are sent to the server.

const USER_COLOR = "#f5c542";   // gold: apart from the aqua ramp and the network colours
let coordTable = new Map();     // id -> [{lat, lon}, ...]: an ASV can be sampled at several places

const validLat = (x) => Number.isFinite(x) && x >= -90 && x <= 90;
const validLon = (x) => Number.isFinite(x) && x >= -180 && x <= 180;

// Records of the text area as the server reads them: identifier = first word
// of the FASTA header (query_1, query_2, ... for plain sequences); sequence
// in upper case, U as T, without spaces and gaps.
function queryRecords(text) {
  const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l);
  if (!lines.length) return [];
  const clean = (s) => s.toUpperCase().replace(/[\s.-]/g, "").replace(/U/g, "T");
  if (!lines[0].startsWith(">")) return lines.map((l, k) => ({ id: `query_${k + 1}`, seq: clean(l) }));
  const out = [];
  for (const l of lines) {
    if (l.startsWith(">")) out.push({ id: l.slice(1).trim().split(/\s+/)[0], seq: "" });
    else if (out.length) out[out.length - 1].seq += clean(l);
  }
  return out;
}

const queryIds = (text) => queryRecords(text).map((r) => r.id);

// Each ASV must appear once: the same sequence under two identifiers is an
// error (an ASV found at several places has several rows in the table of
// coordinates instead). Returns the sampling places of each identifier.
function checkQueries(text) {
  const records = queryRecords(text);
  const bySeq = new Map();
  for (const r of records) {
    if (!r.seq) continue;
    if (!bySeq.has(r.seq)) bySeq.set(r.seq, []);
    bySeq.get(r.seq).push(r.id);
  }
  const repeated = [...bySeq.values()].filter((ids) => ids.length > 1);
  if (repeated.length) {
    const list = repeated.slice(0, 3).map((ids) => ids.join(" = ")).join("; ");
    throw new Error(`Each sequence must appear once; identical sequences: ${list}` +
      (repeated.length > 3 ? ` and ${repeated.length - 3} more` : "") +
      ". For an ASV found at several places, give one identifier and several rows in the coordinate table.");
  }
  const coords = new Map();
  for (const r of records) if (coordTable.has(r.id)) coords.set(r.id, coordTable.get(r.id));
  return coords;
}

// Table of coordinates: a header with exactly the columns asv_id, latitude
// and longitude (in any order), then one row per sampling place (an ASV found
// at several places has several rows). Comma, semicolon or tab separated.
// Coordinates are decimal degrees in WGS84 (e.g. 47.393030, as GPS and web
// maps give them), with a point or, in semicolon or tab separated tables, a
// comma as decimal separator.
const DECIMAL = /^[+-]?\d+(\.\d+)?$/;
const COORD_COLUMNS = ["asv_id", "latitude", "longitude"];

function readCoordTable(text) {
  const rows = text.replace(/^\uFEFF/, "").split(/\r?\n/).filter((l) => l.trim());
  if (!rows.length) throw new Error("The coordinate table is empty.");
  const sep = rows[0].includes("\t") ? "\t" : rows[0].includes(";") ? ";" : ",";
  const split = (l) => l.split(sep).map((x) => x.trim().replace(/^"(.*)"$/, "$1").trim());
  const decimal = (x) => {
    const t = String(x ?? "").replace(",", ".");
    return DECIMAL.test(t) ? Number(t) : NaN;
  };
  const head = split(rows[0]).map((h) => h.toLowerCase());
  if (head.length !== 3 || !COORD_COLUMNS.every((c) => head.includes(c))) {
    throw new Error("The coordinate table must have a header with exactly the three columns " +
      `asv_id, latitude and longitude; found: ${split(rows[0]).join(", ")}.`);
  }
  const cols = COORD_COLUMNS.map((c) => head.indexOf(c));
  const start = 1;
  const table = new Map();
  for (const [k, l] of rows.slice(start).entries()) {
    const f = split(l);
    const where = `Coordinate table, row ${k + start + 1}`;
    const id = f[cols[0]];
    if (!id) throw new Error(`${where}: no identifier.`);
    const lat = decimal(f[cols[1]]), lon = decimal(f[cols[2]]);
    if (!Number.isFinite(lat) || !Number.isFinite(lon)) {
      throw new Error(`${where} (${id}): coordinates must be decimal degrees (WGS84), e.g. 46.0712 and 11.1217, ` +
        `not degrees-minutes-seconds; got "${f[cols[1]] ?? ""}", "${f[cols[2]] ?? ""}".`);
    }
    if (!validLat(lat) || !validLon(lon)) {
      const swapped = validLat(lon) && validLon(lat) ? " Are latitude and longitude swapped?" : "";
      throw new Error(`${where} (${id}): latitude must be between −90 and 90 and longitude between −180 and 180.${swapped}`);
    }
    // One row per place: an ASV can have several (repeated rows count once).
    if (!table.has(id)) table.set(id, []);
    const list = table.get(id);
    if (!list.some((x) => x.lat === lat && x.lon === lon)) list.push({ lat, lon });
  }
  if (!table.size) throw new Error("The coordinate table has no rows.");
  return table;
}

// Great-circle distance in km.
function km(a, b) {
  const rad = Math.PI / 180, R = 6371.0088;
  const dlat = (b.lat - a.lat) * rad, dlon = (b.lon - a.lon) * rad;
  const h = Math.sin(dlat / 2) ** 2 + Math.cos(a.lat * rad) * Math.cos(b.lat * rad) * Math.sin(dlon / 2) ** 2;
  return 2 * R * Math.asin(Math.min(1, Math.sqrt(h)));
}

const kmText = (d) => (d < 1 ? "< 1 km" : d < 100 ? `${d.toFixed(0)} km` : `${fmt(Math.round(d / 10) * 10)} km`);

// Sampling locations of a query (possibly none).
const samplesOf = (id) => (current && current.coords && current.coords.get(id)) || [];

// Nearest matching record of a query from one sampling location:
// {d (km), o (occurrence)}, or null without georeferenced matches.
function nearestFrom(q, s) {
  let best = null;
  for (const h of q.hits) {
    for (const o of current.references[h.reference].occurrences) {
      const d = km(s, o);
      if (best === null || d < best.d) best = { d, o };
    }
  }
  return best;
}

// Distance from the nearest sampling location of a query to its nearest
// matching record, or null.
function nearestKm(q) {
  let best = null;
  for (const s of samplesOf(q.id)) {
    const n = nearestFrom(q, s);
    if (n && (best === null || n.d < best)) best = n.d;
  }
  return best;
}

// A star marker, gold with a dark outline (shape and colour both set it apart
// from the records).
// Without `fills`, the star is gold with a dark outline. With `fills` (one
// colour, or two for a star split in halves, e.g. the groups of the network),
// it is filled with them and keeps a gold border, the mark of the user's ASVs.
let starIds = 0;

function starSvg(size, fills) {
  const c = size / 2, R = size / 2 - 2, r = R * 0.45;
  const pts = [];
  for (let k = 0; k < 10; k++) {
    const a = -Math.PI / 2 + (k * Math.PI) / 5, rr = k % 2 ? r : R;
    pts.push(`${(c + rr * Math.cos(a)).toFixed(2)},${(c + rr * Math.sin(a)).toFixed(2)}`);
  }
  const poly = (attrs) => `<polygon points="${pts.join(" ")}" stroke-linejoin="round" ${attrs}/>`;
  let body;
  if (!fills) {
    body = poly(`fill="${USER_COLOR}" stroke="#03101b" stroke-width="1.5"`);
  } else {
    let fill = fills[0], defs = "";
    if (fills.length > 1) {
      const id = `star-split-${++starIds}`;
      defs = `<defs><linearGradient id="${id}"><stop offset="50%" stop-color="${fills[0]}"/>` +
        `<stop offset="50%" stop-color="${fills[1]}"/></linearGradient></defs>`;
      fill = `url(#${id})`;
    }
    // Border widths in proportion to the size (2 px gold at 24 px).
    const gold = (size / 12).toFixed(2), dark = (size / 6).toFixed(2);
    body = defs + poly(`fill="none" stroke="#03101b" stroke-width="${dark}"`) +
      poly(`fill="${fill}" stroke="${USER_COLOR}" stroke-width="${gold}"`);
  }
  return `<svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">${body}</svg>`;
}

function starIcon(size = 22, fills = null) {
  return L.divIcon({ className: "star-icon", html: starSvg(size, fills),
    iconSize: [size, size], iconAnchor: [size / 2, size / 2] });
}

// Stars for the sampling locations of some queries (one per place, shared by
// the ASVs sampled there), added to `target`; `describe(ids, s)` gives the
// lines of the tooltip of place `s` and `fillsOf(ids)`, if given, the fills
// of its star (see starIcon). Returns their coordinates.
function drawSamples(target, ids, describe, fillsOf = null) {
  const at = new Map();
  for (const id of ids) {
    for (const s of samplesOf(id)) {
      const key = s.lat.toFixed(6) + "," + s.lon.toFixed(6);
      if (!at.has(key)) at.set(key, { s, ids: [] });
      at.get(key).ids.push(id);
    }
  }
  const points = [];
  for (const { s, ids: here } of at.values()) {
    const size = fillsOf ? 24 : 22;
    const m = L.marker([s.lat, s.lon], { icon: starIcon(size, fillsOf ? fillsOf(here) : null),
      keyboard: true, zIndexOffset: 2000 });
    m.bindTooltip(el("div", { class: "place-tip" },
      el("b", {}, `Your sample${countryAt(s.lat, s.lon) ? " · " + countryAt(s.lat, s.lon) : ""}`),
      ...describe(here, s).map((line) => el("div", {}, line))),
      { direction: "top", offset: [0, -11] });
    m.addTo(target);
    points.push([s.lat, s.lon]);
  }
  return points;
}

function drawMap() {
  if (!map) return 0;
  layer.clearLayers();
  const places = new Map();
  for (const q of current.queries) {
    if (selected !== null && q.id !== selected) continue;
    for (const h of q.hits) {
      const ref = current.references[h.reference];
      ref.occurrences.forEach((o) => {
        const key = o.lat.toFixed(6) + "," + o.lon.toFixed(6);
        if (!places.has(key)) {
          places.set(key, { lat: o.lat, lon: o.lon, refs: new Map(), best: 0, country: countryAt(o.lat, o.lon) });
        }
        const p = places.get(key);
        p.best = Math.max(p.best, h.identity);
        if (!p.refs.has(h.reference)) p.refs.set(h.reference, { fields: [], queries: [] });
        const r = p.refs.get(h.reference);
        if (!r.fields.includes(o)) r.fields.push(o);
        if (!r.queries.some((x) => x.id === q.id)) r.queries.push({ id: q.id, identity: h.identity });
      });
    }
  }
  // Sampling locations of the ASVs shown (all, or the selected one).
  const shown = current.queries.filter((q) => selected === null || q.id === selected);
  const samples = drawSamples(layer, shown.map((q) => q.id), (ids, s) => ids.map((id) => {
    const n = nearestFrom(current.queries.find((x) => x.id === id), s);
    return n === null ? `${id}: no matching record with coordinates` : `${id}: nearest record ${kmText(n.d)}`;
  }));
  $("legend-user").hidden = samples.length === 0;
  $("legend").hidden = places.size === 0 && samples.length === 0;
  if (places.size === 0) {
    placeBounds = samples.length ? samples : null;
    fitPlaces();
    const q = current.queries.find((x) => x.id === selected);
    const v = current.params.identity;
    const level = v === 100 ? "identical to" : `at least ${v}% similar to`;
    const wider = widerLevel(v) === null ? "" : " Try a lower similarity above the map.";
    showNote(selected === null
      ? (current.queries.some((x) => x.hits.length) ? "The matching records have no coordinates."
        : `No record ${level} your ASVs.${wider}`)
      : q.hits.length ? `The records matching ${selected} have no coordinates.`
      : `No record ${level} ${selected}.${wider}`);
    return 0;
  }
  showNote("");
  const bounds = [];
  for (const p of places.values()) {
    const b = bucket(p.best);
    const m = L.circleMarker([p.lat, p.lon], {
      radius: 6, weight: 2, color: colorOf("map-sea"), fillColor: colorOf(b), fillOpacity: 1,
    });
    const n = [...p.refs.values()].reduce((a, r) => a + r.fields.length, 0);
    const names = [...new Set([...p.refs.values()].flatMap((r) => r.queries.map((x) => x.id)))];
    m.bindTooltip(el("div", { class: "place-tip" },
      el("b", {}, p.country || `${p.lat.toFixed(2)}, ${p.lon.toFixed(2)}`),
      el("div", {}, `${plural(n, "record")} · ${p.best >= 100 - 1e-9 ? "identical" : "best " + pct(p.best)}`),
      el("div", { class: "muted" }, names.length > 3 ? `${names.slice(0, 3).join(", ")} +${names.length - 3}`
        : names.join(", "))),
      { direction: "top", offset: [0, -6] });
    m.bindPopup(() => popup(p), { maxWidth: 360 });
    m.addTo(layer);
    bounds.push([p.lat, p.lon]);
  }
  // From each sampling location of the selected ASV to its nearest record.
  if (selected !== null) {
    const q = current.queries.find((x) => x.id === selected);
    for (const s of samplesOf(selected)) {
      const n = nearestFrom(q, s);
      if (!n) continue;
      L.polyline([[s.lat, s.lon], [n.o.lat, n.o.lon]], {
        color: USER_COLOR, weight: 1.2, opacity: 0.6, dashArray: "4 5", interactive: false,
      }).addTo(layer);
    }
  }
  placeBounds = bounds.concat(samples);
  fitPlaces();
  return places.size;
}

// Never closer than a region, so that the place can be recognised.
let placeBounds = null;

function fitPlaces() {
  if (!map || !placeBounds) return;
  if (placeBounds.length === 1) map.setView(placeBounds[0], 5);
  else map.fitBounds(placeBounds, { padding: [40, 40], maxZoom: 6 });
}

// The places stay in view when the window changes size.
let resizeTimer = null;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (!map || $("results-view").hidden) return;
    map.invalidateSize();
    fitPlaces();
  }, 200);
});

// For records whose study has no site coordinates (see the description of the
// reference, tab "Reference data"): what the point on the map stands for.
function locationNote(fields) {
  const sources = dataset && dataset.about && Array.isArray(dataset.about.sources) ? dataset.about.sources : [];
  const d = (fields.find((f) => f.name.toLowerCase() === "dataset") || {}).value;
  const s = d && sources.find((x) => x.dataset === d);
  if (!s) return "";
  return s.coordinates === "approximate" ? "Approximate location: random point in the study area"
    : s.coordinates === "study area" ? "Approximate location: one point for the whole study area" : "";
}

function popup(p) {
  const box = el("div", {},
    p.country ? el("div", { class: "popup-country" }, p.country) : null,
    el("div", { class: "popup-coords" }, `${p.lat.toFixed(5)}, ${p.lon.toFixed(5)}`));
  for (const [id, r] of p.refs) {
    const ref = current.references[id];
    const [name, i] = deepest(ref.lineage);
    const ul = el("ul");
    for (const q of r.queries) ul.append(el("li", {}, el("span", {}, q.id), el("span", {}, pct(q.identity))));
    const part = el("div", { class: "popup-ref" },
      el("h3", {}, taxonNode(name, i), " ", el("span", { class: "ref-id" }, id)), ul);
    for (const o of r.fields) {
      if (o.fields.length) {
        part.append(el("div", { class: "fields" },
          o.fields.map((f) => `${f.name}: ${f.value}`).join(" · ")));
      }
      const note = locationNote(o.fields);
      if (note) part.append(el("div", { class: "fields approx" }, note));
    }
    box.append(part);
  }
  return box;
}

function identityCell(identity) {
  return el("span", { class: "identity-cell" }, el("i", { class: "dot " + bucket(identity) }), pct(identity));
}

// Matching reference ASVs, records, places and countries of some queries.
function stats(queries) {
  const refs = new Set(), places = new Set(), countries = new Set();
  let records = 0;
  for (const q of queries) for (const h of q.hits) refs.add(h.reference);
  for (const id of refs) {
    for (const o of current.references[id].occurrences) {
      records++;
      places.add(o.lat.toFixed(6) + "," + o.lon.toFixed(6));
      const c = countryAt(o.lat, o.lon);
      if (c) countries.add(c);
    }
  }
  return { refs: refs.size, records, places: places.size, countries };
}

const plural = (n, word) => `${fmt(n)} ${word}${n === 1 ? "" : "s"}`;

// The answer above the map, for all ASVs or the selected one.
function updateSummary() {
  const v = current.params.identity;
  const qs = selected === null ? current.queries : current.queries.filter((q) => q.id === selected);
  const st = stats(qs);
  const found = qs.filter((q) => q.hits.length).length;
  const box = $("map-summary");
  box.className = "";
  box.replaceChildren();
  if (!st.records) {
    const who = selected !== null ? selected : qs.length === 1 ? qs[0].id : "Your ASVs";
    box.append(`${who}: no ${levelText(v)} record.`);
    const w = widerLevel(v);
    if (w !== null) box.append(" ", tryButton(w));
    return;
  }
  const who = selected !== null ? selected : qs.length === 1 ? qs[0].id
    : `${fmt(found)} of your ${fmt(qs.length)} ASVs`;
  const where = st.countries.size ? ` in ${countryList(st.countries)}` : "";
  box.append(el("b", {}, who), ` ${selected === null && qs.length > 1 ? "are" : "is"} found at ` +
    `${plural(st.places, "place")}${where}: ${plural(st.records, "record")} ${levelText(v)}.`);
}

function tryButton(level) {
  const b = el("button", { type: "button", class: "link-button" }, `Try ${level}%`);
  b.addEventListener("click", (ev) => {
    ev.stopPropagation();
    setIdentity(level);
    showLevel(level);
  });
  return b;
}

// Taxa of the best matches of a query.
function bestTaxa(q) {
  const best = q.hits.length ? q.hits[0].identity : null;
  const out = new Map();
  for (const h of q.hits) {
    if (h.identity !== best) continue;
    const [name, i] = deepest(current.references[h.reference].lineage);
    if (name && !out.has(name)) out.set(name, i);
  }
  return out;
}

function networkButton(q) {
  const k = current.ranks.indexOf("species");
  const name = k < 0 || !q.hits.length ? "" : current.references[q.hits[0].reference].lineage[k];
  if (!name || !(dataset && dataset.networks)) return null;
  const b = el("button", { type: "button", class: "ghost net-button", title: `Haplotype network of ${name}` },
    "Network →");
  b.addEventListener("click", (ev) => {
    ev.stopPropagation();
    openNetwork(name, q.id);
  });
  return b;
}

// One line per ASV: what it is, where it is found, how similar.
function renderResults() {
  const list = $("asv-list");
  list.replaceChildren();
  const v = current.params.identity;
  for (const q of current.queries) {
    const st = stats([q]);
    const li = el("li", { class: q.hits.length ? "asv" : "asv nohit", "data-id": q.id });
    const main = el("button", { type: "button", class: "asv-main", "aria-expanded": "false" });
    const taxa = el("span", { class: "asv-taxa" });
    if (q.hits.length) {
      const t = bestTaxa(q);
      if (!t.size) taxa.append(el("span", { class: "muted" }, "no taxonomy"));
      [...t].forEach(([name, i], k) => { if (k) taxa.append(" | "); taxa.append(taxonNode(name, i)); });
    } else {
      taxa.append(`No ${levelText(v)} record`);
    }
    const where = el("span", { class: "asv-where" });
    if (q.hits.length) {
      where.append(st.places
        ? `${plural(st.places, "place")}` + (st.countries.size ? ` in ${countryList(st.countries)}` : "") +
          ` · ${plural(st.records, "record")}`
        : `${plural(st.records, "record")}, without coordinates`);
    }
    const d = nearestKm(q);
    const ns = samplesOf(q.id).length;
    if (d !== null) {
      where.append(el("span", { class: "asv-distance" }, el("i", { class: "star" }),
        `nearest ${kmText(d)} from ${ns === 1 ? "your sample" : `your ${ns} sampling places`}`));
    } else if (ns) {
      if (q.hits.length) where.append(el("span", { class: "asv-distance none" }, "no matching record with coordinates"));
    } else if (current.tableLoaded) {
      // A table of coordinates was given, but not for this ASV: say so.
      where.append(el("span", { class: "asv-distance none" }, "no sampling place in your table"));
    }
    main.append(el("span", { class: "asv-name" }, q.id), taxa, where,
      el("span", { class: "asv-identity" }, q.hits.length ? identityCell(q.hits[0].identity) : ""));
    main.addEventListener("click", () => toggle(q, li));
    const actions = el("span", { class: "asv-actions" });
    const nb = networkButton(q);
    if (nb) actions.append(nb);
    if (!q.hits.length && widerLevel(v) !== null) actions.append(tryButton(widerLevel(v)));
    li.append(el("div", { class: "asv-row" }, main, actions));
    list.append(li);
  }
  $("edit-label").textContent = `Edit ${plural(current.queries.length, "ASV")}`;
}

function toggle(q, li) {
  const open = li.classList.contains("open");
  document.querySelectorAll("#asv-list li.open").forEach((x) => {
    x.classList.remove("open");
    x.querySelector(".asv-main").setAttribute("aria-expanded", "false");
    x.querySelector(".asv-details")?.remove();
  });
  if (open && selected === q.id) {
    selected = null;
  } else {
    selected = q.id;
    li.classList.add("open");
    li.querySelector(".asv-main").setAttribute("aria-expanded", "true");
    li.append(details(q));
  }
  $("show-all").hidden = selected === null;
  drawMap();
  updateSummary();
}

// Matches of a query grouped by taxon; the alignment details are folded.
function details(q) {
  const box = el("div", { class: "asv-details" });
  if (!q.hits.length) {
    box.append(el("p", { class: "details-note" }, `No reference ASV ${levelText(current.params.identity)}` +
      (q.exhaustive ? "." : " among the candidates examined.")));
    return box;
  }
  const groups = new Map();
  for (const h of q.hits) {
    const ref = current.references[h.reference];
    const [name, i] = deepest(ref.lineage);
    if (!groups.has(name)) groups.set(name, { name, i, lineage: ref.lineage, hits: [], records: 0 });
    const g = groups.get(name);
    g.hits.push(h);
    g.records += ref.occurrences.length;
  }
  box.append(el("h3", {}, "Matches by taxon"));
  for (const g of groups.values()) {
    const ids = g.hits.map((h) => h.identity);
    const hi = Math.max(...ids), lo = Math.min(...ids);
    const range = hi === lo ? pct(hi) : `${pct(lo)}–${pct(hi)}`;
    box.append(el("div", { class: "match-group" },
      el("div", { class: "match-head" },
        g.name ? taxonNode(g.name, g.i) : el("span", { class: "muted" }, "no taxonomy"),
        el("span", { class: "match-count" },
          `${plural(g.hits.length, "reference ASV")} · ${range} · ${plural(g.records, "record")}`)),
      el("div", { class: "lineage" }, g.lineage.filter((x) => x).join(" › "))));
  }
  const body = el("tbody");
  for (const h of q.hits) {
    const ref = current.references[h.reference];
    body.append(el("tr", {},
      el("td", {}, h.reference),
      el("td", {}, taxonNode(...deepest(ref.lineage))),
      el("td", { class: "num" }, identityCell(h.identity)),
      el("td", { class: "num" }, pct(h.coverage * 100)),
      el("td", {}, h.strand === "-" ? "reverse" : "forward"),
      el("td", { class: "num" }, `${h.reference_start}–${h.reference_end}`),
      el("td", { class: "num" }, ref.occurrences.length)));
  }
  const table = el("table", {},
    el("thead", {}, el("tr", {}, ...[["Reference ASV"], ["Taxon"], ["Similarity", "num"], ["Coverage", "num"],
      ["Strand"], ["Position on reference", "num"], ["Records", "num"]]
      .map(([h, c]) => el("th", c ? { class: c } : {}, h)))),
    body);
  box.append(el("details", { class: "tech" },
    el("summary", {}, "Technical details"),
    el("p", { class: "tech-note" }, "Coverage: share of your ASV in the alignment. Strand: reverse when " +
      "your ASV was read on the other strand. Position: aligned bases of the reference ASV."),
    el("div", { class: "table-wrap" }, table)));
  return box;
}

// Text cells starting with = + - @ would run as formulas in spreadsheet
// programs (CSV/TSV injection): they get a leading apostrophe. Numbers are
// written as numbers.
function tsv(rows) {
  const clean = (x) => {
    if (x === null || x === undefined) return "";
    if (typeof x === "number") return String(Math.round(x * 1e6) / 1e6);
    const s = String(x).replace(/[\t\r\n]/g, " ");
    return /^[=+\-@]/.test(s) ? "'" + s : s;
  };
  return rows.map((r) => r.map(clean).join("\t")).join("\n") + "\n";
}

function download(kind) {
  const ranks = current.ranks;
  let rows;
  if (kind === "summary") {
    rows = [["query", "length", "n_hits", "best_identity", "best_references", "best_taxa", "n_georeferenced",
      "n_sampling_places", "sampling_places", "nearest_record_km"]];
    for (const q of current.queries) {
      const best = q.hits.length ? q.hits[0].identity : null;
      rows.push([q.id, q.length, q.hits.length, best,
        q.hits.filter((h) => h.identity === best).map((h) => h.reference).join(";"),
        q.best_taxa.join(" | "),
        q.hits.filter((h) => current.references[h.reference].occurrences.length).length,
        samplesOf(q.id).length, samplesOf(q.id).map((s) => `${s.lat} ${s.lon}`).join("; "),
        nearestKm(q) === null ? "" : Math.round(nearestKm(q) * 10) / 10]);
    }
  } else if (kind === "hits") {
    rows = [["query", "hit", "reference", "identity", "coverage", "strand", "query_start", "query_end",
      "reference_start", "reference_end", "cigar", "taxon", ...ranks, "n_occurrences"]];
    for (const q of current.queries) {
      q.hits.forEach((h, k) => {
        const ref = current.references[h.reference];
        rows.push([q.id, k + 1, h.reference, h.identity, h.coverage, h.strand, h.query_start, h.query_end,
          h.reference_start, h.reference_end, h.cigar, ref.taxon, ...ref.lineage, ref.occurrences.length]);
      });
    }
  } else {
    rows = [["query", "reference", "identity", "taxon", "latitude", "longitude", "fields"]];
    for (const q of current.queries) {
      for (const h of q.hits) {
        const ref = current.references[h.reference];
        for (const o of ref.occurrences) {
          rows.push([q.id, h.reference, h.identity, deepest(ref.lineage)[0], o.lat, o.lon,
            o.fields.map((f) => `${f.name}=${f.value}`).join("; ")]);
        }
      }
    }
  }
  const blob = new Blob([tsv(rows)], { type: "text/tab-separated-values" });
  const a = el("a", { href: URL.createObjectURL(blob), download: `asvatlas_${kind}.tsv` });
  document.body.append(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

// Before a search the page shows only the form; after it, the map and the
// list of ASVs ("Edit" goes back to the form).
function showView(results) {
  $("search-form").hidden = results;
  $("results-view").hidden = !results;
  if (results) {
    initMap();
    if (map) map.invalidateSize();
  } else {
    $("sequences").focus();
  }
}

let searching = false;
let full = null;      // the search at the lowest level, of which `current` is a view
let levelToken = 0;   // the last level asked for (older answers are ignored)

const lowestLevel = () => Math.min(...((dataset && dataset.identity_levels) || [100]));

async function postSearch(body) {
  const res = await fetch("api/search", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || res.statusText);
  return data;
}

// The results at similarity `v`, from the search at the lowest level: the
// matches of each ASV are filtered in the page. An ASV that reached the
// largest number of matches may lack some at `v` (the search stopped there):
// those ASVs alone are searched again at `v`, so the results are those of a
// search at `v`.
async function levelResults(v) {
  if (full.levels.has(v)) return full.levels.get(v);
  const cap = full.params.max_hits;
  const redo = v > full.params.identity ? full.queries.filter((q) => q.hits.length >= cap) : [];
  let again = new Map(), refs = full.references;
  if (redo.length) {
    const seqs = new Map(queryRecords(full.sequences).map((r) => [r.id, r.seq]));
    const data = await postSearch({ ...full.request, identity: v,
      sequences: redo.map((q) => `>${q.id}\n${seqs.get(q.id)}`).join("\n") });
    again = new Map(data.queries.map((q) => [q.id, q]));
    refs = { ...full.references, ...data.references };
  }
  const queries = full.queries.map((q) => {
    if (again.has(q.id)) return again.get(q.id);
    const hits = q.hits.filter((h) => h.identity >= v - 1e-9);
    return { ...q, hits, best_taxa: hits.length ? q.best_taxa : [] };
  });
  const out = { ...full, queries, references: refs, params: { ...full.params, identity: v } };
  full.levels.set(v, out);
  return out;
}

// Shows the results at similarity `v`; a search of the server is needed only
// for ASVs at the limit of matches (see levelResults).
async function showLevel(v, first = false) {
  const token = ++levelToken;
  const slow = !full.levels.has(v) && v > full.params.identity &&
    full.queries.some((q) => q.hits.length >= full.params.max_hits);
  if (slow && !first) {
    $("map-summary").className = "busy";
    $("map-summary").textContent = `Searching at ${levelText(v)}…`;
    $("map-wrap").classList.add("busy");
  }
  try {
    const out = await levelResults(v);
    if (token !== levelToken) return;
    current = out;
    selected = null;
    $("show-all").hidden = true;
    showView(true);
    renderResults();
    drawMap();
    updateSummary();
  } catch (e) {
    if (token !== levelToken) return;
    if (first) throw e;
    $("map-summary").className = "error";
    $("map-summary").textContent = e.message;
  } finally {
    if (token === levelToken) $("map-wrap").classList.remove("busy");
  }
}

// Searches the ASVs of the text area at the chosen similarity. Errors are
// shown under the button of the form, or above the map in the results.
async function search() {
  if (searching) return;
  const fromStart = !$("search-form").hidden;
  let coords;
  try {
    coords = checkQueries($("sequences").value);
  } catch (e) {
    if (fromStart || !current) setStatus(e.message, true);
    return;
  }
  // One search at the lowest similarity of the page: every higher level is a
  // filter of it, applied in the page (see showLevel).
  const body = {
    sequences: $("sequences").value,
    identity: lowestLevel(),
    query_cov: Number($("coverage").value) / 100,
    max_hits: Number($("max-hits").value),
    strand: $("both-strands").checked ? "both" : "plus",
  };
  searching = true;
  $("submit").disabled = true;
  if (fromStart) {
    setStatus("Searching the reference…");
  } else {
    $("map-summary").className = "busy";
    $("map-summary").textContent = "Searching…";
    $("map-wrap").classList.add("busy");
  }
  try {
    const data = await postSearch(body);
    // The countries need the world map; the results are shown without it.
    await (worldData || Promise.resolve()).catch(() => null);
    full = data;
    full.request = body;
    full.sequences = body.sequences;   // for the networks of the species found
    full.coords = coords;              // sampling places, kept in the page
    full.tableLoaded = coordTable.size > 0;
    full.levels = new Map();           // results of each level, once computed
    selected = null;
    $("show-all").hidden = true;
    setStatus("");
    await showLevel(chosenIdentity(), true);
  } catch (e) {
    if (fromStart || !current) {
      setStatus(e.message, true);
    } else {
      $("map-summary").className = "error";
      $("map-summary").textContent = e.message;
    }
  } finally {
    searching = false;
    $("submit").disabled = false;
    $("map-wrap").classList.remove("busy");
  }
}

$("sequences").addEventListener("input", () => { countSequences(); coordMatch(); });

// How many of the ASVs in the text area have coordinates in the table.
function coordMatch() {
  if (!coordTable.size) return;
  const records = queryRecords($("sequences").value);
  const n = records.filter((r) => coordTable.has(r.id)).length;
  $("coord-match").textContent = records.length
    ? `${fmt(n)} of your ${plural(records.length, "ASV")} have coordinates in the table.` : "";
}
$("example").addEventListener("click", () => {
  if (!dataset || !dataset.example) return;
  $("sequences").value = dataset.example;
  countSequences();
  coordMatch();
  updateClear();
  setStatus("Example loaded: press “Find in the world”.");
});
$("edit").addEventListener("click", () => showView(false));

// Empties the form: sequences, table of coordinates and messages.
function clearInput() {
  $("sequences").value = "";
  countSequences();
  coordTable = new Map();
  $("coord-count").textContent = "No table loaded";
  $("coord-match").textContent = "";
  setStatus("");
  updateClear();
}

// A new search: the form empty, the results and the network of the previous
// ASVs gone, the similarity back to identical records.
function newSearch() {
  clearInput();
  current = null;
  full = null;
  levelToken++;
  selected = null;
  placeBounds = null;
  if (layer) layer.clearLayers();
  $("asv-list").replaceChildren();
  $("show-all").hidden = true;
  setIdentity(100);
  resetNetwork();
  showView(false);
}

// "Clear" is offered once there is something to clear.
function updateClear() {
  $("clear").hidden = !$("sequences").value.trim() && !coordTable.size;
}

$("new-search").addEventListener("click", newSearch);
$("clear").addEventListener("click", () => { clearInput(); $("sequences").focus(); });
$("sequences").addEventListener("input", updateClear);
$("coord-file").addEventListener("change", (ev) => {
  const f = ev.target.files[0];
  ev.target.value = "";
  if (!f) return;
  f.text().then((t) => {
    try {
      coordTable = readCoordTable(t);
      const nplaces = [...coordTable.values()].reduce((a, l) => a + l.length, 0);
      $("coord-count").textContent = `${f.name}: ${plural(nplaces, "sampling place")} of ${plural(coordTable.size, "ASV")}`;
      coordMatch();
      updateClear();
      setStatus("");
    } catch (e) {
      coordTable = new Map();
      $("coord-count").textContent = "No table loaded";
      $("coord-match").textContent = "";
      setStatus(e.message, true);
    }
  });
});
$("file").addEventListener("change", (ev) => {
  const f = ev.target.files[0];
  if (f) readFile(f);
  ev.target.value = "";
});
// A FASTA file can also be dropped on the text area.
const drop = $("drop");
drop.addEventListener("dragover", (ev) => { ev.preventDefault(); drop.classList.add("dragging"); });
drop.addEventListener("dragleave", () => drop.classList.remove("dragging"));
drop.addEventListener("drop", (ev) => {
  ev.preventDefault();
  drop.classList.remove("dragging");
  const f = ev.dataTransfer.files[0];
  if (f) readFile(f);
});
$("show-all").addEventListener("click", () => {
  const open = document.querySelector("#asv-list li.open");
  if (open) toggle(current.queries.find((q) => q.id === open.dataset.id), open);
});
document.querySelectorAll("#download-menu button[data-table]").forEach((b) =>
  b.addEventListener("click", () => {
    $("download-menu").open = false;
    if (current) download(b.dataset.table);
  }));
// The download menu closes when clicking elsewhere.
document.addEventListener("click", (ev) => {
  if (!$("download-menu").contains(ev.target)) $("download-menu").open = false;
});
$("search-form").addEventListener("submit", (ev) => {
  ev.preventDefault();
  if (!countSequences()) {
    setStatus("Paste or load at least one sequence.", true);
    return;
  }
  search();
});

// The world map is fetched at once: it also gives the countries of the places.
function loadWorld() {
  return fetch("world.geojson").then((res) => res.json())
    .then((data) => { buildCountryIndex(data); return data; });
}
worldData = loadWorld();
loadDataset();
