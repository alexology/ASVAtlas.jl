"use strict";

// Haplotype network tab: a map of the places, with the network in a small
// translucent panel at its top right. A slider groups the haplotypes within a
// number of mutations of a centre (by default the haplotype of the user's
// ASV); places and nodes are coloured by group. Uses the helpers of app.js
// (el, $, fmt, plural, taxonNode, makeMap, tsv, current, dataset).

// Two colours, checked for the dark map (all pairs, ΔE 31.8 normal vision,
// 26.8 colour-blind, ≥ 3:1 against land and sea): warm for the group of the
// centre, cool for the others.
const IN_COLOR = "#d95926";
const OUT_COLOR = "#3987e5";
const INFERRED_COLOR = "#5d6b78";
const QUERY_RING = "#3cc8f0";
const SVG_NS = "http://www.w3.org/2000/svg";

let netMap = null;
let netLayer = null;
let network = null;       // last network response
let centre = 1;           // centre haplotype (1-based)
let distance = 0;         // slider value
let wantedCentre = null;  // query ASV to centre on after the next build
let baseBox = null;       // viewBox that fits the whole network
let box = null;           // current viewBox [x, y, w, h]

function svgEl(tag, attrs, ...children) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs || {})) node.setAttribute(k, v);
  for (const c of children) node.append(c instanceof Node ? c : String(c));
  return node;
}

function netStatus(text, error) {
  $("net-status").textContent = text;
  $("net-status").className = error ? "error" : "";
}

function netNote(text) {
  $("net-map-note").textContent = text;
  $("net-map-note").hidden = !text;
}

// After a new search the network of the previous ASVs is removed.
function resetNetwork() {
  network = null;
  baseBox = box = null;
  if (netLayer) netLayer.clearLayers();
  $("net-svg").replaceChildren();
  for (const id of ["net-overlay", "net-legend", "net-group", "net-results"]) $(id).hidden = true;
  $("net-title").textContent = "Haplotypes on the map";
  netStatus("");
  netNote("Choose a taxon and build its haplotype network.");
  fillSpecies();
}

// ---- Tabs -------------------------------------------------------------------

function showTab(name) {
  for (const t of ["search", "network", "about"]) {
    $("tab-" + t).hidden = t !== name;
    $("tab-btn-" + t).setAttribute("aria-selected", String(t === name));
  }
  window.scrollTo(0, 0);
  if (name === "network") {
    if (!netMap && typeof L !== "undefined") {
      ({ map: netMap, layer: netLayer } = makeMap("net-map", netNote, "topleft"));
    }
    if (netMap) setTimeout(() => netMap.invalidateSize(), 0);
    fillSpecies();
  } else if (name === "about") {
    fillAboutMenus();
  } else if (map) {
    setTimeout(() => map.invalidateSize(), 0);
  }
}

// ---- Species of the user's ASVs -------------------------------------------------

// Species matched by the ASVs of the last search, with those ASVs.
function speciesOfQueries() {
  const out = new Map();
  if (!current) return out;
  const k = current.ranks.indexOf("species");
  if (k < 0) return out;
  for (const q of current.queries) {
    for (const h of q.hits) {
      const name = current.references[h.reference].lineage[k];
      if (!name) continue;
      if (!out.has(name)) out.set(name, new Set());
      out.get(name).add(q.id);
    }
  }
  return out;
}

function fillSpecies() {
  const sel = $("net-species");
  const before = sel.value;
  const species = speciesOfQueries();
  sel.replaceChildren();
  for (const [name, ids] of [...species].sort((a, b) => b[1].size - a[1].size || a[0].localeCompare(b[0]))) {
    sel.append(el("option", { value: name }, `${name} · ${ids.size} of your ASV${ids.size === 1 ? "" : "s"}`));
  }
  const ready = species.size > 0 && dataset && dataset.networks;
  sel.disabled = !ready;
  $("net-submit").disabled = !ready;
  if ([...sel.options].some((o) => o.value === before)) sel.value = before;
  $("net-species-note").textContent = !(dataset && dataset.networks)
    ? "The reference has no species names: networks are not available."
    : !current ? "Search your ASVs first (tab “ASV search”): networks are built for the species they match, with your ASVs placed in them."
    : species.size ? `Species matched by your ASVs at ≥ ${current.params.identity}% similarity. Your ASVs matching the species are placed in its network.`
    : `None of your ASVs matches a species at ≥ ${current.params.identity}% similarity.`;
}

// Opens the tab on a species, from the search results, centred on the query.
function openNetwork(species, queryId) {
  showTab("network");
  $("net-species").value = species;
  wantedCentre = queryId;
  buildNetwork();
}

// ---- Building ------------------------------------------------------------------

async function buildNetwork(ev) {
  if (ev) ev.preventDefault();
  const name = $("net-species").value;
  const ids = speciesOfQueries().get(name);
  if (!name || !ids) {
    netStatus("Choose a species matched by your ASVs.", true);
    return;
  }
  const body = {
    taxon: name,
    method: document.querySelector('input[name="net-method"]:checked').value,
    epsilon: Number($("net-epsilon").value) || 0,
    identity: current.params.identity,
    include: [...ids],
    sequences: current.sequences,
  };
  $("net-submit").disabled = true;
  netStatus("Aligning the sequences and building the network…");
  const t0 = performance.now();
  try {
    const res = await fetch("api/network", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || res.statusText);
    network = data;
    fillCentres();
    drawNetwork();
    renderHaplotypes();
    distance = 0;
    setCentre(initialCentre(), true);
    const parts = data.components > 1 ? ` The network has ${data.components} separate parts (TCS connection limit).` : "";
    const rev = data.reversed.length
      ? ` Reverse-complemented to match the others: ${short(data.reversed, 5)}.` : "";
    const out = data.excluded.length
      ? ` Left out (no 8-mer shared with the others): ${short(data.excluded, 5)}.` : "";
    netStatus(`Built in ${((performance.now() - t0) / 1000).toFixed(1)} s.${parts}${rev}${out}`);
  } catch (e) {
    netStatus(e.message, true);
  } finally {
    $("net-submit").disabled = false;
    wantedCentre = null;
  }
}

// ---- Centre and grouping ------------------------------------------------------------

// The user's ASVs come first: the haplotype of an ASV is the natural centre.
function queryHaplotypes() {
  const out = [];
  network.haplotypes.forEach((hap, i) => hap.queries.forEach((q) => out.push([q, i + 1])));
  return out;
}

// The ASV the tab was opened for, else the first of the user's ASVs, else H1.
function initialCentre() {
  const sel = $("net-centre");
  const wanted = [...sel.options].find((o) => o.value === `q:${wantedCentre}`);
  sel.value = wanted ? wanted.value : sel.options[0].value;
  return Number(sel.selectedOptions[0].dataset.h);
}

function fillCentres() {
  const sel = $("net-centre");
  sel.replaceChildren();
  const qs = queryHaplotypes();
  if (qs.length) {
    const g = el("optgroup", { label: "Your ASVs" });
    for (const [q, h] of qs) g.append(el("option", { value: `q:${q}`, "data-h": h }, `${q} (${network.haplotypes[h - 1].name})`));
    sel.append(g);
  }
  const g = el("optgroup", { label: "Haplotypes" });
  network.haplotypes.forEach((hap, i) => {
    const t = hap.taxa.length ? ` · ${hap.taxa.map((x) => x.name).join(", ")}` : "";
    g.append(el("option", { value: `h:${i + 1}`, "data-h": i + 1 }, `${hap.name} · ${plural(hap.occurrences, "occurrence")}${t}`));
  });
  sel.append(g);
}

const distFrom = (h) => network.distances[centre - 1][h - 1];
const inGroup = (h) => h >= 1 && distFrom(h) <= distance;

// Sets the centre (and the slider range), then redraws everything.
function setCentre(h, fit) {
  centre = h;
  const sel = $("net-centre");
  if (Number(sel.selectedOptions[0]?.dataset.h) !== h) sel.value = `h:${h}`;
  const max = Math.max(1, ...network.distances[h - 1]);
  const slider = $("net-distance");
  slider.max = max;
  distance = Math.min(distance, max);
  slider.value = distance;
  $("net-distance-max").textContent = max;
  $("net-group").hidden = false;
  regroup(fit);
}

function regroup(fit) {
  $("net-distance-value").textContent = distance;
  const n = network.haplotypes.length;
  const members = [...Array(n).keys()].map((i) => i + 1).filter(inGroup);
  const places = network.sites.filter((s) => s.counts.some((c) => inGroup(c.haplotype))).length;
  const occ = members.reduce((a, h) => a + network.haplotypes[h - 1].occurrences, 0);
  const centreName = $("net-centre").selectedOptions[0]?.value.startsWith("q:")
    ? $("net-centre").value.slice(2) : network.haplotypes[centre - 1].name;
  $("net-group-note").replaceChildren(
    el("span", { class: "group-dot in" }),
    `${plural(members.length, "haplotype")} within ${plural(distance, "mutation")} of ${centreName}: ` +
    `${plural(occ, "occurrence")} at ${plural(places, "place")}.`);
  $("net-kpi-in").textContent = `${fmt(members.length)} / ${fmt(n)}`;
  $("net-kpi-sites").textContent = `${fmt(places)} / ${fmt(network.sites.length)}`;
  colourNetwork();
  colourTable();
  drawLegend(centreName);
  drawSites(fit);
}

// ---- Network drawing -----------------------------------------------------------------

// Circles and labels have a constant size on screen (in pixels): zooming in
// spreads the nodes apart instead of enlarging them. Sizes are set again
// whenever the view or the panel changes (resizeNodes).
function drawNetwork() {
  const svg = $("net-svg");
  svg.replaceChildren();
  // Coordinates are in mutational steps.
  const nodes = network.nodes;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const n of nodes) {
    x0 = Math.min(x0, n.x); x1 = Math.max(x1, n.x);
    y0 = Math.min(y0, n.y); y1 = Math.max(y1, n.y);
  }
  const pad = Math.max(1.5, 0.08 * Math.max(x1 - x0, y1 - y0));
  const w = Math.max(x1 - x0 + 2 * pad, 8), h = Math.max(y1 - y0 + 2 * pad, 5);
  baseBox = [(x0 + x1) / 2 - w / 2, (y0 + y1) / 2 - h / 2, w, h];

  const edges = svgEl("g", { class: "net-edges" });
  for (const e of network.edges) {
    const a = nodes[e.source - 1], b = nodes[e.target - 1];
    edges.append(svgEl("line", { x1: a.x, y1: a.y, x2: b.x, y2: b.y }));
    if (e.steps > 1) {
      edges.append(svgEl("text", { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, class: "net-steps" }, e.steps));
    }
  }
  svg.append(edges);

  const g = svgEl("g", { class: "net-nodes" });
  nodes.forEach((n) => {
    const hp = n.haplotype;
    // Radius on screen: 2.5 px for inferred nodes, 4–14 px with the frequency.
    const px = hp === 0 ? 2.5 : 4 + 10 * Math.min(1, n.r / 0.8);
    if (hp === 0) {
      g.append(svgEl("circle", { cx: n.x, cy: n.y, fill: INFERRED_COLOR, class: "inferred", "data-px": px },
        svgEl("title", {}, `${n.name} (inferred)`)));
      return;
    }
    const hap = network.haplotypes[hp - 1];
    const node = svgEl("g", { class: "hap", "data-h": hp, tabindex: 0, role: "button" });
    node.append(svgEl("circle", { cx: n.x, cy: n.y, class: "centre-ring", "data-px": px + 5 }));
    if (hap.queries.length) {
      node.append(svgEl("circle", { cx: n.x, cy: n.y, class: "query-ring", stroke: QUERY_RING, "data-px": px + 2.5 }));
    }
    node.append(svgEl("circle", { cx: n.x, cy: n.y, class: "hap-dot", "data-px": px }));
    node.append(svgEl("text", { x: n.x, "data-y": n.y, "data-px": px, class: "hap-label" }, hap.name));
    node.append(svgEl("title", {}, hapTitle(hap)));
    const pick = () => setCentre(hp, false);
    node.addEventListener("click", pick);
    node.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter" || ev.key === " ") { ev.preventDefault(); pick(); }
    });
    g.append(node);
  });
  svg.append(g);
  $("net-overlay").hidden = false;
  $("net-overlay-title").textContent = `${network.method.toUpperCase()} · ${network.taxon}`;
  $("net-title").textContent = `Haplotypes of ${network.taxon} on the map`;
  box = baseBox.slice();
  setBox();
}

// Network units per screen pixel for the current view.
function unitsPerPixel() {
  const rect = $("net-svg").getBoundingClientRect();
  if (!rect.width || !rect.height) return box[2] / 340;
  return Math.max(box[2] / rect.width, box[3] / rect.height);
}

function resizeNodes() {
  if (!box) return;
  const u = unitsPerPixel();
  for (const c of document.querySelectorAll("#net-svg circle")) c.setAttribute("r", Number(c.dataset.px) * u);
  for (const t of document.querySelectorAll("#net-svg .hap-label")) {
    t.setAttribute("font-size", 11 * u);
    t.setAttribute("y", Number(t.dataset.y) + (Number(t.dataset.px) + 8) * u);
  }
  for (const t of document.querySelectorAll("#net-svg .net-steps")) t.setAttribute("font-size", 9.5 * u);
}

function colourNetwork() {
  const many = network.haplotypes.length > 40;
  document.querySelectorAll("#net-svg .hap").forEach((g) => {
    const h = Number(g.dataset.h);
    const inside = inGroup(h);
    g.classList.toggle("in", inside);
    g.classList.toggle("centre", h === centre);
    g.querySelector(".hap-dot").setAttribute("fill", inside ? IN_COLOR : OUT_COLOR);
    // In large networks only the group is labelled.
    g.querySelector(".hap-label").classList.toggle("hidden", many && !inside);
  });
}

function hapTitle(hap) {
  const lines = [`${hap.name}: ${plural(hap.occurrences, "occurrence")}`];
  if (hap.taxa.length) lines.push(hap.taxa.map((t) => t.name).join(", "));
  lines.push("Reference ASVs: " + (short(hap.references, 6) || "none"));
  if (hap.queries.length) {
    // With a table of coordinates, the ASVs without a place in it are named.
    const none = current && current.tableLoaded ? hap.queries.filter((q) => !samplesOf(q).length) : [];
    lines.push("Your ASVs: " + hap.queries.join(", "));
    if (none.length) lines.push("Not on the map (no sampling place in your table): " + none.join(", "));
  }
  lines.push("Click to make it the centre");
  return lines.join("\n");
}

function drawLegend(centreName) {
  const legend = $("net-legend");
  const dot = (cls) => el("i", { class: "group-dot " + cls });
  legend.replaceChildren(
    el("span", {}, dot("in"), `≤ ${plural(distance, "mutation")} from ${centreName}`),
    el("span", {}, dot("out"), "Farther"));
  if (network.haplotypes.some((hp) => hp.queries.length)) {
    legend.append(el("span", {}, el("i", { class: "ring" }), "Your ASVs (network)"));
  }
  if (network.haplotypes.some((hp) => hp.queries.some((q) => samplesOf(q).length))) {
    // The stars of the legend are those of the map: group colour, gold border.
    const star = (fill) => {
      const i = el("i", { class: "legend-star", "aria-hidden": "true" });
      i.innerHTML = starSvg(18, [fill]);   // fixed markup: numbers and constant colours only
      return i;
    };
    legend.append(el("span", {}, star(IN_COLOR), star(OUT_COLOR), "Your samples, by group"));
  }
  legend.hidden = false;
}

// Frames the haplotypes of the group (and the centre), keeping the aspect of
// the whole network.
function zoomGroup() {
  if (!network) return;
  const nodes = document.querySelectorAll("#net-svg .hap.in .hap-dot, #net-svg .hap.centre .hap-dot");
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const c of nodes) {
    const x = Number(c.getAttribute("cx")), y = Number(c.getAttribute("cy")), r = Number(c.getAttribute("r"));
    x0 = Math.min(x0, x - r); x1 = Math.max(x1, x + r); y0 = Math.min(y0, y - r); y1 = Math.max(y1, y + r);
  }
  const pad = baseBox[2] / 25;
  const scale = Math.max((x1 - x0 + 2 * pad) / baseBox[2], (y1 - y0 + 2 * pad) / baseBox[3], 0.08);
  const w = baseBox[2] * scale, h = baseBox[3] * scale;
  box = [(x0 + x1) / 2 - w / 2, (y0 + y1) / 2 - h / 2, w, h];
  setBox();
}

function setBox() {
  $("net-svg").setAttribute("viewBox", box.join(" "));
  resizeNodes();
}

// Zoom with the wheel around the pointer, pan by dragging.
function initNetworkView() {
  const svg = $("net-svg");
  const toSvg = (ev) => {
    const rect = svg.getBoundingClientRect();
    const s = Math.max(box[2] / rect.width, box[3] / rect.height);
    const ox = box[0] + (box[2] - rect.width * s) / 2, oy = box[1] + (box[3] - rect.height * s) / 2;
    return [ox + (ev.clientX - rect.left) * s, oy + (ev.clientY - rect.top) * s, s];
  };
  svg.addEventListener("wheel", (ev) => {
    if (!box) return;
    ev.preventDefault();
    const [px, py] = toSvg(ev);
    const f = ev.deltaY > 0 ? 1.15 : 1 / 1.15;
    const w = Math.min(Math.max(box[2] * f, baseBox[2] / 20), baseBox[2] * 4);
    const k = w / box[2];
    box = [px - (px - box[0]) * k, py - (py - box[1]) * k, w, box[3] * k];
    setBox();
  }, { passive: false });
  let drag = null;
  svg.addEventListener("pointerdown", (ev) => {
    if (!box || ev.target.closest(".hap")) return;
    drag = { x: ev.clientX, y: ev.clientY, box: box.slice(), s: toSvg(ev)[2] };
    svg.setPointerCapture(ev.pointerId);
  });
  svg.addEventListener("pointermove", (ev) => {
    if (!drag) return;
    box = [drag.box[0] - (ev.clientX - drag.x) * drag.s, drag.box[1] - (ev.clientY - drag.y) * drag.s,
      drag.box[2], drag.box[3]];
    setBox();
  });
  svg.addEventListener("pointerup", () => { drag = null; });
  // The map does not react to the wheel and drags meant for the network.
  if (typeof L !== "undefined") {
    L.DomEvent.disableClickPropagation($("net-overlay"));
    L.DomEvent.disableScrollPropagation($("net-overlay"));
  }
  $("net-fit").addEventListener("click", () => { if (baseBox) { box = baseBox.slice(); setBox(); } });
  $("net-zoom-group").addEventListener("click", zoomGroup);
  // Sizes in pixels follow the panel when it is enlarged or the window resized.
  if (typeof ResizeObserver !== "undefined") new ResizeObserver(resizeNodes).observe(svg);
  $("net-expand").addEventListener("click", () => {
    const big = $("net-overlay").classList.toggle("expanded");
    $("net-expand").setAttribute("aria-pressed", String(big));
    $("net-expand").title = big ? "Shrink the network" : "Enlarge the network";
  });
}

// ---- Map --------------------------------------------------------------------------

// Pie chart of the occurrences of the group and of the others at a place, as
// an SVG string made only of numbers and the two fixed colours.
function pieSvg(nIn, nOut, size) {
  const r = size / 2 - 1.5, cx = size / 2, cy = size / 2;
  const total = nIn + nOut;
  let paths;
  if (nIn === 0 || nOut === 0) {
    paths = `<circle cx="${cx}" cy="${cy}" r="${r}" fill="${nIn ? IN_COLOR : OUT_COLOR}"/>`;
  } else {
    const a0 = -Math.PI / 2, a1 = a0 + 2 * Math.PI * nIn / total;
    const large = a1 - a0 > Math.PI ? 1 : 0;
    const p = (a) => `${(cx + r * Math.cos(a)).toFixed(2)} ${(cy + r * Math.sin(a)).toFixed(2)}`;
    paths = `<path d="M${cx} ${cy}L${p(a0)}A${r} ${r} 0 ${large} 1 ${p(a1)}Z" fill="${IN_COLOR}"/>` +
      `<path d="M${cx} ${cy}L${p(a1)}A${r} ${r} 0 ${1 - large} 1 ${p(a0 + 2 * Math.PI)}Z" fill="${OUT_COLOR}"/>`;
  }
  return `<svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">` +
    `${paths}<circle cx="${cx}" cy="${cy}" r="${r}" fill="none" stroke="#03101b" stroke-width="1.5"/></svg>`;
}

function drawSites(fit) {
  if (!netMap || !network) return;
  netLayer.clearLayers();
  // Sampling locations of the user's ASVs in the network (from the table of
  // coordinates of the search), with their haplotype and group.
  const hapOf = new Map();
  network.haplotypes.forEach((hap, i) => hap.queries.forEach((q) => hapOf.set(q, i + 1)));
  const samples = drawSamples(netLayer, [...hapOf.keys()], (ids) => ids.map((id) => {
    const h = hapOf.get(id);
    return `${id}: ${network.haplotypes[h - 1].name}, ${inGroup(h) ? "in the group" : "farther"} ` +
      `(${plural(distFrom(h), "mutation")})`;
  }), (ids) => {
    // Filled with the colour of the group; split when the ASVs of a place differ.
    const inside = ids.map((id) => inGroup(hapOf.get(id)));
    return inside.every((x) => x) ? [IN_COLOR] : inside.some((x) => x) ? [IN_COLOR, OUT_COLOR] : [OUT_COLOR];
  });
  const sites = network.sites;
  if (!sites.length) {
    netNote("None of these haplotypes has occurrences to show on the map.");
    if (fit && samples.length) netMap.setView(samples[0], 5);
    return;
  }
  netNote("");
  const bounds = samples.slice();
  // Places of the group are drawn last, on top.
  const ranked = sites.map((s) => {
    let nIn = 0, nOut = 0;
    for (const c of s.counts) inGroup(c.haplotype) ? (nIn += c.n) : (nOut += c.n);
    return { s, nIn, nOut };
  }).sort((a, b) => (a.nIn > 0) - (b.nIn > 0));
  for (const { s, nIn, nOut } of ranked) {
    const total = nIn + nOut;
    const size = Math.round(Math.min(44, 14 + 6 * Math.sqrt(total)));
    const icon = L.divIcon({ className: nIn ? "pie-icon in" : "pie-icon out", html: pieSvg(nIn, nOut, size),
      iconSize: [size, size], iconAnchor: [size / 2, size / 2] });
    const m = L.marker([s.lat, s.lon], { icon, keyboard: true, zIndexOffset: nIn ? 1000 : 0 });
    const tip = nIn && nOut ? `${nIn} in the group, ${nOut} farther` : nIn ? `${plural(nIn, "occurrence")} in the group`
      : `${plural(nOut, "occurrence")} farther`;
    m.bindTooltip(el("span", {}, tip), { direction: "top", offset: [0, -size / 2] });
    m.bindPopup(() => sitePopup(s), { maxWidth: 360 });
    m.addTo(netLayer);
    bounds.push([s.lat, s.lon]);
  }
  if (!fit) return;
  if (bounds.length === 1) netMap.setView(bounds[0], 8);
  else netMap.fitBounds(bounds, { padding: [40, 40], maxZoom: 10 });
}

function sitePopup(s) {
  const box = el("div", {}, el("div", { class: "popup-coords" }, `${s.lat.toFixed(5)}, ${s.lon.toFixed(5)}`));
  const ul = el("ul");
  for (const c of [...s.counts].sort((a, b) => distFrom(a.haplotype) - distFrom(b.haplotype))) {
    const hap = network.haplotypes[c.haplotype - 1];
    ul.append(el("li", {},
      el("span", {}, el("i", { class: "group-dot " + (inGroup(c.haplotype) ? "in" : "out") }), hap.name,
        el("span", { class: "muted" }, ` · ${plural(distFrom(c.haplotype), "mutation")}`)),
      el("span", {}, `×${c.n}`)));
  }
  box.append(el("div", { class: "popup-ref" }, ul));
  for (const o of s.occurrences.slice(0, 20)) {
    const fields = o.fields.map((f) => `${f.name}: ${f.value}`).join(" · ");
    box.append(el("div", { class: "fields" }, o.reference + (fields ? " — " + fields : "")));
    const note = locationNote(o.fields);
    if (note) box.append(el("div", { class: "fields approx" }, note));
  }
  if (s.occurrences.length > 20) box.append(el("div", { class: "fields" }, `… ${s.occurrences.length - 20} more`));
  return box;
}

// ---- Table ------------------------------------------------------------------------

function hapTaxa(hap) {
  const span = el("span");
  hap.taxa.forEach((t, k) => {
    if (k) span.append(", ");
    span.append(taxonNode(t.name, rankIndexOf(t.rank)));
  });
  return span;
}

function rankIndexOf(rank) {
  return ((dataset && dataset.ranks) || []).indexOf(rank);
}

function short(list, n) {
  return list.length <= n ? list.join(", ") : `${list.slice(0, n).join(", ")} … (+${list.length - n})`;
}

function placesOf() {
  const places = new Map();
  for (const s of network.sites) for (const c of s.counts) places.set(c.haplotype, (places.get(c.haplotype) || 0) + 1);
  return places;
}

function renderHaplotypes() {
  const places = placesOf();
  const tbody = $("net-table").querySelector("tbody");
  tbody.replaceChildren();
  network.haplotypes.forEach((hap, i) => {
    const h = i + 1;
    const tr = el("tr", { class: "query", "data-h": h },
      el("td", { class: "asv-id" }, el("span", { class: "identity-cell" }, el("i", { class: "group-dot" }), hap.name)),
      el("td", { class: "num dist" }, ""),
      el("td", { class: "num" }, fmt(hap.occurrences)),
      el("td", { class: "num" }, fmt(places.get(h) || 0)),
      el("td", {}, hapTaxa(hap)),
      el("td", { class: "lineage" }, short(hap.references, 4)),
      el("td", {}, hap.queries.length ? short(hap.queries, 4) : "—"));
    tr.addEventListener("click", () => setCentre(h, false));
    tbody.append(tr);
  });
  $("net-kpi-haps").textContent = fmt(network.haplotypes.length);
  $("net-kpi-haps-label").textContent = `Haplotypes from ${plural(network.sequences, "sequence")}`;
  $("net-kpi-var").textContent = fmt(network.variable_sites);
  $("net-kpi-var-label").textContent = `Variable sites (of ${fmt(network.sites_used)} used)`;
  $("net-results").hidden = false;
}

// Group colour, distance and order by distance from the centre.
function colourTable() {
  const tbody = $("net-table").querySelector("tbody");
  const rows = [...tbody.children];
  for (const tr of rows) {
    const h = Number(tr.dataset.h);
    tr.querySelector(".group-dot").className = "group-dot " + (inGroup(h) ? "in" : "out");
    tr.querySelector(".dist").textContent = distFrom(h);
    tr.classList.toggle("selected", h === centre);
  }
  rows.sort((a, b) => distFrom(Number(a.dataset.h)) - distFrom(Number(b.dataset.h)) ||
    Number(a.dataset.h) - Number(b.dataset.h));
  tbody.append(...rows);
}

function downloadHaplotypes() {
  if (!network) return;
  const places = placesOf();
  const rows = [["haplotype", "mutations_from_centre", "in_group", "frequency", "occurrences", "places", "taxa",
    "reference_asvs", "query_asvs"]];
  network.haplotypes.forEach((hap, i) => rows.push([hap.name, distFrom(i + 1), inGroup(i + 1) ? "yes" : "no",
    hap.frequency, hap.occurrences, places.get(i + 1) || 0, hap.taxa.map((t) => t.name).join("; "),
    hap.references.join(";"), hap.queries.join(";")]));
  const blob = new Blob([tsv(rows)], { type: "text/tab-separated-values" });
  const a = el("a", { href: URL.createObjectURL(blob),
    download: `asvatlas_haplotypes_${network.method}_${network.haplotypes[centre - 1].name}_d${distance}.tsv` });
  document.body.append(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

// ---- Wiring -----------------------------------------------------------------------

$("tab-btn-search").addEventListener("click", () => showTab("search"));
$("tab-btn-network").addEventListener("click", () => showTab("network"));
$("tab-btn-about").addEventListener("click", () => showTab("about"));
$("network-form").addEventListener("submit", buildNetwork);
$("net-download").addEventListener("click", downloadHaplotypes);
$("net-centre").addEventListener("change", () => {
  if (network) setCentre(Number($("net-centre").selectedOptions[0].dataset.h), false);
});
$("net-distance").addEventListener("input", () => {
  if (!network) return;
  distance = Number($("net-distance").value);
  regroup(false);
});
initNetworkView();
