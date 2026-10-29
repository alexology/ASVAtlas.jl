"use strict";

// Tab "Reference data": what the reference is made of (marker, primers,
// studies, processing, origin of the coordinates), from the description
// given to load_reference(...; about) and the counts of the server. Uses the
// helpers of app.js (el, $, fmt, plural, dataset). All text is written with
// textContent; links only to http(s) addresses.

const COORD_KINDS = {
  sites: ["Sites", "Coordinates of the sampling sites, from the publication."],
  approximate: ["Approximate", "Sites placed at random around a point of the study area."],
  "study area": ["Study area", "All sites at one point of the study area."],
  none: ["None", "No coordinates: the records are not on the map."],
};

// References served by this page: one for now; the menus are ready for more
// (one per marker and primer pair).
const references = () => (dataset ? [dataset] : []);

function fillAboutMenus() {
  const markers = [...new Set(references().map((r) => r.marker || "Unknown marker"))];
  const ms = $("about-marker");
  ms.replaceChildren(...markers.map((m) => el("option", { value: m }, m)));
  fillAboutDatabases();
}

function fillAboutDatabases() {
  const m = $("about-marker").value;
  const ds = references().filter((r) => (r.marker || "Unknown marker") === m);
  $("about-db").replaceChildren(...ds.map((r) =>
    el("option", { value: r.name }, r.primers ? `${r.name} (${r.primers})` : r.name)));
  renderAbout();
}

function link(text, url) {
  return /^https?:\/\//.test(url || "")
    ? el("a", { href: url, target: "_blank", rel: "noopener noreferrer" }, text)
    : el("span", {}, text);
}

function renderAbout() {
  const ref = references().find((r) => r.name === $("about-db").value);
  const box = $("about-content");
  box.replaceChildren();
  if (!ref) return;
  const about = ref.about || {};
  const sources = Array.isArray(about.sources) ? about.sources : [];
  const counts = ref.sources || {};

  // Overview.
  const chips = el("div", { class: "about-chips" });
  for (const [k, v] of [["Marker", ref.marker], ["Primers", ref.primers]]) {
    if (v) chips.append(el("span", { class: "chip" }, el("b", {}, k), " ", v));
  }
  const countries = new Set(sources.map((s) => s.country).filter((x) => x));
  const figures = el("dl", { class: "about-figures" });
  for (const [k, v] of [["Reference ASVs", ref.references], ["Records", ref.occurrences],
    ["Studies", sources.length || null], ["Countries", countries.size || null]]) {
    if (v !== null && v !== undefined) figures.append(el("div", {}, el("dt", {}, k), el("dd", {}, fmt(v))));
  }
  box.append(el("section", { class: "card about-card" },
    el("h1", {}, ref.name), chips,
    about.description ? el("p", { class: "about-lead" }, about.description) : null,
    figures));

  // Studies.
  if (sources.length) {
    const body = el("tbody");
    for (const s of sources) {
      const c = counts[s.dataset] || {};
      const kind = COORD_KINDS[s.coordinates] || COORD_KINDS.sites;
      const study = el("td", {},
        el("div", { class: "study-label" }, link(s.label || s.dataset, s.url)),
        s.citation ? el("div", { class: "study-citation", title: s.citation }, s.citation) : null,
        el("div", { class: "study-id" }, s.dataset));
      body.append(el("tr", {},
        study,
        el("td", {}, s.country || "—", s.area ? el("div", { class: "muted" }, s.area) : null),
        el("td", {}, s.sample || "—"),
        el("td", {}, s.primers || "—"),
        el("td", { class: "num" }, c.asvs === undefined ? "0" : fmt(c.asvs)),
        el("td", { class: "num" }, c.records === undefined ? "0" : fmt(c.records)),
        el("td", {}, el("span", { class: `coord-kind kind-${(s.coordinates || "sites").replace(" ", "-")}`,
          title: kind[1] }, kind[0]))));
    }
    box.append(el("section", { class: "card" },
      el("div", { class: "card-head" }, el("h2", {}, "Studies"),
        el("span", { class: "muted" }, `${sources.length === 1 ? "1 study" : `${fmt(sources.length)} studies`}; ` +
          "the link of a study opens its raw data")),
      el("div", { class: "table-wrap" }, el("table", { class: "about-table" },
        el("thead", {}, el("tr", {}, ...[["Study"], ["Country"], ["Sample"], ["Primers"], ["ASVs", "num"],
          ["Records", "num"], ["Coordinates"]].map(([h, cl]) => el("th", cl ? { class: cl } : {}, h)))),
        body))));
  }

  // Processing and coordinates.
  const notes = el("div", { class: "about-notes" });
  if (Array.isArray(about.processing) && about.processing.length) {
    notes.append(el("section", { class: "card about-card" }, el("h2", {}, "How the data were processed"),
      el("ol", { class: "about-steps" }, ...about.processing.map((t) => el("li", {}, t)))));
  }
  if (about.coordinates) {
    const legend = el("dl", { class: "coord-legend" });
    for (const [k, [name, text]] of Object.entries(COORD_KINDS)) {
      legend.append(el("dt", {}, el("span", { class: `coord-kind kind-${k.replace(" ", "-")}` }, name)),
        el("dd", {}, text));
    }
    notes.append(el("section", { class: "card about-card" }, el("h2", {}, "Coordinates"),
      el("p", {}, about.coordinates), legend));
  }
  if (notes.childNodes.length) box.append(notes);
  if (!ref.about) {
    box.append(el("p", { class: "muted about-empty" }, "No description was given for this reference."));
  }
}

$("about-marker").addEventListener("change", fillAboutDatabases);
$("about-db").addEventListener("change", renderAbout);
