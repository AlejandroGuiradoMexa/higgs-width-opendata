// Geant4 calorimeter on the project page: a to-scale schematic of the simulated detector and an event display
// of single simulated showers. Shared by index.html and paper.html; colours come from each page's CSS variables.
// Geometry from g4calo/include/Geometry.hh and g4calo/src/Detector.cc.
(function () {
  const root = document.documentElement;
  const cssv = (...names) => { for (const n of names) { const v = getComputedStyle(root).getPropertyValue(n).trim(); if (v) return v; } return "#888"; };
  const pal = () => ({ ink: cssv("--ink", "--text"), muted: cssv("--ink-2", "--muted"), line: cssv("--hair", "--border"),
                       tint: cssv("--tint", "--surface-2"), sheet: cssv("--sheet", "--surface"), accent: cssv("--c-h", "--accent") });
  const NS = "http://www.w3.org/2000/svg";
  const el = (tag, attrs = {}, text) => { const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]); if (text !== undefined) e.textContent = text; return e; };

  // Layers of the readout, in the order of the output files (eta x phi windows centred on the impact point).
  const LAYERS = [
    { key: "ps", name: "Presampler", nx: 5, ny: 3, deta: 0.025, dphi: 0.1 },
    { key: "s1", name: "Strips (S1)", nx: 41, ny: 3, deta: 0.003125, dphi: 0.1 },
    { key: "s2", name: "Middle (S2)", nx: 11, ny: 11, deta: 0.025, dphi: 0.0245 },
    { key: "s3", name: "Back (S3)", nx: 5, ny: 11, deta: 0.05, dphi: 0.0245 },
    { key: "had", name: "Hadronic", nx: 3, ny: 3, deta: 0.1, dphi: 0.1 },
  ];

  // ---------------------------------------------------------------- detector schematic
  function detector(container) {
    const c = pal();
    const W = 1000, H = 330, y0 = 92, y1 = 212;           // block band
    const s = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Side view of the simulated calorimeter" });
    s.style.fontFamily = "inherit";
    const defs = el("defs");
    const hatch = el("pattern", { id: "g4hatch", width: 7, height: 7, patternUnits: "userSpaceOnUse", patternTransform: "rotate(45)" });
    hatch.appendChild(el("line", { x1: 0, y1: 0, x2: 0, y2: 7, stroke: c.muted, "stroke-width": 1.2 }));
    defs.appendChild(hatch); s.appendChild(defs);
    const T = (x, y, t, o = {}) => s.appendChild(el("text", Object.assign({ x, y, "font-size": 13, fill: c.ink, "text-anchor": "middle" }, o), t));
    const box = (x, w, fill) => s.appendChild(el("rect", { x, y: y0, width: w, height: y1 - y0, fill, stroke: c.ink, "stroke-width": 1 }));
    const periods = (x, w, n, col, sw) => { for (let i = 0; i < n; i++) { const xi = x + (i + 0.3) * w / n;
      s.appendChild(el("line", { x1: xi, y1: y0, x2: xi, y2: y1, stroke: col, "stroke-width": sw })); } };

    // Beam
    s.appendChild(el("line", { x1: 8, y1: (y0 + y1) / 2, x2: 58, y2: (y0 + y1) / 2, stroke: c.accent, "stroke-width": 2.5 }));
    s.appendChild(el("path", { d: `M58 ${(y0 + y1) / 2 - 6} L70 ${(y0 + y1) / 2} L58 ${(y0 + y1) / 2 + 6} Z`, fill: c.accent }));
    T(36, (y0 + y1) / 2 - 14, "beam", { "font-style": "italic" });
    T(36, (y0 + y1) / 2 + 26, "1.5 m", { "font-size": 11, fill: c.muted });

    // EM layers drawn to scale in X0 (17.9 px per X0); the others are not to scale.
    const px = 17.9, xS1 = 214, wS1 = 4.46 * px, wS2 = 16.05 * px, wS3 = 2.08 * px;
    box(82, 52, "url(#g4hatch)");                                     // upstream Al
    box(165, 14, c.tint);                                             // presampler (LAr)
    box(xS1, wS1 + wS2 + wS3, c.tint);
    periods(xS1, wS1, 15, c.ink, 1.1); periods(xS1 + wS1, wS2, 54, c.ink, 1.1); periods(xS1 + wS1 + wS2, wS3, 7, c.ink, 1.1);
    [xS1 + wS1, xS1 + wS1 + wS2].forEach(x => s.appendChild(el("line", { x1: x, y1: y0 - 8, x2: x, y2: y1 + 8, stroke: c.ink, "stroke-width": 1.6 })));
    const xH = 690, wH = 290;
    box(xH, wH, c.tint); periods(xH, wH, 20, c.muted, 5);
    // axis break on the hadronic section
    [y0, y1].forEach(y => s.appendChild(el("path", { d: `M${xH + wH / 2 - 9} ${y + 7} l6 -14 M${xH + wH / 2 + 1} ${y + 7} l6 -14`, stroke: c.ink, "stroke-width": 1.2 })));

    // Names (above)
    T(108, y0 - 30, "Upstream"); T(108, y0 - 14, "material", { "font-size": 12, fill: c.muted });
    T(172, y0 - 30, "PS");
    T(xS1 + wS1 / 2, y0 - 30, "S1"); T(xS1 + wS1 / 2, y0 - 14, "strips", { "font-size": 12, fill: c.muted });
    T(xS1 + wS1 + wS2 / 2, y0 - 30, "S2"); T(xS1 + wS1 + wS2 / 2, y0 - 14, "middle", { "font-size": 12, fill: c.muted });
    T(xS1 + wS1 + wS2 + wS3 / 2 + 6, y0 - 30, "S3"); T(xS1 + wS1 + wS2 + wS3 / 2 + 6, y0 - 14, "back", { "font-size": 12, fill: c.muted });
    T(xH + wH / 2, y0 - 30, "Hadronic calorimeter"); T(xH + wH / 2, y0 - 14, "not to scale", { "font-size": 12, fill: c.muted });

    // Details (below)
    const d = (x, lines) => lines.forEach((t, i) => T(x, y1 + 22 + 16 * i, t, { "font-size": 12, fill: i ? c.muted : c.ink }));
    d(108, ["Al, 0–6 X₀", "2 X₀ nominal", "ID, solenoid,", "cryostat"]);
    d(xS1 + wS1 / 2 - 4, ["15 × Pb/LAr", "4.5 X₀", "Δη×Δφ =", "0.0031×0.1"]);
    d(xS1 + wS1 + wS2 / 2, ["54 × (1.5 mm Pb + 4.2 mm LAr)", "16.1 X₀", "Δη×Δφ = 0.025×0.0245"]);
    d(xS1 + wS1 + wS2 + wS3 / 2 + 20, ["7 ×", "2.1 X₀", "0.05×", "0.0245"]);
    d(xH + wH / 2, ["20 × (40 mm Fe + 3 mm scint.)", "4.8 λ", "Δη×Δφ = 0.1×0.1"]);
    T(174, y1 + 22, "LAr", { "font-size": 12 }); T(174, y1 + 38, "11 mm", { "font-size": 12, fill: c.muted });

    // X0 scale under the EM layers
    const ya = H - 22;
    s.appendChild(el("line", { x1: xS1, y1: ya, x2: xS1 + 22.6 * px, y2: ya, stroke: c.ink, "stroke-width": 1 }));
    for (let t = 0; t <= 20; t += 5) { const x = xS1 + t * px;
      s.appendChild(el("line", { x1: x, y1: ya - 5, x2: x, y2: ya, stroke: c.ink, "stroke-width": 1 })); T(x, ya + 15, String(t), { "font-size": 11 }); }
    T(xS1 + 22.6 * px + 8, ya + 4, "X₀ (Pb/LAr, X₀ = 19.2 mm)", { "font-size": 11, "text-anchor": "start", fill: c.muted });
    container.replaceChildren(s);
  }

  // ---------------------------------------------------------------- event display
  const VIRIDIS = ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"];
  const lerp = (a, b, t) => Math.round(a + (b - a) * t);
  function viridis(t) {
    t = Math.max(0, Math.min(1, t)); const k = Math.min(3, Math.floor(t * 4)), u = t * 4 - k;
    const h = x => [1, 3, 5].map(i => parseInt(x.slice(i, i + 2), 16));
    const [a, b] = [h(VIRIDIS[k]), h(VIRIDIS[k + 1])];
    return `rgb(${lerp(a[0], b[0], u)},${lerp(a[1], b[1], u)},${lerp(a[2], b[2], u)})`;
  }
  const LABEL = { "e-": "e⁻", gamma: "γ", pi0: "π⁰", "pi-": "π⁻" };
  const LMIN = -1;   // 0.1 MeV: lower end of the colour scale

  function events(container, samples) {
    if (!samples || !samples.length) { container.textContent = "Event data not available."; return { draw() {} }; }
    let si = 0, ei = 0;
    container.classList.add("g4ev");
    const ctl = document.createElement("div"); ctl.className = "g4ev-ctl";
    const sel = document.createElement("select"); sel.setAttribute("aria-label", "Simulated sample");
    sel.innerHTML = samples.map((s, i) => `<option value="${i}">${LABEL[s.particle] || s.particle}, ${s.E} GeV, ${s.x0} X₀ upstream</option>`).join("");
    const seg = document.createElement("span"); seg.className = "seg";
    seg.innerHTML = samples[0].events.map((_, i) => `<button data-i="${i}"${i ? "" : ' class="on"'}>Event ${i + 1}</button>`).join("");
    ctl.append(sel, seg);
    const grid = document.createElement("div"); grid.className = "g4ev-grid";
    const foot = document.createElement("div"); foot.className = "g4ev-foot";
    container.replaceChildren(ctl, grid, foot);
    sel.onchange = () => { si = +sel.value; draw(); };
    seg.querySelectorAll("button").forEach(b => b.onclick = () => {
      seg.querySelectorAll("button").forEach(x => x.classList.remove("on")); b.classList.add("on"); ei = +b.dataset.i; draw(); });

    function draw() {
      const c = pal(), ev = samples[si].events[ei];
      const tot = LAYERS.reduce((a, L) => a + ev[L.key].flat().reduce((x, y) => x + y, 0), 0);
      const lmax = Math.log10(Math.max(1, ...LAYERS.map(L => Math.max(...ev[L.key].flat()))));
      const R = 0.16;   // eta and phi half-range of every map
      grid.replaceChildren(...LAYERS.map(L => {
        const fig = document.createElement("figure"); fig.className = "g4ev-map";
        const s = el("svg", { viewBox: `${-R} ${-R} ${2 * R} ${2 * R}`, role: "img", "aria-label": L.name });
        s.appendChild(el("rect", { x: -R, y: -R, width: 2 * R, height: 2 * R, fill: c.sheet, stroke: c.line, "stroke-width": 0.002 }));
        const cells = ev[L.key];
        for (let ix = 0; ix < L.nx; ix++) for (let iy = 0; iy < L.ny; iy++) {
          const e = cells[ix][iy], x = (ix - (L.nx - 1) / 2 - 0.5) * L.deta, y = -((iy - (L.ny - 1) / 2) + 0.5) * L.dphi;
          const on = e > Math.pow(10, LMIN);
          s.appendChild(el("rect", { x, y, width: L.deta, height: L.dphi, fill: on ? viridis((Math.log10(e) - LMIN) / (lmax - LMIN)) : "none",
            stroke: c.line, "stroke-width": L.key === "s1" ? 0.0006 : 0.0015 }));
        }
        s.appendChild(el("line", { x1: -R, y1: 0, x2: R, y2: 0, stroke: c.muted, "stroke-width": 0.0015, "stroke-dasharray": "0.006 0.006" }));
        s.appendChild(el("line", { x1: 0, y1: -R, x2: 0, y2: R, stroke: c.muted, "stroke-width": 0.0015, "stroke-dasharray": "0.006 0.006" }));
        const share = tot > 0 ? 100 * cells.flat().reduce((a, b) => a + b, 0) / tot : 0;
        const cap = document.createElement("figcaption");
        cap.innerHTML = `<b>${L.name}</b><span>${L.nx}×${L.ny} cells · ${share.toFixed(share < 1 ? 1 : 0)}%</span>`;
        fig.append(s, cap); return fig;
      }), depthChart(ev.depth, c));
      foot.innerHTML = `<span class="g4ev-bar" style="background:linear-gradient(90deg,${VIRIDIS.join(",")})"></span>
        <span>visible energy per cell, log scale: 0.1 MeV – ${Math.pow(10, lmax) >= 1000 ? (Math.pow(10, lmax) / 1000).toFixed(1) + " GeV" : Math.pow(10, lmax).toFixed(0) + " MeV"}</span>
        <span>Each map covers Δη × Δφ = 0.32 × 0.32 around the impact point (dashed lines); % = share of the visible energy in the windows.</span>`;
    }
    function depthChart(depth, c) {
      const fig = document.createElement("figure"); fig.className = "g4ev-map g4ev-depth";
      depth = depth.slice(0, 23);  // the EM layers end at 22.6 X0
      const W = 100, H = 100, n = depth.length, m = Math.max(1, ...depth);
      const s = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Longitudinal profile" });
      s.appendChild(el("rect", { x: 0, y: 0, width: W, height: H, fill: c.sheet, stroke: c.line, "stroke-width": 0.4 }));
      depth.forEach((v, i) => { const h = 88 * v / m; s.appendChild(el("rect", { x: 4 + i * (92 / n), y: 96 - h, width: 92 / n - 0.6, height: h, fill: c.accent })); });
      [0, 5, 10, 15, 20].forEach(t => s.appendChild(el("line", { x1: 4 + t * 92 / n, y1: 96, x2: 4 + t * 92 / n, y2: 99, stroke: c.ink, "stroke-width": 0.5 })));
      const cap = document.createElement("figcaption");
      cap.innerHTML = `<b>Depth profile</b><span>S1–S3 energy per X₀, 0–23 X₀</span>`;
      fig.append(s, cap); return fig;
    }
    draw();
    return { draw };
  }

  // ---------------------------------------------------------------- 3D view of one event, drawn by Geant4
  // The JSON files come from the VRML 2.0 driver of Geant4 (g4calo/macros/vrml.mac, g4calo/vis/wrl2json.py):
  // the bounding boxes of the placed volumes and the trajectories, coloured by particle. three.js is loaded only
  // when the viewer asks for the 3D view (import map "three" in the page).
  const TRACK_COL = { "e-": "#d62728", "e+": "#1f77b4", gamma: "#2ca02c", other: "#8c8c8c" };
  const TRACK_LAB = { "e-": "e⁻", "e+": "e⁺", gamma: "γ", other: "other" };
  const CUT = 150;  // mm: volumes drawn within ±15 cm of the beam axis

  function view3d(container, samples) {
    container.classList.add("g4v3");
    const ctl = document.createElement("div"); ctl.className = "g4ev-ctl";
    const sel = document.createElement("select"); sel.setAttribute("aria-label", "Simulated event");
    sel.innerHTML = samples.map((s, i) => `<option value="${i}">${s.label}</option>`).join("");
    const views = document.createElement("span"); views.className = "seg";
    views.innerHTML = '<button data-v="side" class="on">Side</button><button data-v="persp">3D</button><button data-v="beam">Along the beam</button>';
    const show = document.createElement("span"); show.className = "seg";
    show.innerHTML = Object.keys(TRACK_COL).map(k => `<button data-k="${k}"${k === "other" ? "" : ' class="on"'}><span class="g4v3-dot" style="background:${TRACK_COL[k]}"></span>${TRACK_LAB[k]}</button>`).join("");
    ctl.append(sel, views, show);
    const stage = document.createElement("div"); stage.className = "g4v3-stage";
    const start = document.createElement("button"); start.className = "g4v3-start";
    start.textContent = "Load the 3D view (three.js, ≈0.2 MB per event)";
    stage.append(start);
    const info = document.createElement("div"); info.className = "g4ev-foot";
    info.textContent = "Drag to rotate, scroll or pinch to zoom, right-drag to pan. Volumes are cut to ±15 cm around the beam axis.";
    container.replaceChildren(ctl, stage, info);

    let THREE, Orbit, renderer, scene, camera, controls, groups = {}, current = null, built = false;
    const cache = {};
    const load = i => cache[i] || (cache[i] = fetch(samples[i].file).then(r => r.json()));

    async function init() {
      start.disabled = true; start.textContent = "Loading…";
      try {
        THREE = await import("three");
        ({ OrbitControls: Orbit } = await import("three/addons/controls/OrbitControls.js"));
      } catch (e) { start.textContent = "Could not load three.js (" + e.message + ")"; return; }
      renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
      renderer.setPixelRatio(Math.min(2, devicePixelRatio || 1));
      stage.replaceChildren(renderer.domElement);
      scene = new THREE.Scene();
      camera = new THREE.PerspectiveCamera(35, 1, 1, 20000);
      controls = new Orbit(camera, renderer.domElement);
      controls.target.set(0, 0, 330);
      controls.enableDamping = true;
      new ResizeObserver(resize).observe(stage); resize();
      setView(views.querySelector(".on").dataset.v);
      built = true;
      await draw();
      (function loop() { controls.update(); renderer.render(scene, camera); requestAnimationFrame(loop); })();
    }
    function resize() {
      if (!renderer) return;
      const w = stage.clientWidth, h = stage.clientHeight;
      renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix();
    }
    function setView(v) {
      const t = controls.target;
      if (v === "side") camera.position.set(-1500, 120, t.z);  // beam from left to right
      if (v === "persp") camera.position.set(-1000, 650, -450);
      if (v === "beam") camera.position.set(-0.01, 0, -1100);
      camera.up.set(0, 1, 0); camera.lookAt(t); controls.update();
    }
    async function draw() {
      const d = await load(+sel.value), c = pal(), q = d.units_mm;
      if (current) { scene.remove(current); current.traverse(o => { if (o.geometry) o.geometry.dispose(); }); }
      current = new THREE.Group(); groups = {};
      // Volumes: faint absorber plates, and one outline per readout layer.
      const box = (x0, x1, y0, y1, z0, z1) => { const g = new THREE.BoxGeometry(x1 - x0, y1 - y0, z1 - z0);
        g.translate((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2); return g; };
      const clamp = v => Math.max(-CUT, Math.min(CUT, v));
      const plate = new THREE.MeshBasicMaterial({ color: c.muted, transparent: true, opacity: 0.07, depthWrite: false, side: THREE.DoubleSide });
      const edge = new THREE.LineBasicMaterial({ color: c.ink, transparent: true, opacity: 0.55 });
      const env = {};
      d.solids.forEach(([name, x0, x1, y0, y1, z0, z1]) => {
        const layer = { Pb: null, LAr_1: "S1", LAr_2: "S2", LAr_3: "S3", Fe: "Had", Scint: "Had", PS: "PS", Upstream: "Upstream" }[name];
        if (name === "Pb" || name === "Fe" || name === "Upstream") current.add(new THREE.Mesh(box(clamp(x0), clamp(x1), clamp(y0), clamp(y1), z0, z1), plate));
        if (layer) { const e = env[layer] || (env[layer] = [z0, z1]); e[0] = Math.min(e[0], z0); e[1] = Math.max(e[1], z1); }
      });
      // A period is Pb then LAr, so each layer starts where the previous one ends.
      if (env.S1 && env.PS) env.S1[0] = Math.min(env.S1[0], ...d.solids.filter(s => s[0] === "Pb").map(s => s[5]));
      if (env.S2 && env.S1) env.S2[0] = env.S1[1];
      if (env.S3 && env.S2) env.S3[0] = env.S2[1];
      Object.values(env).forEach(([z0, z1]) => current.add(new THREE.LineSegments(new THREE.EdgesGeometry(box(-CUT, CUT, -CUT, CUT, z0, z1)), edge)));
      // Trajectories
      Object.entries(d.tracks).forEach(([k, list]) => {
        const pos = [];
        list.forEach(t => { for (let i = 0; i + 5 < t.length; i += 3) pos.push(t[i] * q, t[i + 1] * q, t[i + 2] * q, t[i + 3] * q, t[i + 4] * q, t[i + 5] * q); });
        const g = new THREE.BufferGeometry(); g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
        const m = new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: TRACK_COL[k] || TRACK_COL.other, transparent: true, opacity: k === "gamma" ? 0.3 : 1 }));
        m.visible = show.querySelector(`[data-k="${k}"]`).classList.contains("on");
        groups[k] = m; current.add(m);
      });
      scene.add(current);
      info.textContent = `${d.label}: ${d.n_tracks_kept} of ${d.n_tracks_in} trajectories shown (longer than 3 mm, simplified to 0.5 mm). ` +
        "Neutrons, protons and other particles (grey) are hidden by default. Drag to rotate, scroll or pinch to zoom, right-drag to pan. " +
        "Volumes are cut to ±15 cm around the beam axis.";
    }
    start.onclick = init;
    sel.onchange = () => { if (built) draw(); };
    views.querySelectorAll("button").forEach(b => b.onclick = () => {
      views.querySelectorAll("button").forEach(x => x.classList.remove("on")); b.classList.add("on"); if (built) setView(b.dataset.v); });
    show.querySelectorAll("button").forEach(b => b.onclick = () => { b.classList.toggle("on"); if (groups[b.dataset.k]) groups[b.dataset.k].visible = b.classList.contains("on"); });
    return { draw() { if (built) draw(); } };
  }

  window.G4View = { detector, events, view3d };
})();
