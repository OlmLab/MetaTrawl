const DATA = "/data/";
const state = {
  catalog: null, entry: null, manifest: null, samples: null, stats: null,
  clusters: null, dendrogram: null, network: null, distributions: null,
  ani: null, genomeSort: "samples", selectedSample: null, dendrogramScale: 1,
};

const $ = (selector, root = document) => root.querySelector(selector);
const $$ = (selector, root = document) => [...root.querySelectorAll(selector)];
const fmt = new Intl.NumberFormat("en-US", { maximumFractionDigits: 2 });
const esc = value => String(value ?? "").replace(/[&<>"']/g, char => ({
  "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
}[char]));

async function getJSON(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`${response.status} while loading ${path}`);
  return response.json();
}

async function boot() {
  try {
    state.catalog = await getJSON(`${DATA}catalog.json`);
    renderCatalog();
    bindUI();
    const requested = decodeURIComponent(location.hash.replace(/^#\/?/, ""));
    const entry = state.catalog.genomes.find(item => item.genome === requested || item.path === requested)
      || state.catalog.genomes[0];
    if (!entry) throw new Error("The catalog contains no completed genome views.");
    await selectGenome(entry);
    $("#loading").hidden = true;
    $("#app").hidden = false;
  } catch (error) {
    $("#loading").hidden = true;
    $("#error").hidden = false;
    $("#error pre").textContent = error.stack || error.message;
  }
}

function bindUI() {
  $("#genome-search").addEventListener("input", renderCatalog);
  $("#sort-genomes").addEventListener("click", () => {
    state.genomeSort = state.genomeSort === "samples" ? "name" : "samples";
    $("#sort-genomes").textContent = state.genomeSort === "samples" ? "Samples ↓" : "Name ↑";
    renderCatalog();
  });
  $$(".tabs button").forEach(button => button.addEventListener("click", () => openTab(button.dataset.tab)));
  $$("[data-go]").forEach(button => button.addEventListener("click", () => openTab(button.dataset.go)));
  $("#distribution-select").addEventListener("change", renderDistribution);
  $("#sample-search").addEventListener("input", renderSampleTable);
  $("#rail-toggle").addEventListener("click", () => document.body.classList.toggle("rail-open"));
  $("#copy-link").addEventListener("click", async event => {
    await navigator.clipboard.writeText(location.href);
    event.currentTarget.textContent = "Copied";
    setTimeout(() => event.currentTarget.textContent = "Copy link", 1200);
  });
  $$("[data-zoom]").forEach(button => button.addEventListener("click", () => {
    const action = button.dataset.zoom;
    state.dendrogramScale = action === "reset" ? 1 : Math.max(.55, Math.min(2.5, state.dendrogramScale + (action === "in" ? .15 : -.15)));
    const svg = $("#dendrogram-wrap svg");
    if (svg) svg.style.transform = `scale(${state.dendrogramScale})`;
  }));
  $("#ani-filter").addEventListener("input", event => {
    $("#ani-output").textContent = Number(event.target.value).toFixed(2);
    renderNetwork();
  });
  $("#network-reset").addEventListener("click", () => {
    state.selectedSample = null;
    renderNetwork();
    renderInspector(null);
  });
  addEventListener("resize", debounce(() => {
    if ($("#tab-heatmap").classList.contains("active")) renderHeatmap();
    if ($("#tab-network").classList.contains("active")) renderNetwork();
  }, 120));
}

function renderCatalog() {
  if (!state.catalog) return;
  const query = ($("#genome-search")?.value || "").trim().toLowerCase();
  const entries = state.catalog.genomes
    .filter(item => String(item.genome).toLowerCase().includes(query))
    .sort((a, b) => state.genomeSort === "samples"
      ? (b.sample_count || 0) - (a.sample_count || 0) || String(a.genome).localeCompare(String(b.genome))
      : String(a.genome).localeCompare(String(b.genome)));
  $("#genome-count").textContent = `${fmt.format(entries.length)} genomes`;
  $("#catalog-date").textContent = `Updated ${new Date(state.catalog.generated_at).toLocaleString()}`;
  $("#genome-list").innerHTML = entries.map(entry => `
    <button class="genome-item ${state.entry?.path === entry.path ? "active" : ""}" data-path="${esc(entry.path)}">
      <b title="${esc(entry.genome)}">${esc(entry.genome)}</b><span>${fmt.format(entry.sample_count || 0)}</span>
      <em>${fmt.format(entry.strain_cluster_count || 0)} strains · ${fmt.format(entry.neighbor_edge_count || 0)} edges</em>
    </button>`).join("") || `<div class="empty">No matching genomes</div>`;
  $$(".genome-item").forEach(button => button.addEventListener("click", () => {
    const entry = state.catalog.genomes.find(item => item.path === button.dataset.path);
    selectGenome(entry);
  }));
}

async function selectGenome(entry) {
  if (!entry || entry.path === state.entry?.path) return;
  state.entry = entry;
  state.ani = null;
  state.selectedSample = null;
  renderCatalog();
  document.body.classList.remove("rail-open");
  $("#workspace").classList.add("busy");
  const base = `${DATA}${encodeURIComponent(entry.path)}/`;
  const manifest = await getJSON(base + "manifest.json");
  const files = manifest.files || {};
  const pathFor = (key, fallback) => base + (files[key]?.path || fallback);
  const [samples, stats, clusters, dendrogram, network, distributions] = await Promise.all([
    getJSON(pathFor("samples", "samples.json")),
    getJSON(pathFor("sample_stats", "sample_stats.json")),
    getJSON(pathFor("clusters", "clusters.json")),
    getJSON(pathFor("dendrogram", "dendrogram.json")),
    getJSON(pathFor("neighbor_network", "neighbor_network.json")),
    getJSON(pathFor("distributions", "distributions.json")),
  ]);
  Object.assign(state, { manifest, samples, stats, clusters, dendrogram, network, distributions, base });
  location.hash = `/${encodeURIComponent(entry.genome)}`;
  $("#genome-title").textContent = entry.genome;
  $("#download-stats").href = pathFor("sample_stats_parquet", "sample_stats.parquet");
  $("#clustermap-preview").src = pathFor("clustermap_preview", "clustermap.png");
  renderMetrics();
  renderDistributionOptions();
  renderDistribution();
  renderClusters();
  renderDendrogram();
  renderNetwork();
  renderSampleTable();
  $("#workspace").classList.remove("busy");
  if ($("#tab-heatmap").classList.contains("active")) await renderHeatmap();
}

function renderMetrics() {
  const m = state.manifest;
  const values = [
    ["Samples", m.sample_count],
    ["Strain clusters", m.strain_cluster_count],
    ["Clonal clusters", m.clonal_cluster_count],
    ["Neighbor edges", m.neighbor_edge_count],
    ["Minimum overlap", m.options?.min_comp_len],
  ];
  $("#metrics").innerHTML = values.map(([label, value]) =>
    `<div class="metric"><span>${label}</span><strong>${fmt.format(value ?? 0)}</strong></div>`).join("");
}

function renderDistributionOptions() {
  const labels = {
    genome_ani: "ANI (95–100%)", genome_ani_full: "ANI (full range)",
    total_positions: "Compared positions", coverage: "Coverage", breadth: "Breadth",
    ber: "BER", ref_ani: "Reference ANI", sylph_abundance: "Sylph abundance"
  };
  $("#distribution-select").innerHTML = Object.keys(state.distributions.histograms)
    .filter(key => state.distributions.histograms[key].counts?.length)
    .map(key => `<option value="${key}">${labels[key] || key}</option>`).join("");
}

function renderDistribution() {
  const key = $("#distribution-select").value;
  const hist = state.distributions.histograms[key];
  const target = $("#distribution-chart");
  if (!hist?.counts?.length) { target.innerHTML = `<div class="empty">No finite values</div>`; return; }
  const width = 760, height = 245, pad = { l: 36, r: 8, t: 15, b: 28 };
  const max = Math.max(...hist.counts, 1), barW = (width - pad.l - pad.r) / hist.counts.length;
  const bars = hist.counts.map((count, i) => {
    const h = (height - pad.t - pad.b) * count / max;
    return `<rect class="bar" x="${pad.l + i * barW}" y="${height - pad.b - h}" width="${Math.max(.8, barW + .15)}" height="${h}"><title>${fmt.format(hist.bin_edges[i])}–${fmt.format(hist.bin_edges[i + 1])}: ${fmt.format(count)}</title></rect>`;
  }).join("");
  const start = hist.bin_edges[0], end = hist.bin_edges.at(-1);
  target.innerHTML = `<svg viewBox="0 0 ${width} ${height}" preserveAspectRatio="none">${bars}
    <line x1="${pad.l}" y1="${height-pad.b}" x2="${width-pad.r}" y2="${height-pad.b}" stroke="rgba(17,17,17,.2)"/>
    <text class="axis-label" x="${pad.l}" y="${height-7}">${fmt.format(start)}</text>
    <text class="axis-label" x="${width-pad.r}" y="${height-7}" text-anchor="end">${fmt.format(end)}</text>
    <text class="axis-label" x="${pad.l}" y="10">${fmt.format(hist.value_count)} values</text></svg>`;
}

function renderClusters() {
  const groups = state.clusters.clusters;
  const block = (label, values, klass) => {
    const total = values.reduce((sum, item) => sum + item.sample_count, 0) || 1;
    const sorted = [...values].sort((a,b) => b.sample_count - a.sample_count);
    return `<div class="cluster-row"><header><b>${label}</b><span>${values.length} groups</span></header>
      <div class="cluster-strip ${klass}">${sorted.map(item => `<i style="width:${100*item.sample_count/total}%" title="Cluster ${item.cluster_id}: ${item.sample_count}"></i>`).join("")}</div>
      <p class="cluster-list">${sorted.slice(0, 5).map(item => `#${item.cluster_id} · ${fmt.format(item.sample_count)}`).join(" &nbsp; / &nbsp; ")}</p></div>`;
  };
  $("#cluster-summary").innerHTML = block("Strain clusters", groups.strain, "strain") + block("Clonal clusters", groups.clonal, "clonal");
}

function openTab(name) {
  $$(".tabs button").forEach(button => button.classList.toggle("active", button.dataset.tab === name));
  $$(".tab-panel").forEach(panel => panel.classList.toggle("active", panel.id === `tab-${name}`));
  if (name === "heatmap") renderHeatmap();
  if (name === "network") renderNetwork();
}

async function loadANI() {
  if (state.ani) return state.ani;
  if (!("DecompressionStream" in window)) throw new Error("This browser cannot decode gzip matrices. Use a current Chrome, Firefox, or Safari release.");
  const path = state.manifest.files?.similarity_matrix?.path || "similarity_ani.condensed.f32.gz";
  const response = await fetch(state.base + path);
  if (!response.ok) throw new Error(`Unable to load ANI matrix (${response.status})`);
  const decompressed = response.body.pipeThrough(new DecompressionStream("gzip"));
  const buffer = await new Response(decompressed).arrayBuffer();
  state.ani = new Float32Array(buffer);
  return state.ani;
}

function condensedIndex(n, i, j) {
  if (i === j) return -1;
  if (i > j) [i, j] = [j, i];
  return n * i - i * (i + 1) / 2 + j - i - 1;
}

async function renderHeatmap() {
  if (!state.samples) return;
  const canvas = $("#heatmap-canvas"), wrap = $("#heatmap-wrap"), tip = $("#heatmap-tip");
  sizeCanvas(canvas, wrap);
  const ctx = canvas.getContext("2d"), n = state.samples.sample_count;
  ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, canvas.width, canvas.height);
  try {
    const ani = await loadANI();
    const order = state.samples.leaf_order;
    const side = Math.min(canvas.width, canvas.height), x0 = (canvas.width-side)/2, y0 = (canvas.height-side)/2;
    const imageSide = Math.min(side, Math.max(2, n));
    const off = document.createElement("canvas"); off.width = imageSide; off.height = imageSide;
    const image = off.getContext("2d").createImageData(imageSide, imageSide);
    for (let y=0; y<imageSide; y++) for (let x=0; x<imageSide; x++) {
      const i = order[Math.min(n-1, Math.floor(y*n/imageSide))];
      const j = order[Math.min(n-1, Math.floor(x*n/imageSide))];
      const value = i === j ? 100 : ani[condensedIndex(n,i,j)];
      const color = aniColor(value);
      const k = (y*imageSide+x)*4;
      image.data[k]=color[0]; image.data[k+1]=color[1]; image.data[k+2]=color[2]; image.data[k+3]=255;
    }
    off.getContext("2d").putImageData(image,0,0);
    ctx.imageSmoothingEnabled = false; ctx.drawImage(off,x0,y0,side,side);
    canvas.onmousemove = event => {
      const p = canvasPoint(canvas,event), x=(p.x-x0)/side, y=(p.y-y0)/side;
      if (x<0||y<0||x>=1||y>=1) { tip.style.display="none"; return; }
      const ri=Math.min(n-1,Math.floor(y*n)), rj=Math.min(n-1,Math.floor(x*n));
      const i=order[ri], j=order[rj], value=i===j?100:ani[condensedIndex(n,i,j)];
      showTip(tip,event,`<b>${esc(state.samples.samples[i].sample_id)}</b><br>${esc(state.samples.samples[j].sample_id)}<br>ANI <strong>${Number(value).toFixed(4)}%</strong>`);
    };
    canvas.onmouseleave = () => tip.style.display="none";
  } catch (error) {
    ctx.fillStyle="#676762"; ctx.font="14px DINish"; ctx.textAlign="center";
    ctx.fillText(error.message,canvas.width/2,canvas.height/2);
  }
}

function aniColor(value) {
  const t = Math.max(0, Math.min(1, (Number(value)-95)/5));
  const stops = [[250,249,242],[213,160,33],[240,74,35]];
  const q=t*2, a=stops[Math.min(1,Math.floor(q))], b=stops[Math.min(2,Math.ceil(q))], f=q-Math.floor(q);
  return a.map((v,i)=>Math.round(v+(b[i]-v)*f));
}

function renderDendrogram() {
  const data=state.dendrogram, n=data.samples.length, row=Math.max(2.5, Math.min(18, 6400/n));
  const height=Math.max(580,n*row+50), width=1100, left=215, right=35;
  const sampleIndex=new Map(data.samples.map((sample,index)=>[sample,index])), nodes=new Map();
  data.ordered_samples.forEach((sample,rank)=>{
    const id=sampleIndex.get(sample), y=25+rank*row;
    nodes.set(id,{x:width-right,y});
  });
  const maxDist=Math.max(...data.linkage.rows.map(row=>row[2]),.001);
  let paths="";
  data.linkage.rows.forEach((merge,index)=>{
    const [l,r,d]=merge, a=nodes.get(l),b=nodes.get(r),x=width-right-(width-left-right)*d/maxDist;
    paths+=`<path d="M${a.x},${a.y}H${x}V${b.y}H${b.x}" fill="none" stroke="#1c4e80" stroke-width="1"/>`;
    nodes.set(n+index,{x,y:(a.y+b.y)/2});
  });
  const labels=n<=500?data.ordered_samples.map((sample,rank)=>`<text x="${width-right+7}" y="${29+rank*row}" font-family="Azeret Mono" font-size="${Math.max(5,Math.min(10,row*.7))}" fill="#676762">${esc(sample)}</text>`).join(""):"";
  const thresholds=[["clonal",state.manifest.options.clonal_cluster_threshold,"#f04a23"],["strain",state.manifest.options.strain_cluster_threshold,"#d5a021"]]
    .map(([label,ani,color])=>{const d=1-ani/100,x=width-right-(width-left-right)*d/maxDist;return `<line x1="${x}" y1="10" x2="${x}" y2="${height-15}" stroke="${color}" stroke-dasharray="4 4"/><text x="${x+4}" y="17" font-size="9" fill="${color}">${label} ${ani}%</text>`}).join("");
  $("#dendrogram-wrap").innerHTML=`<svg viewBox="0 0 ${width+(n<=500?260:0)} ${height}" style="height:${height}px">${paths}${thresholds}${labels}</svg>`;
}

function renderNetwork() {
  if (!state.network) return;
  const canvas=$("#network-canvas"),wrap=$("#network-wrap"),tip=$("#network-tip");
  sizeCanvas(canvas,wrap); const ctx=canvas.getContext("2d"), nodes=state.network.nodes, threshold=Number($("#ani-filter").value);
  const edges=state.network.edges.filter(edge=>edge.ani>=threshold), clusters=[...new Set(nodes.map(node=>node.strain_cluster))].sort();
  const clusterIndex=new Map(clusters.map((id,i)=>[id,i])), center={x:canvas.width/2,y:canvas.height/2};
  const radii=Math.min(canvas.width,canvas.height)*.34, positions=new Map();
  const grouped=new Map(); nodes.forEach(node=>{const a=grouped.get(node.strain_cluster)||[];a.push(node);grouped.set(node.strain_cluster,a)});
  grouped.forEach((group,cluster)=>{
    const ci=clusterIndex.get(cluster), angle=2*Math.PI*ci/clusters.length-Math.PI/2;
    const cx=center.x+Math.cos(angle)*radii*.55,cy=center.y+Math.sin(angle)*radii*.55;
    group.forEach((node,i)=>{const a=2*Math.PI*i/group.length+angle,r=Math.min(90,8+Math.sqrt(group.length)*7);positions.set(node.id,{x:cx+Math.cos(a)*r,y:cy+Math.sin(a)*r,node})});
  });
  ctx.clearRect(0,0,canvas.width,canvas.height); ctx.fillStyle="#ffffff";ctx.fillRect(0,0,canvas.width,canvas.height);
  ctx.lineWidth=1;
  edges.forEach(edge=>{const a=positions.get(edge.source),b=positions.get(edge.target);if(!a||!b)return;ctx.strokeStyle=`rgba(28,78,128,${.08+.45*Math.max(0,(edge.ani-threshold)/(100-threshold||1))})`;ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke()});
  positions.forEach(({x,y,node})=>{const selected=state.selectedSample===node.id;ctx.beginPath();ctx.arc(x,y,selected?6:Math.max(2.2,Math.min(4.5,2+Math.sqrt(node.coverage||0))),0,Math.PI*2);ctx.fillStyle=selected?"#f04a23":clusterColor(clusterIndex.get(node.strain_cluster));ctx.fill();});
  canvas.onmousemove=event=>{const p=canvasPoint(canvas,event);let nearest=null,dist=11;positions.forEach(pos=>{const d=Math.hypot(pos.x-p.x,pos.y-p.y);if(d<dist){nearest=pos;dist=d}});if(nearest)showTip(tip,event,`<b>${esc(nearest.node.id)}</b><br>Strain ${nearest.node.strain_cluster} · Clonal ${nearest.node.clonal_cluster}<br>Coverage ${value(nearest.node.coverage)}`);else tip.style.display="none"};
  canvas.onclick=event=>{const p=canvasPoint(canvas,event);let nearest=null,dist=13;positions.forEach(pos=>{const d=Math.hypot(pos.x-p.x,pos.y-p.y);if(d<dist){nearest=pos;dist=d}});if(nearest){state.selectedSample=nearest.node.id;renderInspector(nearest.node);renderNetwork()}};
  canvas.onmouseleave=()=>tip.style.display="none";
}

function clusterColor(index) {
  const palette=["#1c4e80","#f04a23","#d5a021","#264653","#6a704c","#8e4f3e","#55738c","#8c8055"];
  return palette[index%palette.length];
}

function renderInspector(node) {
  const target=$("#sample-inspector");
  if(!node){target.innerHTML="<span>Select a node</span><h3>Sample details</h3><p>Click any sample to inspect its cluster and profile statistics.</p>";return}
  target.innerHTML=`<span>Selected sample</span><h3>${esc(node.id)}</h3><dl>
    <dt>Strain cluster</dt><dd>${node.strain_cluster}</dd><dt>Clonal cluster</dt><dd>${node.clonal_cluster}</dd>
    <dt>Coverage</dt><dd>${value(node.coverage)}</dd><dt>Breadth</dt><dd>${value(node.breadth)}</dd>
    <dt>BER</dt><dd>${value(node.ber)}</dd><dt>Reference ANI</dt><dd>${value(node.ref_ani)}</dd>
    <dt>Sylph abundance</dt><dd>${value(node.sylph_abundance)}</dd></dl>`;
}

function renderSampleTable() {
  if(!state.stats)return;
  const columns=state.stats.columns, query=$("#sample-search").value.trim().toLowerCase();
  const preferred=["sample_id","coverage","breadth","ber","ref_ani","sylph_abundance","strain_cluster","clonal_cluster","null_fraction"];
  const names=preferred.filter(name=>columns[name]).concat(Object.keys(columns).filter(name=>!preferred.includes(name)).slice(0,5));
  const rows=Array.from({length:state.stats.row_count},(_,i)=>Object.fromEntries(names.map(name=>[name,columns[name][i]])))
    .filter(row=>String(row.sample_id).toLowerCase().includes(query));
  $("#sample-table").innerHTML=`<thead><tr>${names.map(name=>`<th>${esc(name.replaceAll("_"," "))}</th>`).join("")}</tr></thead>
    <tbody>${rows.slice(0,5000).map(row=>`<tr>${names.map(name=>`<td>${name==="sample_id"?esc(row[name]):value(row[name])}</td>`).join("")}</tr>`).join("")}</tbody>`;
}

function sizeCanvas(canvas,wrap) {
  const ratio=Math.min(devicePixelRatio||1,2),rect=wrap.getBoundingClientRect();
  canvas.width=Math.max(1,Math.floor(rect.width*ratio));canvas.height=Math.max(1,Math.floor((rect.height||650)*ratio));
}
function canvasPoint(canvas,event){const rect=canvas.getBoundingClientRect();return{x:(event.clientX-rect.left)*canvas.width/rect.width,y:(event.clientY-rect.top)*canvas.height/rect.height}}
function showTip(tip,event,html){tip.innerHTML=html;tip.style.display="block";const parent=tip.parentElement.getBoundingClientRect();tip.style.left=`${Math.min(parent.width-280,event.clientX-parent.left+14)}px`;tip.style.top=`${Math.max(8,event.clientY-parent.top-28)}px`}
function value(v){return v==null?"—":typeof v==="number"?fmt.format(v):esc(v)}
function debounce(fn,wait){let timer;return(...args)=>{clearTimeout(timer);timer=setTimeout(()=>fn(...args),wait)}}

boot();
