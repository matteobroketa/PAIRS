(() => {
  "use strict";

  const SUPPORTED_SCHEMA = 4;
  const DATA_ROOT = `data/v${SUPPORTED_SCHEMA}`;
  const PAGE_RENDER_SIZE = 30;
  const $ = selector => document.querySelector(selector);
  const $$ = selector => [...document.querySelectorAll(selector)];

  const els = {
    q: $("#query"),
    search: $("#searchBtn"),
    suggestions: $("#suggestions"),
    textMode: $("#textModeBtn"),
    sequenceMode: $("#sequenceModeBtn"),
    textPanel: $("#textSearchPanel"),
    sequencePanel: $("#sequenceSearchPanel"),
    heavySequence: $("#heavySequenceQuery"),
    lightSequence: $("#lightSequenceQuery"),
    sequenceSearch: $("#sequenceSearchBtn"),
    sequenceType: $("#sequenceType"),
    nearMatches: $("#nearMatches"),
    batchQuery: $("#batchSequenceQuery"),
    batchSearch: $("#batchSearchBtn"),
    batchClear: $("#batchClearBtn"),
    batchCsv: $("#batchCsvBtn"),
    batchProgress: $("#batchProgress"),
    batchResults: $("#batchResults"),
    main: $("#main"),
    targetName: $("#targetName"),
    targetMeta: $("#targetMeta"),
    results: $("#results"),
    summary: $("#summaryGrid"),
    loadMore: $("#loadMore"),
    sort: $("#sortSelect"),
    status: $("#statusStrip"),
    warning: $("#sourceWarning"),
    footer: $("#footerSnapshot"),
    modal: $("#sourceModal"),
    sourceList: $("#sourceList"),
    browseQuery: $("#browseQuery"),
    browseSort: $("#browseSort"),
    browseMeta: $("#browseMeta"),
    targetGrid: $("#targetGrid"),
    browseGrid: $("#browseGridBtn"),
    browseField: $("#browseFieldBtn"),
    bindingField: $("#bindingField"),
    bindingFieldSvg: $("#bindingFieldSvg"),
    fieldTooltip: $("#fieldTooltip"),
    heroSvg: $("#heroAmbientSvg"),
    browseMore: $("#browseMore"),
    filters: $("#filters"),
    browseFacets: $("#browseFacets"),
    copyUrl: $("#copyUrlBtn"),
    csv: $("#csvBtn"),
    fasta: $("#fastaBtn"),
  };

  const state = {
    manifest: null,
    targets: [],
    selected: null,
    mode: "target",
    rawResults: [],
    filtered: [],
    shown: PAGE_RENDER_SIZE,
    filter: "all",
    activeSuggestion: -1,
    suggestionItems: [],
    suggestionToken: 0,
    antibodySearchCache: new Map(),
    antibodyShardCache: new Map(),
    sequenceSearchCache: new Map(),
    cdrBucketCache: new Map(),
    targetPageCache: new Map(),
    loadedTargetPages: new Set(),
    loadedFunctionalPages: new Set(),
    loadedNegativePages: new Set(),
    browseShown: 60,
    browseView: "grid",
    fieldObserver: null,
    heroObserver: null,
    sequenceQueryLabel: "",
    batchRows: [],
  };

  const evidenceRank = {
    STRUCTURE: 6,
    MEASURED: 5,
    CURATED: 4,
    PUBLICATION: 3,
    LITERATURE_METADATA: 2,
    PATENT: 1,
    METADATA: 0,
  };

  const relationLabel = relationship =>
    ({
      binds: "binds",
      does_not_bind: "does not bind",
      neutralizes: "neutralizes",
      does_not_neutralize: "does not neutralize",
      protects: "protects",
      does_not_protect: "does not protect",
      targets: "curated target",
      mentioned_with: "literature mention",
    })[relationshipValue(relationship)] || relationshipValue(relationship).replaceAll("_", " ");

  const esc = value =>
    String(value ?? "").replace(
      /[&<>"']/g,
      character =>
        ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[character],
    );

  const residueGroup = residue => {
    const groups = {
      A: "hydrophobic",
      I: "hydrophobic",
      L: "hydrophobic",
      M: "hydrophobic",
      F: "hydrophobic",
      W: "hydrophobic",
      V: "hydrophobic",
      Y: "hydrophobic",
      K: "basic",
      R: "basic",
      H: "basic",
      D: "acidic",
      E: "acidic",
      S: "polar",
      T: "polar",
      N: "polar",
      Q: "polar",
      G: "glycine",
      P: "proline",
      C: "cysteine",
    };
    return groups[residue] || "other";
  };

  // Display-only coloring: copied/exported sequences continue to use the
  // original plain value from data-copy-seq and never include markup.
  const coloredSequence = value =>
    [...String(value || "")]
      .map(
        residue => `<span class="residue residue-${residueGroup(residue)}">${esc(residue)}</span>`,
      )
      .join("");

  const stateGlyph = kind =>
    `<span class="state-glyph state-glyph-${kind}" aria-hidden="true">${kind === "error" ? "!" : kind === "empty" ? "∅" : kind === "loading" ? "⋯" : "✓"}</span>`;
  const loadingState = (title, detail = "") =>
    `<div class="empty state-panel state-loading" role="status" aria-live="polite">${stateGlyph("loading")}<strong>${esc(title)}</strong>${detail ? `<span>${esc(detail)}</span>` : ""}<span class="skeleton-line"></span><span class="skeleton-line short"></span></div>`;
  const emptyState = (title, detail = "") =>
    `<div class="empty state-panel">${stateGlyph("empty")}<strong>${esc(title)}</strong>${detail ? `<span>${esc(detail)}</span>` : ""}</div>`;
  const errorState = (title, detail = "") =>
    `<div class="empty state-panel state-error" role="alert">${stateGlyph("error")}<strong>${esc(title)}</strong>${detail ? `<span>${esc(detail)}</span>` : ""}</div>`;

  const norm = value =>
    String(value ?? "")
      .toLowerCase()
      .normalize("NFKD")
      .replace(/[^a-z0-9]+/g, " ")
      .trim();

  const compact = value => norm(value).replaceAll(" ", "");
  const fmt = value => new Intl.NumberFormat().format(value || 0);
  const prefersReducedMotion = () =>
    Boolean(window.matchMedia?.("(prefers-reduced-motion: reduce)")?.matches);
  const relationshipValue = relationship =>
    typeof relationship === "string" ? relationship : relationship?.relationship || "";
  const isCompleteSnapshot = () => {
    if (!state.manifest) return false;
    const sources = Object.values(state.manifest.sources || {});
    const okay = state.manifest.sources_ok ?? sources.filter(source => source.ok).length;
    const expected = state.manifest.sources_expected ?? sources.length;
    return expected > 0 && okay === expected;
  };
  const isNegative = relationship =>
    relationshipValue(relationship).includes("does_not") ||
    relationshipValue(relationship).includes("not_");
  const isLiteratureContext = relationship =>
    ["mentioned_with", "literature", "literature_mention"].includes(
      relationshipValue(relationship),
    );
  const isPrimaryPositive = relationship =>
    ["binds", "targets"].includes(relationshipValue(relationship));
  const isFunctionalPositive = relationship =>
    ["neutralizes", "protects"].includes(relationshipValue(relationship));

  function targetNames(values) {
    return (values || [])
      .map(value =>
        typeof value === "string" ? value : value?.name || value?.target_name || value?.label,
      )
      .filter(Boolean);
  }

  function targetEntries(values) {
    return (values || [])
      .map(value => {
        if (typeof value === "string") return { name: value, id: "" };
        return {
          id: value?.id || value?.target_id || "",
          name: value?.name || value?.target_name || value?.label || "",
        };
      })
      .filter(value => value.name);
  }

  function directTargetNames(antibody) {
    return targetNames(
      antibody.direct_targets || antibody.directTargets || antibody.targets_supported || [],
    );
  }

  function functionalActivityNames(antibody) {
    return targetNames(
      antibody.functional_activity ||
        antibody.functionalActivity ||
        antibody.functional_targets ||
        antibody.functionalTargets ||
        [],
    );
  }

  function negativeEvidenceNames(antibody) {
    return targetNames(antibody.negative_evidence || antibody.negativeEvidence || []);
  }

  function literatureMentionNames(antibody) {
    return targetNames(antibody.literature_mentions || antibody.literatureMentions || []);
  }

  function legacyTargetNames(antibody) {
    return targetNames(antibody.targets || []);
  }

  function structureEntries(values, defaultTier = "unknown") {
    return (values || [])
      .map(value => {
        if (typeof value === "string") return { id: value, tier: defaultTier };
        const id = value?.pdb_id || value?.pdb || value?.structure_id || value?.id || "";
        const rawTier = String(
          value?.tier || value?.structure_tier || value?.identity_tier || value?.status || "",
        ).toLowerCase();
        const identity = Number(value?.sequence_identity ?? value?.identity ?? NaN);
        const tier =
          rawTier.includes("exact") || rawTier === "100" || identity >= 100
            ? "exact"
            : rawTier.includes("homolog") ||
                rawTier.includes("99") ||
                rawTier.includes("95") ||
                identity > 0
              ? "homologous"
              : defaultTier;
        return { id, tier, identity: Number.isFinite(identity) ? identity : null };
      })
      .filter(value => value.id);
  }

  function structureCollections(antibody) {
    const tiered = Object.entries(antibody.structure_tiers || {}).flatMap(([tier, values]) => {
      const identity = Number.parseFloat(tier);
      const defaultTier =
        identity === 100 ? "exact" : Number.isFinite(identity) ? "homologous" : "unknown";
      return structureEntries(values, defaultTier).map(entry => ({
        ...entry,
        identity: Number.isFinite(identity) ? identity : entry.identity,
        identityLabel: tier,
      }));
    });
    // When tier data exists it is authoritative. The legacy `structures` and
    // `exact_structures` fields in older builds may contain homologous IDs.
    const hasTierData = tiered.length > 0;
    const exact = hasTierData
      ? []
      : structureEntries(
          antibody.exact_structures || antibody.structures_exact || antibody.exactStructures || [],
          "exact",
        );
    const homologous = hasTierData
      ? []
      : structureEntries(
          antibody.homologous_structures ||
            antibody.structures_homologous ||
            antibody.homologousStructures ||
            [],
          "homologous",
        );
    const records = structureEntries(
      antibody.structure_records || antibody.structureRecords || antibody.structure_evidence || [],
    );
    const legacy = hasTierData ? [] : structureEntries(antibody.structures || []);
    const all = [...exact, ...homologous, ...tiered, ...records, ...legacy];
    const byId = new Map();
    for (const entry of all) {
      const previous = byId.get(entry.id);
      if (!previous || (previous.tier === "unknown" && entry.tier !== "unknown")) {
        byId.set(entry.id, entry);
      }
    }
    return {
      exact: [...byId.values()].filter(entry => entry.tier === "exact"),
      homologous: [...byId.values()].filter(entry => entry.tier === "homologous"),
      unknown: [...byId.values()].filter(entry => entry.tier === "unknown"),
    };
  }

  const hasExactStructure = antibody => structureCollections(antibody).exact.length > 0;

  const jsonCache = new Map();
  async function getJSON(url) {
    if (jsonCache.has(url)) return jsonCache.get(url);
    const request = fetch(url, { cache: "no-cache" }).then(async response => {
      if (!response.ok) {
        if (response.status === 404 && url.startsWith(DATA_ROOT)) {
          throw new Error(
            "A versioned dataset file is unavailable. The site may have been updated while this tab was open; reload the page.",
          );
        }
        throw new Error(`${url}: HTTP ${response.status}`);
      }
      return response.json();
    });
    jsonCache.set(url, request);
    try {
      return await request;
    } catch (error) {
      jsonCache.delete(url);
      throw error;
    }
  }

  function levenshteinBounded(left, right, maximum = 2) {
    if (Math.abs(left.length - right.length) > maximum) return maximum + 1;
    const previous = Array.from({ length: right.length + 1 }, (_, index) => index);
    for (let row = 1; row <= left.length; row += 1) {
      const current = [row];
      let rowMinimum = current[0];
      for (let column = 1; column <= right.length; column += 1) {
        const cost = left[row - 1] === right[column - 1] ? 0 : 1;
        current[column] = Math.min(
          current[column - 1] + 1,
          previous[column] + 1,
          previous[column - 1] + cost,
        );
        rowMinimum = Math.min(rowMinimum, current[column]);
      }
      if (rowMinimum > maximum) return maximum + 1;
      previous.splice(0, previous.length, ...current);
    }
    return previous[right.length];
  }

  function textScore(value, query) {
    const normalizedValue = norm(value);
    const normalizedQuery = norm(query);
    if (!normalizedValue || !normalizedQuery) return -1;
    if (normalizedValue === normalizedQuery) return 10000;
    if (normalizedValue.startsWith(normalizedQuery)) return 8000;
    if (normalizedValue.includes(normalizedQuery)) return 5000;
    if (normalizedQuery.length >= 4 && normalizedValue.length >= 4) {
      const distance = levenshteinBounded(normalizedValue, normalizedQuery, 2);
      if (distance <= 2) return 3200 - distance * 450;
    }
    return -1;
  }

  function scoreTarget(target, query) {
    let score = textScore(target.name, query);
    for (const alias of target.aliases || []) score = Math.max(score, textScore(alias, query));
    return score >= 0 ? score + Math.min(target.result_count || target.count || 0, 999) : -1;
  }

  function matchingTargets(query, limit = 12) {
    return state.targets
      .map(target => [scoreTarget(target, query), target])
      .filter(([score]) => score >= 0)
      .sort((left, right) => right[0] - left[0])
      .slice(0, limit)
      .map(([, target]) => target);
  }

  function antibodyBucket(query) {
    const value = compact(query);
    return value.length >= 2 ? value.slice(0, 2) : `${value}_`;
  }

  async function searchAntibodyNames(query, limit = 10) {
    const normalizedQuery = norm(query);
    if (!normalizedQuery) return [];
    const bucket = antibodyBucket(query);
    let index = state.antibodySearchCache.get(bucket);
    if (!index) {
      try {
        index = await getJSON(`${DATA_ROOT}/antibody-search/${encodeURIComponent(bucket)}.json`);
      } catch {
        index = [];
      }
      state.antibodySearchCache.set(bucket, index);
    }
    return index
      .map(antibody => {
        let score = textScore(antibody.name, query);
        for (const alias of antibody.aliases || [])
          score = Math.max(score, textScore(alias, query));
        return [score, antibody];
      })
      .filter(([score]) => score >= 0)
      .sort((left, right) => right[0] - left[0])
      .slice(0, limit)
      .map(([, antibody]) => antibody);
  }

  function exactTargetMatch(query) {
    const normalized = norm(query);
    if (!normalized) return null;
    return (
      state.targets.find(target =>
        [target.id, target.name, ...(target.aliases || [])].some(
          value => norm(value) === normalized,
        ),
      ) || null
    );
  }

  async function exactSearchResolution(query) {
    const target = exactTargetMatch(query);
    if (target) return { kind: "target", target };
    const normalized = norm(query);
    const antibodies = await searchAntibodyNames(query, 100);
    const antibody = antibodies.find(item =>
      [item.id, item.name, ...(item.aliases || [])].some(value => norm(value) === normalized),
    );
    return antibody ? { kind: "antibody", antibody } : null;
  }

  function setActiveSuggestion(index) {
    const suggestions = $$(".suggestion[data-index]");
    if (!suggestions.length) {
      state.activeSuggestion = -1;
      return;
    }
    state.activeSuggestion = Math.max(0, Math.min(index, suggestions.length - 1));
    suggestions.forEach((suggestion, itemIndex) => {
      const active = itemIndex === state.activeSuggestion;
      suggestion.classList.toggle("active", active);
      suggestion.setAttribute("aria-selected", active ? "true" : "false");
      if (active) suggestion.scrollIntoView({ block: "nearest" });
    });
  }

  function renderNoSuggestion(query) {
    const issue = feedbackIssueUrl(query);
    const report = issue
      ? `<a class="secondary" href="${esc(issue)}" target="_blank" rel="noopener">Report missing query</a>`
      : `<button class="secondary" data-copy-missing="${esc(query)}">Copy missing-query report</button>`;
    return `<div class="suggestion no-hit"><div class="s-main"><strong>No local match</strong><span>Try a synonym, gene symbol, antigen domain, therapeutic name, or browse all targets.</span></div>${report}</div>`;
  }

  async function showSuggestions() {
    const query = els.q.value.trim();
    const token = ++state.suggestionToken;
    if (query.length < 2) {
      closeSuggestions();
      return;
    }

    const targets = matchingTargets(query, 8);
    let antibodies = [];
    if (query.length >= 3) {
      try {
        antibodies = await searchAntibodyNames(query, 6);
      } catch {
        antibodies = [];
      }
    }
    if (token !== state.suggestionToken || query !== els.q.value.trim()) return;

    state.suggestionItems = [];
    const html = [];
    for (const target of targets) {
      const index = state.suggestionItems.push({ kind: "target", target }) - 1;
      const aliases = (target.aliases || [])
        .filter(alias => norm(alias) !== norm(target.name))
        .slice(0, 6)
        .join(" · ");
      html.push(
        `<div class="suggestion" role="option" data-index="${index}"><div class="s-main"><strong>${esc(target.name)}</strong><span>${esc(aliases || target.sources.join(" · "))}</span></div><span class="count">${fmt(target.result_count)} antibodies</span></div>`,
      );
    }
    for (const antibody of antibodies) {
      const index = state.suggestionItems.push({ kind: "antibody", antibody }) - 1;
      const directTargets = directTargetNames(antibody);
      const associatedTargets = legacyTargetNames(antibody);
      const targetSummary = directTargets.length
        ? `Direct target · ${directTargets.slice(0, 4).join(" · ")}`
        : associatedTargets.length
          ? `Associated annotation · ${associatedTargets.slice(0, 4).join(" · ")}`
          : `Antibody · ${(antibody.sources || []).join(" · ")}`;
      html.push(
        `<div class="suggestion" role="option" data-index="${index}"><div class="s-main"><strong>${esc(antibody.name)}</strong><span>${esc(targetSummary)}</span></div><span class="count">${antibody.paired ? "VH + VL" : "sequence record"}</span></div>`,
      );
    }

    els.suggestions.innerHTML = html.join("") || renderNoSuggestion(query);
    $$("#suggestions .suggestion").forEach((item, index) =>
      item.style.setProperty("--suggestion-index", index),
    );
    els.suggestions.classList.add("open");
    state.activeSuggestion = -1;
  }

  function closeSuggestions() {
    els.suggestions.classList.remove("open");
    state.activeSuggestion = -1;
  }

  function activateSuggestion(index) {
    const item = state.suggestionItems[index];
    if (!item) return;
    closeSuggestions();
    if (item.kind === "target") selectTarget(item.target, true);
    else openStandaloneAntibody(item.antibody.id, true);
  }

  function setSearchMode(mode) {
    const sequence = mode === "sequence";
    els.textPanel.hidden = sequence;
    els.sequencePanel.hidden = !sequence;
    els.textMode.classList.toggle("active", !sequence);
    els.sequenceMode.classList.toggle("active", sequence);
    els.textMode.setAttribute("aria-selected", sequence ? "false" : "true");
    els.sequenceMode.setAttribute("aria-selected", sequence ? "true" : "false");
    const activePanel = sequence ? els.sequencePanel : els.textPanel;
    activePanel?.classList.remove("panel-enter");
    requestAnimationFrame(() => activePanel?.classList.add("panel-enter"));
    const activeTab = sequence ? els.sequenceMode : els.textMode;
    if (activeTab && els.textMode?.parentElement) {
      els.textMode.parentElement.style.setProperty("--tab-left", `${activeTab.offsetLeft}px`);
      els.textMode.parentElement.style.setProperty("--tab-width", `${activeTab.offsetWidth}px`);
    }
    closeSuggestions();
    if (sequence) els.heavySequence.focus();
    else els.q.focus();
  }

  function looksLikeSequence(value) {
    const candidate = value.replace(/\s|-/g, "").toUpperCase();
    return candidate.length >= 25 && /^[ACDEFGHIKLMNPQRSTVWYX*]+$/.test(candidate);
  }

  function normalizePastedSequence(value) {
    const withoutHeaders = value
      .split(/\r?\n/)
      .filter(line => !line.trim().startsWith(">"))
      .join("");
    const sequence = withoutHeaders.replace(/\s/g, "").toUpperCase();
    if (!sequence) throw new Error("Empty sequence.");
    const invalid = sequence.match(/[^ACDEFGHIKLMNPQRSTVWYX]/);
    if (invalid) throw new Error(`Invalid amino-acid character '${invalid[0]}'.`);
    return sequence;
  }

  async function sha256Hex(value) {
    if (!globalThis.crypto?.subtle) {
      throw new Error("Exact sequence lookup requires Web Crypto. Use HTTPS or localhost.");
    }
    const bytes = new TextEncoder().encode(value);
    const digestBuffer = await crypto.subtle.digest("SHA-256", bytes);
    return [...new Uint8Array(digestBuffer)]
      .map(byte => byte.toString(16).padStart(2, "0"))
      .join("");
  }

  async function lookupSequence(sequence) {
    const hash = await sha256Hex(sequence);
    const bucket = hash.slice(0, 2);
    let payload = state.sequenceSearchCache.get(bucket);
    if (!payload) {
      try {
        payload = await getJSON(`${DATA_ROOT}/sequence-search/${bucket}.json`);
      } catch {
        payload = {};
      }
      state.sequenceSearchCache.set(bucket, payload);
    }
    return payload[hash] || [];
  }

  function cdrThreshold(length) {
    return length <= 8 ? 1 : 2;
  }
  async function lookupCdr(sequence, field, near = true) {
    const lengths = [sequence.length - 1, sequence.length, sequence.length + 1].filter(
      length => length > 0,
    );
    const buckets = await Promise.all(
      lengths.map(length => {
        const key = `${field}:${length}`;
        if (!state.cdrBucketCache.has(key))
          state.cdrBucketCache.set(
            key,
            getJSON(`${DATA_ROOT}/sequence/${field}/${String(length).padStart(2, "0")}.json`).catch(
              () => [],
            ),
          );
        return state.cdrBucketCache.get(key);
      }),
    );
    const maximum = cdrThreshold(sequence.length);
    return buckets
      .flat()
      .map(candidate => ({
        ...candidate,
        distance: levenshteinBounded(sequence, candidate.sequence, near ? maximum : 0),
      }))
      .filter(candidate => candidate.distance <= (near ? maximum : 0))
      .sort((a, b) => a.distance - b.distance || a.sequence.localeCompare(b.sequence));
  }

  function combineSequenceMatches(first, second = null) {
    const byAntibody = new Map();
    if (!second) {
      for (const match of first) {
        const item = byAntibody.get(match.id) || { ...match, match_fields: [] };
        if (!item.match_fields.includes(match.field)) item.match_fields.push(match.field);
        byAntibody.set(match.id, item);
      }
      return [...byAntibody.values()];
    }

    const firstById = new Map();
    for (const match of first.filter(item => ["heavy", "cdrh3"].includes(item.field))) {
      const fields = firstById.get(match.id) || [];
      fields.push(match.field);
      firstById.set(match.id, fields);
    }
    for (const match of second.filter(item => ["light", "cdrl3"].includes(item.field))) {
      if (!firstById.has(match.id)) continue;
      const base = byAntibody.get(match.id) || { ...match, match_fields: [] };
      base.match_fields = [...new Set([...firstById.get(match.id), match.field])];
      byAntibody.set(match.id, base);
    }
    return [...byAntibody.values()];
  }

  function antibodyShardFromId(antibodyId) {
    return antibodyId.startsWith("ab_") ? antibodyId.slice(3, 5) : "";
  }

  async function fetchAntibody(antibodyId, shard = antibodyShardFromId(antibodyId)) {
    if (!shard) return null;
    let payload = state.antibodyShardCache.get(shard);
    if (!payload) {
      payload = await getJSON(`${DATA_ROOT}/antibodies/${shard}.json`);
      state.antibodyShardCache.set(shard, payload);
    }
    return payload[antibodyId] || null;
  }

  async function fetchFullRecords(antibodyIds) {
    const ids = [...new Set(antibodyIds)];
    const byShard = new Map();
    for (const antibodyId of ids) {
      const shard = antibodyShardFromId(antibodyId);
      if (!byShard.has(shard)) byShard.set(shard, []);
      byShard.get(shard).push(antibodyId);
    }

    await Promise.all(
      [...byShard.keys()].map(async shard => {
        if (state.antibodyShardCache.has(shard)) return;
        const payload = await getJSON(`${DATA_ROOT}/antibodies/${shard}.json`);
        state.antibodyShardCache.set(shard, payload);
      }),
    );

    const records = new Map();
    for (const [shard, shardIds] of byShard) {
      const payload = state.antibodyShardCache.get(shard) || {};
      for (const antibodyId of shardIds) {
        if (payload[antibodyId]) records.set(antibodyId, payload[antibodyId]);
      }
    }
    return records;
  }

  function summaryFromFull(antibody) {
    return {
      id: antibody.id,
      name: antibody.name || antibody.id,
      aliases: (antibody.aliases || []).slice(0, 6),
      organism: antibody.organism || "",
      format: antibody.format || "",
      has_heavy: Boolean(antibody.heavy),
      has_light: Boolean(antibody.light),
      heavy_length: (antibody.heavy || "").length,
      light_length: (antibody.light || "").length,
      structures: (antibody.structures || []).slice(0, 12),
      exact_structures: (antibody.exact_structures || antibody.structures_exact || []).slice(0, 12),
      homologous_structures: (
        antibody.homologous_structures ||
        antibody.structures_homologous ||
        []
      ).slice(0, 12),
      structure_tiers: antibody.structure_tiers || {},
      structure_records: (
        antibody.structure_records ||
        antibody.structureRecords ||
        antibody.structure_evidence ||
        []
      ).slice(0, 12),
      therapeutic_status: antibody.therapeutic_status || "",
      shard: antibodyShardFromId(antibody.id),
    };
  }

  async function runSequenceSearch(firstRaw, secondRaw = "") {
    els.sequenceSearch.disabled = true;
    const original = els.sequenceSearch.textContent;
    els.sequenceSearch.textContent = "Searching…";
    els.results.innerHTML = loadingState(
      "Searching local sequence index…",
      "Your sequence stays in this browser.",
    ).replace("state-loading", "state-loading sequence-scan");
    try {
      const first = normalizePastedSequence(firstRaw);
      const second = secondRaw ? normalizePastedSequence(secondRaw) : "";
      const type = els.sequenceType.value;
      if (["cdrh3", "cdrl3"].includes(type)) {
        const hits = await lookupCdr(first, type, els.nearMatches.checked);
        const records = await fetchFullRecords([...new Set(hits.flatMap(hit => hit.antibody_ids))]);
        state.mode = "sequence";
        state.selected = null;
        state.rawResults = hits.flatMap(hit =>
          hit.antibody_ids.map(id => ({
            antibody: summaryFromFull(records.get(id) || { id }),
            relationships: [],
            evidence: [],
            sources: records.get(id)?.sources || [],
            interactions: [],
            match_fields: [
              hit.distance
                ? `${type.toUpperCase()} near match · edit distance ${hit.distance}`
                : `Exact ${type.toUpperCase()} match`,
            ],
          })),
        );
        state.sequenceQueryLabel = `${type.toUpperCase()} local similarity search`;
        state.filter = "all";
        state.shown = PAGE_RENDER_SIZE;
        resetFilterButtons();
        updateFilterAvailability();
        els.targetName.textContent = `${type.toUpperCase()} sequence matches`;
        els.targetMeta.textContent = `${fmt(state.rawResults.length)} public records · similarity is retrieval evidence, not biological identity`;
        els.main.classList.add("active");
        apply();
        if (!state.rawResults.length) renderSequenceEmpty(first.length, false);
        history.pushState({}, "", location.pathname);
        scrollToResults();
        return;
      }
      const [firstMatches, secondMatches] = await Promise.all([
        lookupSequence(first),
        second ? lookupSequence(second) : Promise.resolve(null),
      ]);
      const combined = combineSequenceMatches(firstMatches, secondMatches);
      const records = await fetchFullRecords(combined.map(match => match.id));
      const matchesById = new Map(combined.map(match => [match.id, match]));
      state.mode = "sequence";
      state.selected = null;
      state.rawResults = [...records.values()].map(antibody => ({
        antibody: summaryFromFull(antibody),
        relationships: [],
        evidence: [],
        sources: antibody.sources || [],
        interactions: [],
        match_fields: matchesById.get(antibody.id)?.match_fields || [],
      }));
      state.filter = "all";
      state.shown = PAGE_RENDER_SIZE;
      state.sequenceQueryLabel = second
        ? `Exact paired sequence · ${first.length} aa + ${second.length} aa`
        : `Exact sequence · ${first.length} aa`;
      resetFilterButtons();
      updateFilterAvailability();
      els.targetName.textContent = "Exact sequence matches";
      els.targetMeta.textContent = `${fmt(state.rawResults.length)} public antibody records · ${state.sequenceQueryLabel} · searched locally in your browser`;
      els.main.classList.add("active");
      apply();
      if (!state.rawResults.length) renderSequenceEmpty(first.length, Boolean(second));
      const url = new URL(location.href);
      url.search = "";
      history.pushState({}, "", url);
      scrollToResults();
    } catch (error) {
      showSearchError(error.message);
    } finally {
      els.sequenceSearch.disabled = false;
      els.sequenceSearch.textContent = original;
    }
  }

  function renderSequenceEmpty(length, paired) {
    const complete = isCompleteSnapshot();
    els.results.innerHTML = `<div class="empty state-panel">${stateGlyph("empty")}<strong>${complete ? "No exact public sequence match" : "No exact match in this partial snapshot"}</strong><span>${complete ? "PAIRS found no" : "PAIRS cannot conclude absence from an incomplete dataset; it found no"} exact ${paired ? "paired " : ""}match for the normalized ${length}-aa query.</span><div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a></div></div>`;
    els.loadMore.hidden = true;
  }

  function showSearchError(message) {
    state.mode = "sequence";
    state.rawResults = [];
    state.filtered = [];
    els.main.classList.add("active");
    els.targetName.textContent = "Sequence search";
    els.targetMeta.textContent = "Exact local sequence lookup";
    els.summary.innerHTML = "";
    els.results.innerHTML = errorState("Could not run sequence search", message);
    els.loadMore.hidden = true;
    scrollToResults();
  }

  function resetResultState() {
    state.rawResults = [];
    state.filtered = [];
    state.shown = PAGE_RENDER_SIZE;
    state.filter = "all";
    state.targetPageCache.clear();
    state.loadedTargetPages.clear();
    state.loadedFunctionalPages.clear();
    state.loadedNegativePages.clear();
    resetFilterButtons();
  }

  function resetFilterButtons() {
    $$(".filter").forEach(button =>
      button.classList.toggle("active", button.dataset.filter === "all"),
    );
  }

  function updateFilterAvailability() {
    const negative = $('.filter[data-filter="negative"]');
    if (negative)
      negative.hidden =
        state.mode !== "target" ||
        (state.selected?.negative_page_count == null && !(state.selected?.stats?.negative > 0));
    const functional = $('.filter[data-filter="functional"]');
    if (functional)
      functional.hidden =
        state.mode !== "target" ||
        (state.selected?.functional_page_count == null && !(state.selected?.stats?.functional > 0));
  }

  function targetPageUrl(target, pageNumber) {
    return `${DATA_ROOT}/targets/${target.dir}/page-${String(pageNumber).padStart(3, "0")}.json`;
  }

  function negativeTargetPageUrl(target, pageNumber) {
    return `${DATA_ROOT}/targets/${target.dir}/negative-page-${String(pageNumber).padStart(3, "0")}.json`;
  }

  function functionalTargetPageUrl(target, pageNumber) {
    return `${DATA_ROOT}/targets/${target.dir}/functional-page-${String(pageNumber).padStart(3, "0")}.json`;
  }

  function targetPageCategory() {
    if (state.filter === "functional") return "functional";
    if (state.filter === "negative") return "negative";
    return "positive";
  }

  function activeTargetPageCategory() {
    const category = targetPageCategory();
    if (category === "negative" && state.selected?.negative_page_count == null) return "positive";
    if (category === "functional" && state.selected?.functional_page_count == null)
      return "positive";
    return category;
  }

  function usesSeparateNegativePages() {
    return activeTargetPageCategory() === "negative";
  }

  function usesSeparateFunctionalPages() {
    return activeTargetPageCategory() === "functional";
  }

  async function loadTargetPage(pageNumber, category = "positive") {
    const pageCount =
      category === "negative"
        ? state.selected?.negative_page_count || 0
        : category === "functional"
          ? state.selected?.functional_page_count || 0
          : state.selected?.page_count || 0;
    if (!state.selected || pageNumber < 1 || pageNumber > pageCount) return [];
    const pages =
      category === "negative"
        ? state.loadedNegativePages
        : category === "functional"
          ? state.loadedFunctionalPages
          : state.loadedTargetPages;
    const key = `${state.selected.id}:${category}:${pageNumber}`;
    if (state.targetPageCache.has(key)) return state.targetPageCache.get(key);
    const rows = await getJSON(
      category === "negative"
        ? negativeTargetPageUrl(state.selected, pageNumber)
        : category === "functional"
          ? functionalTargetPageUrl(state.selected, pageNumber)
          : targetPageUrl(state.selected, pageNumber),
    );
    state.targetPageCache.set(key, rows);
    if (!pages.has(pageNumber)) {
      pages.add(pageNumber);
      state.rawResults.push(...rows);
    }
    return rows;
  }

  async function loadNextTargetPage(category = activeTargetPageCategory()) {
    if (!state.selected) return false;
    if (category === "negative" && state.selected.negative_page_count == null)
      category = "positive";
    if (category === "functional" && state.selected.functional_page_count == null)
      category = "positive";
    const pageCount =
      category === "negative"
        ? state.selected.negative_page_count || 0
        : category === "functional"
          ? state.selected.functional_page_count || 0
          : state.selected.page_count;
    const pages =
      category === "negative"
        ? state.loadedNegativePages
        : category === "functional"
          ? state.loadedFunctionalPages
          : state.loadedTargetPages;
    for (let page = 1; page <= pageCount; page += 1) {
      if (!pages.has(page)) {
        await loadTargetPage(page, category);
        return true;
      }
    }
    return false;
  }

  async function loadAllTargetPages(category = activeTargetPageCategory()) {
    if (!state.selected) return;
    if (category === "negative" && state.selected.negative_page_count == null)
      category = "positive";
    if (category === "functional" && state.selected.functional_page_count == null)
      category = "positive";
    const pageCount =
      category === "negative"
        ? state.selected.negative_page_count || 0
        : category === "functional"
          ? state.selected.functional_page_count || 0
          : state.selected.page_count;
    const pages =
      category === "negative"
        ? state.loadedNegativePages
        : category === "functional"
          ? state.loadedFunctionalPages
          : state.loadedTargetPages;
    const missing = [];
    for (let page = 1; page <= pageCount; page += 1) {
      if (!pages.has(page)) missing.push(page);
    }
    const rows = await Promise.all(
      missing.map(page =>
        getJSON(
          category === "negative"
            ? negativeTargetPageUrl(state.selected, page)
            : category === "functional"
              ? functionalTargetPageUrl(state.selected, page)
              : targetPageUrl(state.selected, page),
        ),
      ),
    );
    rows.forEach((pageRows, index) => {
      const page = missing[index];
      const key = `${state.selected.id}:${category}:${page}`;
      state.targetPageCache.set(key, pageRows);
      pages.add(page);
      state.rawResults.push(...pageRows);
    });
  }

  async function withViewTransition(work) {
    if (prefersReducedMotion() || typeof document.startViewTransition !== "function") return work();
    let transition;
    try {
      transition = document.startViewTransition(() => work());
    } catch {
      // Older/embedded browsers can expose the API but reject a transition;
      // navigation remains fully functional through the normal render path.
      return work();
    }
    await transition.finished.catch(() => {});
  }

  async function selectTarget(target, push = true) {
    return withViewTransition(() => selectTargetContent(target, push));
  }

  async function selectTargetContent(target, push = true) {
    state.mode = "target";
    state.selected = target;
    resetResultState();
    updateFilterAvailability();
    closeSuggestions();
    els.q.value = target.name;
    els.main.classList.add("active");
    els.targetName.textContent = target.name;
    updateTargetMeta();
    els.summary.innerHTML = "";
    els.results.innerHTML = loadingState("Loading target results…");
    els.loadMore.hidden = true;
    try {
      await loadTargetPage(1, "positive");
      apply();
      if (push) {
        const url = new URL(location.href);
        url.search = "";
        url.searchParams.set("target", target.id);
        history.pushState({}, "", url);
      }
      scrollToResults();
    } catch (error) {
      els.results.innerHTML = errorState("Could not load target results", error.message);
    }
  }

  function updateTargetMeta() {
    if (!state.selected) return;
    const stats = state.selected.stats || {};
    const positiveCount =
      stats.positive_results ??
      stats.positive ??
      state.selected.positive_count ??
      (state.selected.negative_page_count != null
        ? state.selected.result_count
        : stats.negative_only != null
          ? Math.max(0, state.selected.result_count - stats.negative_only)
          : null);
    const countLabel =
      positiveCount == null
        ? `${fmt(state.selected.result_count)} indexed antibodies (negative-only rows excluded)`
        : `${fmt(positiveCount)} antibodies with positive evidence`;
    const loadedLabel = usesSeparateNegativePages()
      ? `loaded ${state.loadedNegativePages.size}/${state.selected.negative_page_count} negative pages`
      : usesSeparateFunctionalPages()
        ? `loaded ${state.loadedFunctionalPages.size}/${state.selected.functional_page_count} functional pages`
        : `loaded ${state.loadedTargetPages.size}/${state.selected.page_count} pages`;
    els.targetMeta.textContent = `${countLabel} · ${fmt(state.selected.count)} source-level relationships · ${state.selected.sources.join(" · ")} · ${loadedLabel}`;
  }

  async function openStandaloneAntibody(antibodyId, push = true) {
    return withViewTransition(() => openStandaloneAntibodyContent(antibodyId, push));
  }

  async function openStandaloneAntibodyContent(antibodyId, push = true) {
    closeSuggestions();
    els.main.classList.add("active");
    els.targetName.textContent = "Antibody record";
    els.targetMeta.textContent = "Loading public sequence record…";
    els.results.innerHTML = loadingState("Loading antibody record…");
    els.summary.innerHTML = "";
    els.loadMore.hidden = true;
    try {
      const antibody = await fetchAntibody(antibodyId);
      if (!antibody) throw new Error("Antibody ID is not present in this snapshot.");
      state.mode = "antibody";
      state.selected = null;
      state.filter = "all";
      state.shown = 1;
      state.rawResults = [
        {
          antibody: summaryFromFull(antibody),
          relationships: [],
          evidence: [],
          sources: antibody.sources || [],
          interactions: [],
        },
      ];
      resetFilterButtons();
      updateFilterAvailability();
      els.targetName.textContent = antibody.name || antibody.id;
      els.targetMeta.textContent = `${fmt(antibody.target_count)} target annotations · ${(antibody.sources || []).join(" · ")} · ${antibody.heavy && antibody.light ? "paired VH + VL" : "sequence record"}`;
      apply();
      const card = els.results.querySelector(".card");
      if (card) {
        card.classList.add("open", "highlight");
        await loadCardDetails(card);
      }
      if (push) {
        const url = new URL(location.href);
        url.search = "";
        url.searchParams.set("ab", antibody.id);
        history.pushState({}, "", url);
      }
      scrollToResults();
    } catch (error) {
      els.results.innerHTML = errorState("Could not open antibody", error.message);
    }
  }

  async function search() {
    const query = els.q.value.trim();
    if (!query) return;
    if (looksLikeSequence(query)) {
      setSearchMode("sequence");
      els.heavySequence.value = query;
      await runSequenceSearch(query);
      return;
    }

    const exact = await exactSearchResolution(query);
    if (exact) {
      closeSuggestions();
      if (exact.kind === "target") await selectTarget(exact.target, true);
      else await openStandaloneAntibody(exact.antibody.id, true);
      return;
    }

    await showSuggestions();
    els.suggestions.innerHTML = `<div class="suggestion no-exact"><div class="s-main"><strong>No exact match for “${esc(query)}”</strong><span>${state.suggestionItems.length ? "Closest indexed matches — not selected. Choose one explicitly if it is what you meant." : "No indexed target or antibody matches this query."}</span></div></div>${els.suggestions.innerHTML}`;
    renderTextNoMatch(query);
  }

  function renderTextNoMatch(query) {
    state.mode = "target";
    state.selected = null;
    state.rawResults = [];
    state.filtered = [];
    els.main.classList.add("active");
    const complete = isCompleteSnapshot();
    els.targetName.textContent = complete ? "No local match" : "No match in partial snapshot";
    els.targetMeta.textContent = complete
      ? `No target or antibody matched “${query}” in the complete indexed snapshot.`
      : `PAIRS cannot conclude that “${query}” is absent because this snapshot is incomplete.`;
    els.summary.innerHTML = "";
    const issue = feedbackIssueUrl(query);
    els.results.innerHTML = `<div class="empty"><strong>${complete ? "No indexed match" : "Absence is unknown"}</strong>${complete ? "Try a synonym, browse the target catalogue, or report the missing query so it can be reviewed for a future alias/source update." : "One or more required sources are missing. Use a complete snapshot before interpreting this as a no-hit."}<div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a>${complete ? (issue ? `<a class="secondary" href="${esc(issue)}" target="_blank" rel="noopener">Report missing target</a>` : `<button class="secondary" data-copy-missing="${esc(query)}">Copy missing-query report</button>`) : ""}</div></div>`;
    els.loadMore.hidden = true;
    scrollToResults();
  }

  function passes(result) {
    const antibody = result.antibody;
    if (state.filter === "negative") return result.relationships.some(isNegative);
    if (state.filter === "functional") return result.relationships.some(isFunctionalPositive);
    if (state.mode === "target" && !result.relationships.some(isPrimaryPositive)) return false;
    if (state.filter === "paired") return Boolean(antibody.has_heavy && antibody.has_light);
    if (state.filter === "therapeutic") return Boolean(antibody.therapeutic_status);
    if (state.filter === "structure") return hasExactStructure(antibody);
    return true;
  }

  function rank(result) {
    const evidence = Math.max(...result.evidence.map(item => evidenceRank[item] ?? 0), 0);
    return (
      evidence * 100 +
      result.sources.length * 10 +
      (state.mode === "target" ? 0 : hasExactStructure(result.antibody) ? 5 : 0) +
      (state.mode === "target" ? 0 : result.antibody.therapeutic_status ? 3 : 0)
    );
  }

  function apply() {
    const before = new Map(
      $$("#results .card[data-ab]").map(card => [
        card.dataset.ab,
        card.getBoundingClientRect().top,
      ]),
    );
    state.filtered = state.rawResults.filter(passes);
    if (els.sort.value === "name") {
      state.filtered.sort((left, right) => left.antibody.name.localeCompare(right.antibody.name));
    } else if (els.sort.value === "sources") {
      state.filtered.sort(
        (left, right) =>
          right.sources.length - left.sources.length ||
          left.antibody.name.localeCompare(right.antibody.name),
      );
    } else {
      state.filtered.sort(
        (left, right) =>
          rank(right) - rank(left) || left.antibody.name.localeCompare(right.antibody.name),
      );
    }
    renderSummary();
    renderResults();
    $$("#results .card[data-ab]").forEach(card => {
      const previous = before.get(card.dataset.ab);
      if (previous == null || !card.animate || prefersReducedMotion()) return;
      const delta = previous - card.getBoundingClientRect().top;
      if (Math.abs(delta) < 1) return;
      card.animate([{ transform: `translateY(${delta}px)` }, { transform: "translateY(0)" }], {
        duration: 240,
        easing: "cubic-bezier(.2,.8,.2,1)",
      });
    });
    updateTargetMeta();
  }

  function renderSummary() {
    if (state.mode === "target" && state.selected?.stats) {
      const stats = state.selected.stats;
      const positiveCount =
        stats.positive_results ??
        stats.positive ??
        state.selected.positive_count ??
        (state.selected.negative_page_count != null
          ? state.selected.result_count
          : stats.negative_only != null
            ? Math.max(0, state.selected.result_count - stats.negative_only)
            : null);
      const collectionCount =
        state.filter === "functional"
          ? (stats.functional ?? state.selected.functional_count ?? 0)
          : state.filter === "negative"
            ? (stats.negative ?? state.selected.negative_count ?? 0)
            : positiveCount == null
              ? stats.unique_results
              : positiveCount;
      const collectionLabel =
        state.filter === "functional"
          ? "Functional activity"
          : state.filter === "negative"
            ? "Negative evidence"
            : positiveCount == null
              ? "Indexed antibodies"
              : "Positive evidence";
      const collectionStats =
        state.filter === "functional"
          ? stats.functional_stats || stats
          : state.filter === "negative"
            ? stats.negative_stats || stats
            : stats;
      els.summary.innerHTML = [
        [collectionCount, collectionLabel],
        [collectionStats.paired, "VH + VL pairs"],
        [
          collectionStats.structure_exact ?? collectionStats.exact_structure ?? 0,
          "Exact PDB structures",
        ],
        [stats.negative, "Negative evidence"],
      ]
        .map(
          ([value, label]) => `<div class="metric"><b>${fmt(value)}</b><span>${label}</span></div>`,
        )
        .join("");
      return;
    }

    const results = state.filtered;
    const paired = results.filter(
      result => result.antibody.has_heavy && result.antibody.has_light,
    ).length;
    const structures = results.filter(result => hasExactStructure(result.antibody)).length;
    const therapeutics = results.filter(result => result.antibody.therapeutic_status).length;
    els.summary.innerHTML = [
      [results.length, "Results"],
      [paired, "VH + VL pairs"],
      [structures, "Exact PDB structures"],
      [therapeutics, "Therapeutics"],
    ]
      .map(
        ([value, label]) => `<div class="metric"><b>${fmt(value)}</b><span>${label}</span></div>`,
      )
      .join("");
  }

  function badges(result) {
    const output = [];
    if (result.antibody.has_heavy && result.antibody.has_light) {
      output.push('<span class="badge paired">VH + VL</span>');
    }
    if (state.mode !== "target" && result.antibody.therapeutic_status) {
      output.push('<span class="badge therapeutic">THERAPEUTIC</span>');
    }
    if (state.mode !== "target") {
      const structures = structureCollections(result.antibody);
      if (structures.exact.length) output.push('<span class="badge structure">EXACT PDB</span>');
      else if (structures.homologous.length)
        output.push('<span class="badge structure">HOMOLOGOUS PDB</span>');
      else if (structures.unknown.length)
        output.push('<span class="badge structure">STRUCTURE · TIER UNKNOWN</span>');
    }
    for (const field of result.match_fields || []) {
      output.push(
        `<span class="badge functional sequence-match">EXACT ${esc(field.toUpperCase())}</span>`,
      );
    }
    for (const relationship of result.relationships.slice(0, 4)) {
      const value = relationshipValue(relationship);
      const semanticClass = isNegative(relationship)
        ? "negative"
        : isLiteratureContext(relationship)
          ? "literature"
          : isFunctionalPositive(relationship)
            ? "functional"
            : isPrimaryPositive(relationship)
              ? "direct"
              : "";
      output.push(
        `<span class="badge ${semanticClass}" data-relationship="${esc(value)}">${esc(relationLabel(relationship).toUpperCase())}</span>`,
      );
    }
    return output.join("");
  }

  function provenanceLink(sourceKey, record = {}) {
    const source = state.manifest?.sources?.[sourceKey];
    const sourceName = source?.name || sourceKey || "Source";
    const recordUrl = record.record_url || record.source_record_url || record.exact_url || "";
    const homepage = source?.homepage || record.source_url || "";
    if (recordUrl) {
      return {
        html: `<a href="${esc(recordUrl)}" target="_blank" rel="noopener">${esc(sourceName)} · exact source record</a>`,
        scope: "record",
      };
    }
    if (homepage) {
      return {
        html: `<a href="${esc(homepage)}" target="_blank" rel="noopener">${esc(sourceName)} · source homepage</a>`,
        scope: "source_homepage",
      };
    }
    return { html: esc(`${sourceName} · source reference`), scope: "unlinked" };
  }

  function evidenceRows(interactions) {
    if (!interactions.length)
      return '<div class="meta">Open the target links below to inspect target-specific evidence.</div>';
    return interactions
      .map(interaction => {
        const provenance = provenanceLink(interaction.source, interaction);
        const details = [
          interaction.evidence,
          interaction.source_record_id &&
            `${provenance.scope === "record" ? "Exact source record" : "Source record ID"}: ${interaction.source_record_id}`,
          interaction.assertion_origin === "derived_hierarchy" && "PAIRS-derived parent target",
          interaction.assertion_origin === "source_epitope" && "Derived from source epitope field",
          interaction.epitope && `Epitope: ${interaction.epitope}`,
          interaction.assay && `Assay: ${interaction.assay}`,
          interaction.reference,
          interaction.note,
        ]
          .filter(Boolean)
          .join(" · ");
        return `<div class="evidence-row ${isNegative(interaction.relationship) ? "negative" : ""}"><strong>${esc(relationLabel(interaction.relationship))} · ${provenance.html}</strong><span>${esc(details)}</span></div>`;
      })
      .join("");
  }

  function renderResults() {
    const subset = state.filtered.slice(0, state.shown);
    const negativeView = usesSeparateNegativePages();
    const functionalView = usesSeparateFunctionalPages();
    const targetHasMorePages =
      state.mode === "target" &&
      state.selected &&
      (negativeView
        ? state.loadedNegativePages.size < (state.selected.negative_page_count || 0)
        : functionalView
          ? state.loadedFunctionalPages.size < (state.selected.functional_page_count || 0)
          : state.loadedTargetPages.size < state.selected.page_count);
    els.loadMore.hidden = state.shown >= state.filtered.length && !targetHasMorePages;

    if (!subset.length) {
      els.results.innerHTML = emptyState("No loaded results match this filter.");
      return;
    }

    els.results.innerHTML = subset
      .map(result => {
        const antibody = result.antibody;
        const alias = (antibody.aliases || []).slice(0, 4).join(" · ");
        const directUrl = buildViewUrl({ ab: antibody.id });
        return `<article class="card" data-ab="${esc(antibody.id)}" data-shard="${esc(antibody.shard)}"><div class="card-main"><div><a class="name-link" href="${esc(directUrl)}" data-ab-link="${esc(antibody.id)}">${esc(antibody.name || antibody.id)}</a><div class="meta">${esc(alias || antibody.organism || antibody.format || "Public antibody record")}</div></div><div class="badges">${badges(result)}</div><div class="sources">${esc(result.sources.join(" · ") || "Public source")}<br>${result.evidence.length ? esc(result.evidence.join(" · ")) : esc((result.match_fields || []).length ? "Exact sequence match" : "Sequence/provenance record")}</div><button class="expand" aria-label="Expand antibody details" aria-expanded="false">+</button></div><div class="detail"><div class="detail-grid"><div class="sequence-slot"><div class="meta">Sequence details load only when this record is expanded.</div></div><div><div class="section-title">Evidence for this view</div><div class="evidence-list">${evidenceRows(result.interactions)}</div><div class="full-record-slot"></div></div></div></div></article>`;
      })
      .join("");
    $$("#results .card").forEach((card, index) =>
      card.style.setProperty("--card-index", Math.min(index, 11)),
    );
  }

  function sequenceBlock(label, value) {
    if (!value) return "";
    return `<div class="seq"><div class="seq-head"><span>${esc(label)} · ${fmt(value.length)} aa</span><button class="copy" data-copy-seq="${esc(value)}" aria-label="Copy ${esc(label)} sequence">Copy</button></div><code class="sequence-residues" aria-label="${esc(value)}">${coloredSequence(value)}</code></div>`;
  }

  function renderConstructContext(antibody) {
    const metadata = antibody.metadata || {};
    const constructs = antibody.constructs || [];
    const construct = antibody.construct || constructs[0] || {};
    const arms = construct.arms || antibody.arms || metadata.arms || [];
    const multispecific = Boolean(
      antibody.multispecific ||
      constructs.length ||
      construct.multispecific ||
      metadata.multispecific ||
      arms.length,
    );
    if (!multispecific) return "";
    const constructTargets = targetNames(
      antibody.construct_targets || construct.targets || metadata.construct_targets || [],
    );
    const context =
      metadata.construct_target_context ||
      antibody.construct_target_context ||
      metadata.quarantine_reason ||
      "Arm-specific target assignment is not available unless explicitly mapped below.";
    const armRows = arms
      .map((arm, index) => {
        const name = arm.designation || arm.name || arm.id || `Arm ${index + 1}`;
        const targets = targetNames(arm.direct_targets || arm.targets || arm.target || []);
        const status = arm.target_assignment_status || arm.assignment_status || "";
        const sequence = [arm.heavy || arm.vh, arm.light || arm.vl].filter(Boolean);
        const sequenceCount =
          sequence.length || Number(Boolean(arm.has_heavy)) + Number(Boolean(arm.has_light));
        return `<div class="evidence-row"><strong>${esc(name)}</strong><span>${sequenceCount ? `${sequenceCount === 2 ? "VH + VL" : "Single-sequence"} arm${status ? ` · ${esc(status)}` : ""}` : "Arm sequence unavailable"}${targets.length ? ` · target context: ${esc(targets.join(" · "))}` : " · target assignment unavailable"}</span></div>`;
      })
      .join("");
    return `<div class="section-title" style="margin-top:16px">Construct / arm context</div><div class="meta">${esc(context)}${constructTargets.length ? ` Construct-level targets: ${esc(constructTargets.join(" · "))}.` : ""}</div>${armRows ? `<div class="evidence-list">${armRows}</div>` : ""}`;
  }

  function renderFullRecordSlot(antibody) {
    const renderTargetCollection = (values, label, link = true) => {
      const entries = targetEntries(values).slice(0, 40);
      const content = entries
        .map(target =>
          link && target.id
            ? `<a class="target-pill" href="${esc(buildViewUrl({ target: target.id }))}" data-target-link="${esc(target.id)}">${esc(target.name)}</a>`
            : `<span class="target-pill">${esc(target.name)}</span>`,
        )
        .join("");
      return content
        ? `<div class="section-title" style="margin-top:16px">${esc(label)}</div><div class="target-pills">${content}</div>`
        : "";
    };
    const directTargets = directTargetNames(antibody).length
      ? antibody.direct_targets || antibody.directTargets || antibody.targets_supported
      : [];
    const legacyTargets = legacyTargetNames(antibody);
    const targets = renderTargetCollection(directTargets, "Direct target evidence");
    const functional = renderTargetCollection(
      antibody.functional_activity ||
        antibody.functionalActivity ||
        antibody.functional_targets ||
        antibody.functionalTargets,
      "Functional activity",
    );
    const negative = renderTargetCollection(
      antibody.negative_evidence || antibody.negativeEvidence,
      "Negative evidence",
    );
    const literature = renderTargetCollection(
      antibody.literature_mentions || antibody.literatureMentions,
      "Literature context",
      false,
    );
    const legacy =
      !directTargets.length && legacyTargets.length
        ? renderTargetCollection(legacyTargets, "Legacy associated target annotations", false)
        : "";
    const structures = structureCollections(antibody);
    const renderStructures = (entries, label) => {
      if (!entries.length) return "";
      const links = entries
        .map(entry => {
          const suffix = entry.identityLabel
            ? ` · ${esc(entry.identityLabel)} SI`
            : entry.identity
              ? ` · ${esc(entry.identity)}% SI`
              : "";
          return `<a class="pdb-link" href="https://www.rcsb.org/structure/${encodeURIComponent(entry.id)}" target="_blank" rel="noopener">${esc(entry.id)}${suffix}</a>`;
        })
        .join("");
      return `<div class="section-title" style="margin-top:16px">${esc(label)}</div><div class="target-pills">${links}</div>`;
    };
    const structureHtml =
      renderStructures(structures.exact, "Exact PDB structures") +
      renderStructures(structures.homologous, "Homologous PDB structures") +
      renderStructures(structures.unknown, "Structures · identity tier unavailable");
    const structureContext =
      structureHtml && state.mode === "target"
        ? '<div class="meta" style="margin-top:16px">Structure metadata is sequence-level; target evidence for this view is listed above.</div>'
        : "";
    const provenance = (antibody.source_records || [])
      .map(record => {
        const label = state.manifest?.sources?.[record.source]?.name || record.source;
        const provenance = provenanceLink(record.source, record);
        const scopeLabel =
          record.link_scope === "record" || provenance.scope === "record"
            ? "Exact source record"
            : "Source record ID; source homepage";
        return `<div class="evidence-row"><strong>${provenance.html} · ${esc(record.record_id)}</strong><span>${esc(scopeLabel)} · ${esc(record.reference || "Source record")}</span></div>`;
      })
      .join("");
    return `${renderConstructContext(antibody)}${targets || '<div class="meta" style="margin-top:16px">No direct target evidence stored.</div>'}${functional}${negative}${literature}${legacy}${structureContext}${structureHtml}<div class="section-title" style="margin-top:16px">Record provenance</div><div class="evidence-list">${provenance || '<div class="meta">No source-record details stored.</div>'}</div><div class="detail-actions"><button class="secondary" data-copy-ab-url="${esc(antibody.id)}">Copy antibody URL</button></div>`;
  }

  async function loadCardDetails(card) {
    if (card.dataset.loaded === "true") return;
    const antibody = await fetchAntibody(card.dataset.ab, card.dataset.shard);
    if (!antibody) return;
    const sequenceSlot = card.querySelector(".sequence-slot");
    sequenceSlot.innerHTML = `<div class="section-title">Sequences</div>${sequenceBlock("VH / VHH", antibody.heavy)}${sequenceBlock("VL", antibody.light)}${sequenceBlock("CDRH3", antibody.cdrh3)}${sequenceBlock("CDRL3", antibody.cdrl3)}${!antibody.heavy && !antibody.light ? '<div class="meta">No variable-domain sequence is available in the imported record.</div>' : ""}`;
    card.querySelector(".full-record-slot").innerHTML = renderFullRecordSlot(antibody);
    card.dataset.loaded = "true";
  }

  async function ensureFilteredRows(minimum = PAGE_RENDER_SIZE) {
    if (state.mode !== "target" || !state.selected) return;
    const category = activeTargetPageCategory();
    const loadedPages =
      category === "negative"
        ? state.loadedNegativePages
        : category === "functional"
          ? state.loadedFunctionalPages
          : state.loadedTargetPages;
    const pageCount =
      category === "negative"
        ? state.selected.negative_page_count || 0
        : category === "functional"
          ? state.selected.functional_page_count || 0
          : state.selected.page_count;
    while (state.filtered.length < minimum && loadedPages.size < pageCount) {
      const loaded = await loadNextTargetPage(category);
      if (!loaded) break;
      state.filtered = state.rawResults.filter(passes);
    }
  }

  async function handleFilter(filter) {
    state.filter = filter;
    state.shown = PAGE_RENDER_SIZE;
    $$(".filter").forEach(button =>
      button.classList.toggle("active", button.dataset.filter === filter),
    );
    state.filtered = state.rawResults.filter(passes);
    await ensureFilteredRows(PAGE_RENDER_SIZE);
    apply();
  }

  async function handleSort() {
    if (
      state.mode === "target" &&
      state.selected &&
      els.sort.value !== "evidence" &&
      (usesSeparateNegativePages()
        ? state.loadedNegativePages.size < (state.selected.negative_page_count || 0)
        : usesSeparateFunctionalPages()
          ? state.loadedFunctionalPages.size < (state.selected.functional_page_count || 0)
          : state.loadedTargetPages.size < state.selected.page_count)
    ) {
      const previous = els.sort.disabled;
      els.sort.disabled = true;
      els.targetMeta.textContent = "Loading all target pages for a global sort…";
      try {
        await loadAllTargetPages();
      } finally {
        els.sort.disabled = previous;
      }
    }
    apply();
  }

  async function handleLoadMore() {
    if (state.shown < state.filtered.length) {
      state.shown += PAGE_RENDER_SIZE;
      renderResults();
      return;
    }
    if (state.mode === "target" && state.selected) {
      const category = activeTargetPageCategory();
      const loadedPages =
        category === "negative"
          ? state.loadedNegativePages
          : category === "functional"
            ? state.loadedFunctionalPages
            : state.loadedTargetPages;
      const pageCount =
        category === "negative"
          ? state.selected.negative_page_count || 0
          : category === "functional"
            ? state.selected.functional_page_count || 0
            : state.selected.page_count;
      const before = state.filtered.length;
      while (loadedPages.size < pageCount) {
        const loaded = await loadNextTargetPage(category);
        if (!loaded) break;
        state.filtered = state.rawResults.filter(passes);
        if (state.filtered.length > before) break;
      }
      state.shown += PAGE_RENDER_SIZE;
      apply();
    }
  }

  async function allRowsForCurrentView() {
    if (state.mode === "target" && state.selected) {
      await loadAllTargetPages(activeTargetPageCategory());
      apply();
    }
    return state.rawResults.filter(passes);
  }

  function csvEscape(value) {
    const string = String(value ?? "");
    return /[",\n]/.test(string) ? `"${string.replaceAll('"', '""')}"` : string;
  }

  function downloadBlob(filename, content, type) {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = filename;
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function exportBaseName() {
    const label =
      state.selected?.name || (state.mode === "sequence" ? "sequence-search" : "antibody");
    return `pairs-${label}`
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, "-")
      .replace(/^-|-$/g, "");
  }

  async function withBusyButton(button, busyLabel, callback) {
    const label = button.textContent;
    button.disabled = true;
    button.textContent = busyLabel;
    try {
      await callback();
    } finally {
      button.disabled = false;
      button.textContent = label;
    }
  }

  async function exportCSV() {
    await withBusyButton(els.csv, "Preparing…", async () => {
      const rows = await allRowsForCurrentView();
      const fullRecords = await fetchFullRecords(rows.map(row => row.antibody.id));
      const fields = [
        "name",
        "antibody_id",
        "organism",
        "format",
        "heavy",
        "light",
        "cdrh3",
        "cdrl3",
        "therapeutic_status",
        "exact_structures",
        "homologous_structures",
        "structure_tier_unknown",
        "sources",
        "direct_targets",
        "functional_activity",
        "negative_evidence",
        "literature_mentions",
        "legacy_associated_target_annotations",
        "relationships",
        "evidence",
      ];
      const output = [fields.join(",")];
      for (const row of rows) {
        const antibody = fullRecords.get(row.antibody.id) || {};
        const record = {
          name: antibody.name || row.antibody.name,
          antibody_id: row.antibody.id,
          organism: antibody.organism || "",
          format: antibody.format || "",
          heavy: antibody.heavy || "",
          light: antibody.light || "",
          cdrh3: antibody.cdrh3 || "",
          cdrl3: antibody.cdrl3 || "",
          therapeutic_status: antibody.therapeutic_status || "",
          exact_structures: structureCollections(antibody)
            .exact.map(item => item.id)
            .join(";"),
          homologous_structures: structureCollections(antibody)
            .homologous.map(item =>
              item.identityLabel
                ? `${item.id} (${item.identityLabel} SI)`
                : item.identity
                  ? `${item.id} (${item.identity}% SI)`
                  : item.id,
            )
            .join(";"),
          structure_tier_unknown: structureCollections(antibody)
            .unknown.map(item => item.id)
            .join(";"),
          sources: (antibody.sources || row.sources || []).join(";"),
          direct_targets: directTargetNames(antibody).join(";"),
          functional_activity: functionalActivityNames(antibody).join(";"),
          negative_evidence: negativeEvidenceNames(antibody).join(";"),
          literature_mentions: literatureMentionNames(antibody).join(";"),
          legacy_associated_target_annotations: directTargetNames(antibody).length
            ? ""
            : legacyTargetNames(antibody).join(";"),
          relationships: row.relationships.join(";"),
          evidence: row.evidence.join(";"),
        };
        output.push(fields.map(field => csvEscape(record[field])).join(","));
      }
      downloadBlob(`${exportBaseName()}.csv`, output.join("\n"), "text/csv;charset=utf-8");
    });
  }

  async function exportFASTA() {
    await withBusyButton(els.fasta, "Preparing…", async () => {
      const rows = await allRowsForCurrentView();
      const fullRecords = await fetchFullRecords(rows.map(row => row.antibody.id));
      const output = [];
      for (const row of rows) {
        const antibody = fullRecords.get(row.antibody.id) || {};
        const safeName = String(antibody.name || row.antibody.name || row.antibody.id).replace(
          /\s+/g,
          "_",
        );
        if (antibody.heavy) output.push(`>${safeName}|${row.antibody.id}|VH\n${antibody.heavy}`);
        if (antibody.light) output.push(`>${safeName}|${row.antibody.id}|VL\n${antibody.light}`);
      }
      downloadBlob(`${exportBaseName()}.fasta`, output.join("\n"), "text/plain;charset=utf-8");
    });
  }

  function buildViewUrl(parameters) {
    const url = new URL(location.href);
    url.search = "";
    for (const [key, value] of Object.entries(parameters)) {
      if (value) url.searchParams.set(key, value);
    }
    return `${url.pathname}${url.search}`;
  }

  async function copyText(value, button = null) {
    try {
      await navigator.clipboard.writeText(value);
      if (button) {
        const original = button.textContent;
        const originalLabel = button.getAttribute("aria-label");
        button.classList.add("copied");
        button.textContent = "Copied";
        button.setAttribute("aria-label", "Copied");
        setTimeout(() => {
          button.textContent = original;
          button.classList.remove("copied");
          if (originalLabel) button.setAttribute("aria-label", originalLabel);
          else button.removeAttribute("aria-label");
        }, 900);
      }
    } catch {
      window.prompt("Copy:", value);
    }
  }

  function feedbackIssueUrl(query) {
    if (!location.hostname.endsWith(".github.io")) return "";
    const owner = location.hostname.split(".")[0];
    const repo = location.pathname.split("/").filter(Boolean)[0];
    if (!owner || !repo) return "";
    const title = `Missing PAIRS search result: ${query}`;
    const body = `I searched PAIRS for:\n\n${query}\n\nNo local match was returned in snapshot ${state.manifest?.snapshot || "unknown"}.\n\nPossible synonym/source to review:\n`;
    return `https://github.com/${encodeURIComponent(owner)}/${encodeURIComponent(repo)}/issues/new?title=${encodeURIComponent(title)}&body=${encodeURIComponent(body)}`;
  }

  function fieldHash(value) {
    let hash = 2166136261;
    for (const character of String(value || "")) {
      hash ^= character.charCodeAt(0);
      hash = Math.imul(hash, 16777619);
    }
    return hash >>> 0;
  }

  function targetChannel(target) {
    const stats = target.stats || {};
    if ((target.result_count || 0) === 0) {
      if ((stats.functional || target.functional_count || 0) > 0) return "var(--field-primary)";
      if ((target.negative_count || stats.negative || 0) > 0) return "var(--field-negative)";
    }
    if ((stats.therapeutic || 0) > 0) return "var(--field-therapeutic)";
    const resultCount = target.result_count || stats.unique_results || 0;
    const sequenceEvidence = Math.max(
      stats.paired || 0,
      (stats.structure_exact || 0) + (stats.structure_homologous || 0),
    );
    if (resultCount > 0 && sequenceEvidence >= resultCount / 2) return "var(--field-structural)";
    return "var(--field-primary)";
  }

  function targetRadius(target, scale = 1) {
    const count = target.result_count || target.count || 1;
    const structure = target.stats?.structure_exact || 0;
    return Math.max(
      3.5,
      Math.min(15, (3.5 + Math.sqrt(count) * 0.58 + (structure ? 1.5 : 0)) * scale),
    );
  }

  function yGlyphPath(radius) {
    return `M ${-radius} ${-radius * 0.85} L 0 0 L ${radius} ${-radius * 0.85} M 0 0 V ${radius * 1.22}`;
  }

  function renderHeroField() {
    if (!els.heroSvg) return;
    const targets = [...state.targets]
      .sort(
        (left, right) =>
          (right.result_count || right.count || 0) - (left.result_count || left.count || 0) ||
          left.name.localeCompare(right.name),
      )
      .slice(0, 50);
    const points = targets.map((target, index) => {
      const hash = fieldHash(target.id);
      const x = 24 + ((hash % 1152 || index * 21) % 1152);
      const y = 42 + ((hash >>> 9) % 330 || (index * 47) % 330);
      return {
        target,
        index,
        x,
        y,
        radius: targetRadius(target, 0.48),
        color: targetChannel(target),
      };
    });
    const anchors = points
      .map(
        ({ target, index, x, y, radius, color }) =>
          `<path class="hero-anchor" d="${yGlyphPath(radius)}" transform="translate(${x} ${y})" stroke="${color}" data-hero-index="${index}" />`,
      )
      .join("");
    const mobile = window.matchMedia?.("(max-width: 620px)").matches;
    const decorative = points
      .filter((_, index) => index % Math.max(1, Math.ceil(points.length / (mobile ? 8 : 14))) === 0)
      .slice(0, mobile ? 8 : 14)
      .map(({ target, index, x, y, color }) => {
        const hash = fieldHash(`${target.id}:hero-free`);
        const freeX = 20 + (hash % 1160);
        const freeY = 24 + ((hash >>> 9) % 360);
        return `<path class="hero-ambient-y" d="${yGlyphPath(2.7)}" transform="translate(${freeX} ${freeY})" stroke="${color}" style="--node-color:${color}" data-hero-free-x="${freeX}" data-hero-free-y="${freeY}" data-hero-dock-x="${x}" data-hero-dock-y="${y}" data-hero-index="${index}" />`;
      })
      .join("");
    els.heroSvg.innerHTML = `${anchors}${decorative}`;
    if (state.heroObserver) state.heroObserver.disconnect();
    const animate = () => {
      const reduced = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
      if (!reduced) {
        $$("#heroAmbientSvg .hero-ambient-y").forEach((glyph, index) => {
          const freeX = Number(glyph.dataset.heroFreeX);
          const freeY = Number(glyph.dataset.heroFreeY);
          const dockX = Number(glyph.dataset.heroDockX);
          const dockY = Number(glyph.dataset.heroDockY);
          glyph.animate(
            [
              { transform: `translate(${freeX} ${freeY}) scale(0.72)`, opacity: 0.03 },
              { transform: `translate(${dockX} ${dockY}) scale(1)`, opacity: 0.18, offset: 0.38 },
              { transform: `translate(${dockX} ${dockY}) scale(1.35)`, opacity: 0.3, offset: 0.52 },
              { transform: `translate(${dockX} ${dockY}) scale(1)`, opacity: 0.17, offset: 0.66 },
              { transform: `translate(${freeX} ${freeY}) scale(0.72)`, opacity: 0.03 },
            ],
            {
              duration: 9800 + (index % 7) * 500,
              delay: index * 95,
              iterations: Infinity,
              easing: "ease-in-out",
            },
          );
        });
      }
    };
    if ("IntersectionObserver" in window) {
      state.heroObserver = new IntersectionObserver(
        entries => {
          if (entries.some(entry => entry.isIntersecting)) {
            animate();
            state.heroObserver.disconnect();
          }
        },
        { threshold: 0.1 },
      );
      state.heroObserver.observe(els.heroSvg);
    } else animate();
  }

  function fieldPoint(index, total, id) {
    const columns = Math.max(1, Math.ceil(Math.sqrt(total * 1.45)));
    const rows = Math.max(1, Math.ceil(total / columns));
    const column = index % columns;
    const row = Math.floor(index / columns);
    const hash = fieldHash(id);
    const jitterX = ((hash & 255) / 255 - 0.5) * 0.36;
    const jitterY = (((hash >>> 8) & 255) / 255 - 0.5) * 0.36;
    return {
      x: Math.max(24, Math.min(976, ((column + 0.5 + jitterX) / columns) * 1000)),
      y: Math.max(24, Math.min(536, ((row + 0.5 + jitterY) / rows) * 560)),
    };
  }

  function showFieldTooltip(node) {
    const target = state.targets.find(item => item.id === node.dataset.fieldTarget);
    if (!target || !els.fieldTooltip) return;
    els.fieldTooltip.innerHTML = `<strong>${esc(target.name)}</strong><span>${fmt(target.result_count || target.count || 0)} positive-evidence antibodies · ${fmt((target.sources || []).length)} sources</span>`;
    els.fieldTooltip.style.left = `${Number(node.dataset.fieldX) / 10}%`;
    els.fieldTooltip.style.top = `${Number(node.dataset.fieldY) / 5.6}%`;
    els.fieldTooltip.hidden = false;
  }

  function hideFieldTooltip() {
    if (els.fieldTooltip) els.fieldTooltip.hidden = true;
  }

  function animateFieldNodes() {
    if (!els.bindingFieldSvg || window.matchMedia?.("(prefers-reduced-motion: reduce)").matches)
      return;
    $$("#bindingFieldSvg .field-ambient-y").forEach((node, index) => {
      const freeX = Number(node.dataset.fieldFreeX);
      const freeY = Number(node.dataset.fieldFreeY);
      const dockX = Number(node.dataset.fieldDockX);
      const dockY = Number(node.dataset.fieldDockY);
      node.animate(
        [
          { opacity: 0, transform: `translate(${freeX} ${freeY}) scale(0.35)` },
          { opacity: 0.92, transform: `translate(${dockX} ${dockY}) scale(1)`, offset: 0.38 },
          { opacity: 1, transform: `translate(${dockX} ${dockY}) scale(1.3)`, offset: 0.52 },
          { opacity: 0.82, transform: `translate(${dockX} ${dockY}) scale(1)`, offset: 0.66 },
          { opacity: 0, transform: `translate(${freeX} ${freeY}) scale(0.35)` },
        ],
        {
          duration: 7600 + (index % 9) * 240,
          delay: Math.min(index * 5, 420),
          iterations: Infinity,
          easing: "ease-in-out",
          fill: "both",
        },
      );
    });
  }

  function renderBindingField(targets) {
    if (!els.bindingField || !els.bindingFieldSvg || state.browseView !== "field") return;
    const mobile = window.matchMedia?.("(max-width: 620px)").matches;
    const limit = mobile ? Math.min(targets.length, 72) : targets.length;
    const nodes = targets.slice(0, limit);
    const points = nodes.map((target, index) => ({
      target,
      ...fieldPoint(index, nodes.length, target.id),
    }));
    const links = points
      .filter((_, index) => index > 0 && index % 2 === 0)
      .map((point, index) => {
        const previous = points[index * 2 + 1];
        return previous
          ? `<line class="field-link" x1="${point.x}" y1="${point.y}" x2="${previous.x}" y2="${previous.y}" />`
          : "";
      })
      .join("");
    const decorativeLimit = mobile ? 8 : 14;
    const decorative = points
      .filter((_, index) => index % Math.max(1, Math.ceil(points.length / decorativeLimit)) === 0)
      .slice(0, decorativeLimit)
      .map(({ target, x, y }, index) => {
        const hash = fieldHash(`${target.id}:decorative`);
        const dx = (hash % 44) - 22;
        const dy = ((hash >>> 8) % 30) - 15;
        const dockX = x;
        const dockY = y;
        const freeX = Math.max(18, Math.min(982, x + dx * 5));
        const freeY = Math.max(18, Math.min(542, y + dy * 5));
        const color = targetChannel(target);
        return `<path class="field-ambient-y" d="${yGlyphPath(3.2)}" transform="translate(${freeX} ${freeY})" stroke="${color}" style="--node-color:${color}" data-field-free-x="${freeX}" data-field-free-y="${freeY}" data-field-dock-x="${dockX}" data-field-dock-y="${dockY}" aria-hidden="true" />`;
      })
      .join("");
    const circles = points
      .map(({ target, x, y }) => {
        const radius = targetRadius(target);
        const color = targetChannel(target);
        const label = `${target.name}; ${fmt(target.result_count || target.count || 0)} positive-evidence antibodies`;
        return `<g class="field-node" tabindex="0" role="button" aria-label="${esc(label)}" data-field-target="${esc(target.id)}" data-field-x="${x}" data-field-y="${y}" style="--node-color:${color}" transform="translate(${x} ${y})"><path d="${yGlyphPath(radius)}" /><title>${esc(label)}</title></g>`;
      })
      .join("");
    els.bindingFieldSvg.innerHTML = `${links}${decorative}${circles}`;
    if (state.fieldObserver) state.fieldObserver.disconnect();
    if ("IntersectionObserver" in window) {
      state.fieldObserver = new IntersectionObserver(
        entries => {
          if (entries.some(entry => entry.isIntersecting)) {
            animateFieldNodes();
            state.fieldObserver.disconnect();
          }
        },
        { threshold: 0.15 },
      );
      state.fieldObserver.observe(els.bindingField);
    } else {
      animateFieldNodes();
    }
  }

  function setBrowseView(view) {
    state.browseView = view === "field" ? "field" : "grid";
    const field = state.browseView === "field";
    els.browseGrid?.classList.toggle("active", !field);
    els.browseField?.classList.toggle("active", field);
    els.browseGrid?.setAttribute("aria-pressed", String(!field));
    els.browseField?.setAttribute("aria-pressed", String(field));
    els.targetGrid.hidden = field;
    els.bindingField.hidden = !field;
    const visibleView = field ? els.bindingField : els.targetGrid;
    visibleView.classList.remove("view-enter");
    requestAnimationFrame(() => visibleView.classList.add("view-enter"));
    renderBrowse();
  }

  function renderBrowse() {
    const query = norm(els.browseQuery.value);
    const facets = new Set($$("#browseFacets input:checked").map(input => input.dataset.facet));
    const facetCount = (target, facet) => {
      if (facet === "structure") {
        return (
          target.stats?.structure_exact ??
          target.stats?.exact_structure ??
          target.exact_structure_count ??
          0
        );
      }
      return target.stats?.[facet] || 0;
    };
    let targets = state.targets.filter(target => {
      if ([...facets].some(facet => !(facetCount(target, facet) > 0))) return false;
      if (!query) return true;
      return [target.name, ...(target.aliases || [])].some(value => norm(value).includes(query));
    });

    if (els.browseSort.value === "name") {
      targets.sort((left, right) => left.name.localeCompare(right.name));
    } else if (els.browseSort.value === "sources") {
      targets.sort(
        (left, right) =>
          right.sources.length - left.sources.length ||
          right.result_count - left.result_count ||
          left.name.localeCompare(right.name),
      );
    } else {
      targets.sort(
        (left, right) =>
          right.result_count - left.result_count || left.name.localeCompare(right.name),
      );
    }

    const visible = targets.slice(0, state.browseShown);
    els.browseMeta.textContent = `${fmt(targets.length)} targets shown from ${fmt(state.targets.length)} indexed targets`;
    els.targetGrid.innerHTML = visible
      .map(target => {
        const aliases = (target.aliases || [])
          .filter(alias => norm(alias) !== norm(target.name))
          .slice(0, 3)
          .join(" · ");
        return `<button class="target-card" data-browse-target="${esc(target.id)}"><span><strong>${esc(target.name)}</strong><small>${esc(aliases || target.sources.join(" · "))}</small></span><span class="target-count">${fmt(target.result_count)} Abs<br>${target.sources.length} src</span></button>`;
      })
      .join("");
    $$("#targetGrid .target-card").forEach((card, index) =>
      card.style.setProperty("--card-index", Math.min(index, 11)),
    );
    els.browseMore.hidden = state.browseView === "field" || visible.length >= targets.length;
    renderBindingField(targets);
  }

  function parseBatchFasta(value) {
    const records = [];
    let header = "";
    let lines = [];
    const push = () => {
      if (!header && !lines.length) return;
      const query_id = header || `sequence_${records.length + 1}`;
      try {
        const normalized_sequence = normalizePastedSequence(lines.join("\n"));
        records.push({
          query_id,
          normalized_sequence,
          length: normalized_sequence.length,
          validation_status: "valid",
        });
      } catch (error) {
        records.push({
          query_id,
          normalized_sequence: "",
          length: 0,
          validation_status: error.message,
        });
      }
    };
    for (const line of String(value ?? "").split(/\r?\n/)) {
      if (line.trim().startsWith(">")) {
        push();
        header = line.trim().slice(1).trim() || `sequence_${records.length + 1}`;
        lines = [];
      } else if (line.trim() || header) lines.push(line);
    }
    push();
    return records;
  }
  function renderBatchRows() {
    els.batchResults.innerHTML = state.batchRows
      .map(
        row =>
          `<div class="batch-row"><span>${esc(row.query_id)}</span><strong>${esc(row.match_type)}</strong><span>${row.antibody_id ? `<a href="?ab=${esc(row.antibody_id)}">${esc(row.antibody_name)}</a>` : "—"}</span><span>${esc(row.matched_sequence || row.validation_status || "—")}</span></div>`,
      )
      .join("");
  }
  async function runBatchSearch() {
    const records = parseBatchFasta(els.batchQuery.value),
      type = els.sequenceType.value,
      unique = new Map(
        records
          .filter(row => row.validation_status === "valid")
          .map(row => [row.normalized_sequence, row]),
      );
    els.batchSearch.disabled = true;
    els.batchProgress.textContent = `Parsed ${records.length} sequences; searching ${unique.size} unique sequences…`;
    try {
      const found = new Map();
      await Promise.all(
        [...unique.values()].map(async row => {
          const hits = ["cdrh3", "cdrl3"].includes(type)
            ? await lookupCdr(row.normalized_sequence, type, els.nearMatches.checked)
            : (await lookupSequence(row.normalized_sequence))
                .filter(hit => type === "auto" || hit.field === type)
                .map(hit => ({
                  antibody_ids: [hit.id],
                  sequence: row.normalized_sequence,
                  distance: 0,
                }));
          found.set(
            row.normalized_sequence,
            hits.flatMap(hit =>
              hit.antibody_ids.map(id => ({ id, sequence: hit.sequence, distance: hit.distance })),
            ),
          );
        }),
      );
      const full = await fetchFullRecords([
        ...new Set([...found.values()].flatMap(hits => hits.map(hit => hit.id))),
      ]);
      state.batchRows = records.flatMap(row => {
        if (row.validation_status !== "valid") return [{ ...row, match_type: "INVALID" }];
        const hits = found.get(row.normalized_sequence) || [];
        return hits.length
          ? hits.map(hit => ({
              ...row,
              match_type: hit.distance ? "NEAR" : "EXACT",
              matched_sequence: hit.sequence,
              antibody_id: hit.id,
              antibody_name: full.get(hit.id)?.name || hit.id,
              edit_distance: hit.distance,
              direct_targets: directTargetNames(full.get(hit.id) || {}).join("; "),
              legacy_associated_target_annotations: directTargetNames(full.get(hit.id) || {}).length
                ? ""
                : legacyTargetNames(full.get(hit.id) || {}).join("; "),
              sources: (full.get(hit.id)?.sources || []).join("; "),
            }))
          : [{ ...row, match_type: "NO MATCH" }];
      });
      els.batchProgress.textContent = `${records.length} / ${records.length} complete`;
      els.batchCsv.disabled = false;
      renderBatchRows();
    } catch (error) {
      els.batchProgress.textContent = error.message;
    } finally {
      els.batchSearch.disabled = false;
    }
  }
  function downloadBatchCsv() {
    const columns = [
      "query_id",
      "query_sequence",
      "query_length",
      "query_type",
      "match_type",
      "matched_sequence",
      "edit_distance",
      "antibody_id",
      "antibody_name",
      "direct_targets",
      "legacy_associated_target_annotations",
      "sources",
      "deep_link",
      "snapshot_date",
      "schema_version",
    ];
    const quote = value => `"${String(value ?? "").replaceAll('"', '""')}"`;
    const text = [
      columns.join(","),
      ...state.batchRows.map(row =>
        columns
          .map(column =>
            quote(
              {
                query_id: row.query_id,
                query_sequence: row.normalized_sequence,
                query_length: row.length,
                query_type: els.sequenceType.value,
                match_type: row.match_type,
                matched_sequence: row.matched_sequence,
                edit_distance: row.edit_distance,
                antibody_id: row.antibody_id,
                antibody_name: row.antibody_name,
                direct_targets: row.direct_targets,
                legacy_associated_target_annotations: row.legacy_associated_target_annotations,
                sources: row.sources,
                deep_link: row.antibody_id
                  ? `${location.origin}${location.pathname}?ab=${row.antibody_id}`
                  : "",
                snapshot_date: state.manifest?.snapshot_date || state.manifest?.snapshot,
                schema_version: SUPPORTED_SCHEMA,
              }[column],
            ),
          )
          .join(","),
      ),
    ].join("\n");
    const link = document.createElement("a");
    link.href = URL.createObjectURL(new Blob([text], { type: "text/csv" }));
    link.download = "pairs-v4-screening.csv";
    link.click();
    URL.revokeObjectURL(link.href);
  }

  function renderSourceStatus() {
    if (!state.manifest) return;
    els.sourceList.innerHTML = Object.entries(state.manifest.sources)
      .map(([key, source]) => {
        const details = source.ok
          ? `${fmt(source.records)} records · ${fmt(source.interactions)} evidence rows${source.bytes ? ` · ${formatBytes(source.bytes)}` : ""}${source.last_modified ? ` · upstream last modified ${esc(formatDate(source.last_modified))}` : ""}${source.discovered ? " · current link discovered at build time" : ""}`
          : source.error || "Import failed";
        return `<div class="source-row"><div><strong>${esc(source.name || key)}</strong><div class="meta">${esc(key)}</div></div><div>${source.homepage ? `<a href="${esc(source.homepage)}" target="_blank" rel="noopener">Upstream source</a>` : ""}<div class="meta">${esc(source.license || "")}</div><div class="meta">${esc(details)}</div></div><div class="${source.ok ? "good" : "bad"}">${source.ok ? "Included" : "Failed"}</div></div>`;
      })
      .join("");
  }

  function formatBytes(value) {
    if (!value) return "0 B";
    const units = ["B", "KB", "MB", "GB"];
    let number = value;
    let unit = 0;
    while (number >= 1024 && unit < units.length - 1) {
      number /= 1024;
      unit += 1;
    }
    return `${number.toFixed(unit ? 1 : 0)} ${units[unit]}`;
  }

  function formatDate(value, withTime = false) {
    if (!value) return "unknown";
    const date = new Date(value);
    if (Number.isNaN(date.getTime())) return String(value);
    return new Intl.DateTimeFormat(undefined, {
      dateStyle: "medium",
      ...(withTime ? { timeStyle: "short" } : {}),
    }).format(date);
  }

  function renderStatus() {
    const stats = state.manifest.stats;
    const sourceOkay =
      state.manifest.sources_ok ??
      Object.values(state.manifest.sources).filter(source => source.ok).length;
    const sourceExpected =
      state.manifest.sources_expected ?? Object.keys(state.manifest.sources).length;
    const partial = sourceOkay < sourceExpected;
    const indexed = state.manifest.snapshot_date || state.manifest.snapshot;
    els.status.innerHTML = `<span><span class="status-dot ${partial ? "warning" : ""}"></span>${fmt(stats.antibodies)} antibodies · ${fmt(stats.interactions)} evidence records · ${fmt(stats.targets)} targets</span><span>${sourceOkay}/${sourceExpected} public sources included</span><span>Indexed ${esc(formatDate(indexed))}</span>`;
    els.footer.textContent = `PAIRS indexed ${formatDate(indexed, true)} · schema v${state.manifest.schema_version}`;

    if (partial) {
      els.warning.hidden = false;
      els.warning.innerHTML = `This snapshot is partial: ${sourceOkay} of ${sourceExpected} expected public sources imported successfully. Results may be incomplete. <button id="warningSourceBtn">View source status</button>`;
      $("#warningSourceBtn")?.addEventListener("click", () => els.modal.classList.add("open"));
    } else {
      els.warning.hidden = true;
    }
  }

  function showSchemaError(manifest) {
    els.status.innerHTML = `<span class="bad">This PAIRS frontend supports data schema v${SUPPORTED_SCHEMA}, but the deployed manifest reports v${esc(manifest.schema_version)}. Reload the page; if this persists, deploy the matching frontend and versioned data path together.</span>`;
  }

  function scrollToResults() {
    window.scrollTo({ top: Math.max(0, els.main.offsetTop - 55), behavior: "smooth" });
  }

  async function init() {
    try {
      state.manifest = await getJSON(`${DATA_ROOT}/manifest.json`);
      if (state.manifest.schema_version !== SUPPORTED_SCHEMA) {
        showSchemaError(state.manifest);
        return;
      }
      state.targets = await getJSON(`${DATA_ROOT}/targets.json`);
      renderStatus();
      renderSourceStatus();
      renderHeroField();
      renderBrowse();

      const parameters = new URLSearchParams(location.search);
      const antibodyId = parameters.get("ab");
      const targetId = parameters.get("target");
      const query = parameters.get("q");
      if (antibodyId) {
        await openStandaloneAntibody(antibodyId, false);
      } else if (targetId) {
        const target = state.targets.find(item => item.id === targetId);
        if (target) await selectTarget(target, false);
      } else if (query) {
        els.q.value = query;
        await search();
      }
    } catch (error) {
      els.status.innerHTML = `<span class="bad">Dataset failed to load: ${esc(error.message)} Serve the repository over HTTP rather than opening index.html directly.</span>`;
      console.error(error);
    }
  }

  els.q.addEventListener("input", showSuggestions);
  els.q.addEventListener("keydown", event => {
    if (event.key === "ArrowDown") {
      event.preventDefault();
      if (!els.suggestions.classList.contains("open")) showSuggestions();
      setActiveSuggestion(state.activeSuggestion + 1);
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      setActiveSuggestion(
        state.activeSuggestion <= 0 ? state.suggestionItems.length - 1 : state.activeSuggestion - 1,
      );
    } else if (event.key === "Enter") {
      event.preventDefault();
      if (els.suggestions.classList.contains("open") && state.activeSuggestion >= 0) {
        activateSuggestion(state.activeSuggestion);
      } else {
        search();
      }
    } else if (event.key === "Escape") {
      closeSuggestions();
    }
  });

  els.search.addEventListener("click", search);
  els.suggestions.addEventListener("click", event => {
    const suggestion = event.target.closest(".suggestion[data-index]");
    if (suggestion) activateSuggestion(Number(suggestion.dataset.index));
    const copyMissing = event.target.closest("[data-copy-missing]");
    if (copyMissing)
      copyText(`Missing PAIRS query: ${copyMissing.dataset.copyMissing}`, copyMissing);
  });

  document.addEventListener("click", event => {
    if (!event.target.closest(".search-wrap")) closeSuggestions();
    const copyMissing = event.target.closest("[data-copy-missing]");
    if (copyMissing)
      copyText(`Missing PAIRS query: ${copyMissing.dataset.copyMissing}`, copyMissing);
  });

  els.textMode.addEventListener("click", () => setSearchMode("text"));
  els.sequenceMode.addEventListener("click", () => setSearchMode("sequence"));
  els.sequenceSearch.addEventListener("click", () =>
    runSequenceSearch(els.heavySequence.value, els.lightSequence.value),
  );
  els.batchSearch.addEventListener("click", runBatchSearch);
  els.batchClear.addEventListener("click", () => {
    els.batchQuery.value = "";
    state.batchRows = [];
    els.batchResults.innerHTML = "";
    els.batchProgress.textContent = "";
    els.batchCsv.disabled = true;
  });
  els.batchCsv.addEventListener("click", downloadBatchCsv);

  $$(".chip").forEach(chip =>
    chip.addEventListener("click", () => {
      setSearchMode("text");
      els.q.value = chip.dataset.q;
      search();
    }),
  );

  els.filters.addEventListener("click", event => {
    const filter = event.target.closest(".filter");
    if (filter && !filter.hidden) handleFilter(filter.dataset.filter);
  });
  els.sort.addEventListener("change", handleSort);
  els.loadMore.addEventListener("click", handleLoadMore);

  els.results.addEventListener("click", async event => {
    const expand = event.target.closest(".expand");
    if (expand) {
      const card = expand.closest(".card");
      card.classList.toggle("open");
      expand.setAttribute("aria-expanded", card.classList.contains("open") ? "true" : "false");
      expand.textContent = card.classList.contains("open") ? "−" : "+";
      if (card.classList.contains("open")) await loadCardDetails(card);
      return;
    }

    const sequenceCopy = event.target.closest("[data-copy-seq]");
    if (sequenceCopy) {
      copyText(sequenceCopy.dataset.copySeq, sequenceCopy);
      return;
    }

    const antibodyLink = event.target.closest("[data-ab-link]");
    if (antibodyLink) {
      event.preventDefault();
      openStandaloneAntibody(antibodyLink.dataset.abLink, true);
      return;
    }

    const targetLink = event.target.closest("[data-target-link]");
    if (targetLink) {
      event.preventDefault();
      const target = state.targets.find(item => item.id === targetLink.dataset.targetLink);
      if (target) selectTarget(target, true);
      return;
    }

    const copyAntibodyUrl = event.target.closest("[data-copy-ab-url]");
    if (copyAntibodyUrl) {
      const url = new URL(buildViewUrl({ ab: copyAntibodyUrl.dataset.copyAbUrl }), location.origin);
      copyText(url.href, copyAntibodyUrl);
    }
  });

  els.copyUrl.addEventListener("click", () => copyText(location.href, els.copyUrl));
  els.csv.addEventListener("click", exportCSV);
  els.fasta.addEventListener("click", exportFASTA);

  els.browseQuery.addEventListener("input", () => {
    state.browseShown = 60;
    renderBrowse();
  });
  els.browseSort.addEventListener("change", () => {
    state.browseShown = 60;
    renderBrowse();
  });
  els.browseFacets.addEventListener("change", () => {
    state.browseShown = 60;
    renderBrowse();
  });
  els.browseMore.addEventListener("click", () => {
    state.browseShown += 60;
    renderBrowse();
  });
  els.browseGrid?.addEventListener("click", () => setBrowseView("grid"));
  els.browseField?.addEventListener("click", () => setBrowseView("field"));
  els.targetGrid.addEventListener("click", event => {
    const card = event.target.closest("[data-browse-target]");
    if (!card) return;
    const target = state.targets.find(item => item.id === card.dataset.browseTarget);
    if (target) selectTarget(target, true);
  });
  els.bindingFieldSvg?.addEventListener("click", event => {
    const node = event.target.closest?.("[data-field-target]");
    if (!node) return;
    const target = state.targets.find(item => item.id === node.dataset.fieldTarget);
    if (target) selectTarget(target, true);
  });
  els.bindingFieldSvg?.addEventListener("keydown", event => {
    if (event.key !== "Enter" && event.key !== " ") return;
    const node = event.target.closest?.("[data-field-target]");
    if (!node) return;
    event.preventDefault();
    const target = state.targets.find(item => item.id === node.dataset.fieldTarget);
    if (target) selectTarget(target, true);
  });
  els.bindingFieldSvg?.addEventListener("focusin", event => {
    const node = event.target.closest?.("[data-field-target]");
    if (node) showFieldTooltip(node);
  });
  els.bindingFieldSvg?.addEventListener("mouseover", event => {
    const node = event.target.closest?.("[data-field-target]");
    if (node) showFieldTooltip(node);
  });
  els.bindingFieldSvg?.addEventListener("focusout", hideFieldTooltip);
  els.bindingFieldSvg?.addEventListener("mouseleave", hideFieldTooltip);

  $("#sourceStatusBtn").addEventListener("click", () => els.modal.classList.add("open"));
  $("#closeModal").addEventListener("click", () => els.modal.classList.remove("open"));
  els.modal.addEventListener("click", event => {
    if (event.target === els.modal) els.modal.classList.remove("open");
  });

  window.addEventListener("popstate", async () => {
    const parameters = new URLSearchParams(location.search);
    const antibodyId = parameters.get("ab");
    const targetId = parameters.get("target");
    if (antibodyId) await openStandaloneAntibody(antibodyId, false);
    else if (targetId) {
      const target = state.targets.find(item => item.id === targetId);
      if (target) await selectTarget(target, false);
    }
  });

  requestAnimationFrame(() => {
    document.body.classList.add("page-ready");
    setSearchMode("text");
  });
  init();
})();
