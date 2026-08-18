(() => {
  "use strict";

  const SUPPORTED_SCHEMA = 4;
  const DATA_ROOT = `data/v${SUPPORTED_SCHEMA}`;
  const PAGE_RENDER_SIZE = 30;
  const $ = selector => document.querySelector(selector);
  const $$ = selector => [...document.querySelectorAll(selector)];

  const els = {
    q: $("#query"),
    headerSearchForm: $("#headerSearchForm"),
    headerSearchInput: $("#headerSearchInput"),
    search: $("#searchBtn"),
    suggestions: $("#suggestions"),
    textMode: $("#textModeBtn"),
    sequenceMode: $("#sequenceModeBtn"),
    textPanel: $("#textSearchPanel"),
    sequencePanel: $("#sequenceSearchPanel"),
    sequenceIntent: $("#sequenceIntent"),
    sequenceChain: $("#sequenceChain"),
    heavySequence: $("#heavySequenceQuery"),
    lightSequence: $("#lightSequenceQuery"),
    sequenceSearch: $("#sequenceSearchBtn"),
    sequenceType: $("#sequenceType"),
    nearMatches: $("#nearMatches"),
    nearMatchesControl: $("#nearMatchesControl"),
    batchQuery: $("#batchSequenceQuery"),
    batchSearch: $("#batchSearchBtn"),
    batchClear: $("#batchClearBtn"),
    batchCsv: $("#batchCsvBtn"),
    batchProgress: $("#batchProgress"),
    batchResults: $("#batchResults"),
    main: $("#main"),
    targetName: $("#targetName"),
    targetMeta: $("#targetMeta"),
    entitySummary: $("#entitySummary"),
    results: $("#results"),
    summary: $("#summaryGrid"),
    loadMore: $("#loadMore"),
    sort: $("#sortSelect"),
    status: $("#statusStrip"),
    warning: $("#sourceWarning"),
    footer: $("#footerSnapshot"),
    modal: $("#sourceModal"),
    sourceStatusBtn: $("#sourceStatusBtn"),
    closeModal: $("#closeModal"),
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
    advancedFilters: $("#advancedFilters"),
    resultSourceFilter: $("#resultSourceFilter"),
    resultSpeciesFilter: $("#resultSpeciesFilter"),
    resultSequenceFilter: $("#resultSequenceFilter"),
    resultEvidenceFilter: $("#resultEvidenceFilter"),
    clearAdvancedFilters: $("#clearAdvancedFiltersBtn"),
    browseFacets: $("#browseFacets"),
    copyUrl: $("#copyUrlBtn"),
    csv: $("#csvBtn"),
    fasta: $("#fastaBtn"),
    selectionBar: $("#selectionBar"),
    selectionCount: $("#selectionCount"),
    selectVisible: $("#selectVisibleBtn"),
    selectFiltered: $("#selectFilteredBtn"),
    clearSelection: $("#clearSelectionBtn"),
    compareSelected: $("#compareSelectedBtn"),
    datasetMode: $("#datasetMode"),
    datasetPreview: $("#datasetPreview"),
    selectionExportType: $("#selectionExportType"),
    codingPreset: $("#codingPreset"),
    codingPresetControl: $("#codingPresetControl"),
    selectionExportSheet: $("#selectionExportSheet"),
    confirmExportSelected: $("#confirmExportSelectedBtn"),
    exportSelected: $("#exportSelectedBtn"),
    comparisonPanel: $("#comparisonPanel"),
    comparisonContent: $("#comparisonContent"),
    closeComparison: $("#closeComparisonBtn"),
    recentSearches: $("#recentSearches"),
    exportWorkspace: $("#exportWorkspaceBtn"),
    importWorkspace: $("#importWorkspaceInput"),
    workspaceStatus: $("#workspaceStatus"),
    hero: $(".hero"),
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
    similarityBucketCache: new Map(),
    clusterCache: new Map(),
    targetPageCache: new Map(),
    loadedTargetPages: new Set(),
    loadedFunctionalPages: new Set(),
    loadedNegativePages: new Set(),
    browseShown: 60,
    browseView: "grid",
    fieldObserver: null,
    heroObserver: null,
    headerObserver: null,
    sequenceQueryLabel: "",
    sequenceAlignmentQuery: null,
    batchRows: [],
    selectedIds: new Set(),
    datasetPreviewToken: 0,
    recentViews: [],
    pendingWorkspaceFilter: "all",
    includeDescendants: false,
    descendantRowsLoaded: new Set(),
    modalOpener: null,
  };

  const SELECTION_STORAGE_KEY = "pairs-selection-v1";
  const WORKSPACE_STORAGE_KEY = "pairs-workspace-v1";

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
    return targetNames(antibody.direct_targets || []);
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

  function isExplicitVhh(antibody) {
    const descriptor =
      `${antibody.format || ""} ${antibody.metadata?.domain_type || ""}`.toLowerCase();
    return Boolean(
      antibody.has_heavy !== false &&
      (antibody.heavy || antibody.has_heavy) &&
      !(antibody.light || antibody.has_light) &&
      /(?:\bvhh\b|\bnb\b|nanob(?:ody|odies)|single[- ]domain|single[- ]chain\s+vhh)/.test(
        descriptor,
      ),
    );
  }

  function derivedSequenceQuality(antibody) {
    if (antibody.sequence_quality) return antibody.sequence_quality;
    const heavy = antibody.heavy || "";
    const light = antibody.light || "";
    const hasHeavy = Boolean(heavy || antibody.has_heavy);
    const hasLight = Boolean(light || antibody.has_light);
    return {
      pairing:
        hasHeavy && hasLight
          ? "paired"
          : hasHeavy
            ? "heavy_only"
            : hasLight
              ? "light_only"
              : "no_full_chain",
      heavy_length: heavy.length || antibody.heavy_length || 0,
      light_length: light.length || antibody.light_length || 0,
      explicit_vhh: isExplicitVhh(antibody),
      ambiguous_residues: [...new Set(`${heavy}${light}`)]
        .filter(residue => "*BJOUXZ".includes(residue))
        .sort(),
      source_format_quarantined: Boolean(antibody.metadata?.sequence_quarantine),
      completeness: antibody.metadata?.sequence_completeness || "unknown_not_inferred",
      source_nucleotide_available: Boolean(
        antibody.vh_nt_source ||
        antibody.vl_nt_source ||
        antibody.heavy_nt_source ||
        antibody.light_nt_source,
      ),
    };
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
      const targetSummary = directTargets.length
        ? `Direct target · ${directTargets.slice(0, 4).join(" · ")}`
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

  function setSearchMode(mode, options = {}) {
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
    els.textMode?.setAttribute("tabindex", sequence ? "-1" : "0");
    els.sequenceMode?.setAttribute("tabindex", sequence ? "0" : "-1");
    if (activeTab && els.textMode?.parentElement) {
      els.textMode.parentElement.style.setProperty("--tab-left", `${activeTab.offsetLeft}px`);
      els.textMode.parentElement.style.setProperty("--tab-width", `${activeTab.offsetWidth}px`);
    }
    closeSuggestions();
    if (options.focusInput !== false) {
      if (sequence) els.heavySequence?.focus();
      else els.q?.focus();
    }
  }

  function syncSearchInputs(value) {
    const query = String(value ?? "");
    if (els.q && els.q.value !== query) els.q.value = query;
    if (els.headerSearchInput && els.headerSearchInput.value !== query)
      els.headerSearchInput.value = query;
  }

  function setupHeaderSearch() {
    if (!els.headerSearchForm) return;
    syncSearchInputs(els.q?.value || "");
    const setVisible = visible => {
      els.headerSearchForm.hidden = !visible;
      els.headerSearchForm.setAttribute("aria-hidden", String(!visible));
    };
    setVisible(false);
    if (!els.hero || !("IntersectionObserver" in window)) {
      setVisible(true);
      return;
    }
    state.headerObserver = new IntersectionObserver(
      entries => entries.forEach(entry => setVisible(!entry.isIntersecting)),
      { threshold: 0.05 },
    );
    state.headerObserver.observe(els.hero);
  }

  function setExpandButtonState(button, expanded) {
    if (!button) return;
    button.setAttribute("aria-expanded", String(expanded));
    button.setAttribute("aria-label", expanded ? "Hide antibody details" : "Show antibody details");
    button.innerHTML = `${expanded ? "Hide details" : "Details"} <span class="expand-chevron" aria-hidden="true">${expanded ? "⌃" : "⌄"}</span>`;
  }

  function syncSequenceType() {
    if (!els.sequenceIntent || !els.sequenceChain || !els.sequenceType) return;
    const similarity = els.sequenceIntent.value === "similarity";
    const similarityChains = ["heavy", "light", "paired"];
    if (similarity && !similarityChains.includes(els.sequenceChain.value)) {
      els.sequenceChain.value = "heavy";
    }
    els.sequenceType.value = similarity
      ? `${els.sequenceChain.value}_similarity`
      : els.sequenceChain.value;
    if (els.nearMatchesControl) {
      els.nearMatchesControl.hidden =
        similarity || !["auto", "cdrh3", "cdrl3"].includes(els.sequenceChain.value);
    }
    for (const option of els.sequenceChain.options) {
      option.disabled = similarity && !similarityChains.includes(option.value);
    }
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
      sequence_quality: derivedSequenceQuality(antibody),
      shard: antibodyShardFromId(antibody.id),
    };
  }

  function fnv1a32(value) {
    let hash = 2166136261;
    for (const character of value) hash = Math.imul(hash ^ character.charCodeAt(0), 16777619) >>> 0;
    return hash;
  }

  function sequenceSignature(sequence, k = 5, size = 32) {
    const hashes = new Set();
    for (let index = 0; index <= sequence.length - k; index += 1) {
      const peptide = sequence.slice(index, index + k);
      if (/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(peptide)) hashes.add(fnv1a32(peptide));
    }
    return [...hashes]
      .sort((left, right) => left - right)
      .slice(0, size)
      .map(value => value.toString(16).padStart(8, "0"));
  }

  async function similarityCandidates(sequence, field, limit = 140) {
    const signature = sequenceSignature(sequence);
    if (!signature.length) return [];
    const buckets = [...new Set(signature.map(value => value[0]))];
    await Promise.all(
      buckets.map(async bucket => {
        const key = `${field}:${bucket}`;
        if (state.similarityBucketCache.has(key)) return;
        const payload = await getJSON(`${DATA_ROOT}/similarity/${field}/${bucket}.json`);
        state.similarityBucketCache.set(key, payload);
      }),
    );
    const counts = new Map();
    for (const signatureHash of signature) {
      const payload = state.similarityBucketCache.get(`${field}:${signatureHash[0]}`) || {};
      for (const antibodyId of payload[signatureHash] || []) {
        counts.set(antibodyId, (counts.get(antibodyId) || 0) + 1);
      }
    }
    return [...counts]
      .map(([id, sharedSignatureHashes]) => ({ id, sharedSignatureHashes }))
      .sort((left, right) => right.sharedSignatureHashes - left.sharedSignatureHashes)
      .slice(0, limit);
  }

  async function alignedSimilarityCandidates(sequence, field) {
    const candidates = await similarityCandidates(sequence, field);
    const exact = await lookupSequence(sequence);
    const candidateById = new Map(candidates.map(candidate => [candidate.id, candidate]));
    for (const match of exact.filter(item => item.field === field)) {
      candidateById.set(match.id, { id: match.id, sharedSignatureHashes: 32 });
    }
    const records = await fetchFullRecords([...candidateById.keys()]);
    const aligned = [];
    for (const candidate of candidateById.values()) {
      const antibody = records.get(candidate.id);
      const targetSequence = antibody?.[field] || "";
      if (!targetSequence) continue;
      const lengthRatio =
        Math.min(sequence.length, targetSequence.length) /
        Math.max(sequence.length, targetSequence.length);
      if (lengthRatio < 0.65) continue;
      const alignment = alignProteins(sequence, targetSequence);
      if (alignment.identity < 45 || alignment.coverage < 70) continue;
      aligned.push({
        id: candidate.id,
        antibody,
        field,
        shared_signature_hashes: candidate.sharedSignatureHashes,
        identity: alignment.identity,
        coverage: alignment.coverage,
        score: alignment.score,
        aligned_query: alignment.reference,
        aligned_target: alignment.query,
      });
    }
    return aligned.sort(
      (left, right) =>
        right.identity * right.coverage - left.identity * left.coverage ||
        right.shared_signature_hashes - left.shared_signature_hashes ||
        left.id.localeCompare(right.id),
    );
  }

  async function runFullChainSimilarity(first, second, type, excludeId = "") {
    if (first.length < 40 || (type === "paired_similarity" && second.length < 40)) {
      throw new Error(
        "Full-chain similarity requires at least 40 amino acids per requested chain.",
      );
    }
    const paired = type === "paired_similarity";
    const firstField = type === "light_similarity" ? "light" : "heavy";
    const [firstMatches, secondMatches] = await Promise.all([
      alignedSimilarityCandidates(first, firstField),
      paired ? alignedSimilarityCandidates(second, "light") : Promise.resolve([]),
    ]);
    const secondById = new Map(secondMatches.map(match => [match.id, match]));
    const results = firstMatches
      .filter(match => (!paired || secondById.has(match.id)) && match.id !== excludeId)
      .map(match => {
        const secondMatch = secondById.get(match.id);
        const metrics = secondMatch ? [match, secondMatch] : [match];
        const similarityScore =
          metrics.reduce((sum, item) => sum + (item.identity * item.coverage) / 100, 0) /
          metrics.length;
        return {
          antibody: summaryFromFull(match.antibody),
          relationships: [],
          evidence: [],
          sources: match.antibody.sources || [],
          interactions: [],
          match_fields: [],
          similarity: Object.fromEntries(metrics.map(item => [item.field, item])),
          similarity_score: similarityScore,
        };
      });
    state.mode = "sequence";
    state.selected = null;
    state.rawResults = results;
    state.sequenceAlignmentQuery = {
      heavy: firstField === "heavy" ? first : "",
      light: firstField === "light" ? first : paired ? second : "",
    };
    state.sequenceQueryLabel = `${paired ? "Paired VH + VL" : firstField === "heavy" ? "VH / VHH" : "VL"} aligned similarity · indexed candidate retrieval`;
    state.filter = "all";
    state.shown = PAGE_RENDER_SIZE;
    resetFilterButtons();
    updateFilterAvailability();
    els.targetName.textContent = excludeId
      ? "Related public sequences"
      : "Full-chain similarity candidates";
    els.targetMeta.textContent = `${fmt(results.length)} aligned public records · minimum 45% identity and 70% query coverage · similarity does not establish specificity`;
    els.main.classList.add("active");
    apply();
    if (!results.length) renderSequenceEmpty(first.length, paired);
    history.pushState({}, "", location.pathname);
    scrollToResults();
  }

  async function findRelatedSequences(antibodyId, button) {
    await withBusyButton(button, "Finding…", async () => {
      const antibody = await fetchAntibody(antibodyId);
      if (!antibody?.heavy && !antibody?.light) return;
      setSearchMode("sequence");
      els.heavySequence.value = antibody.heavy || antibody.light;
      els.lightSequence.value = antibody.heavy && antibody.light ? antibody.light : "";
      const type =
        antibody.heavy && antibody.light
          ? "paired_similarity"
          : antibody.heavy
            ? "heavy_similarity"
            : "light_similarity";
      if (els.sequenceIntent && els.sequenceChain) {
        els.sequenceIntent.value = "similarity";
        els.sequenceChain.value = type.replace("_similarity", "");
        syncSequenceType();
      } else {
        els.sequenceType.value = type;
      }
      await runFullChainSimilarity(
        antibody.heavy || antibody.light,
        antibody.heavy && antibody.light ? antibody.light : "",
        type,
        antibodyId,
      );
    });
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
      if (type.endsWith("_similarity")) {
        await runFullChainSimilarity(first, second, type);
        return;
      }
      if (["cdrh3", "cdrl3"].includes(type)) {
        const hits = await lookupCdr(first, type, els.nearMatches.checked);
        const records = await fetchFullRecords([...new Set(hits.flatMap(hit => hit.antibody_ids))]);
        state.mode = "sequence";
        state.selected = null;
        state.sequenceAlignmentQuery = null;
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
      state.sequenceAlignmentQuery = null;
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
    if (state.sequenceAlignmentQuery) {
      els.results.innerHTML = `<div class="empty state-panel">${stateGlyph("empty")}<strong>No aligned similarity candidates</strong><span>The indexed candidate step found no public ${paired ? "paired " : ""}record meeting 45% identity and 70% query coverage. This is not evidence that no related sequence exists outside the PAIRS snapshot.</span><div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a></div></div>`;
      els.loadMore.hidden = true;
      return;
    }
    const complete = isCompleteSnapshot();
    els.results.innerHTML = `<div class="empty state-panel">${stateGlyph("empty")}<strong>${complete ? "No exact public sequence match" : "No exact match in this partial snapshot"}</strong><span>${complete ? "PAIRS found no" : "PAIRS cannot conclude absence from an incomplete dataset; it found no"} exact ${paired ? "paired " : ""}match for the normalized ${length}-aa query.</span><div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a></div></div>`;
    els.loadMore.hidden = true;
  }

  function showSearchError(message) {
    state.mode = "sequence";
    state.sequenceAlignmentQuery = null;
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
    state.sequenceAlignmentQuery = null;
    state.targetPageCache.clear();
    state.loadedTargetPages.clear();
    state.loadedFunctionalPages.clear();
    state.loadedNegativePages.clear();
    state.descendantRowsLoaded.clear();
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
      mergeRowsIntoState(rows);
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
      mergeRowsIntoState(pageRows);
    });
  }

  function mergeRowsIntoState(rows) {
    const byId = new Map(state.rawResults.map(row => [row.antibody.id, row]));
    for (const row of rows) {
      const existing = byId.get(row.antibody.id);
      if (!existing) {
        state.rawResults.push(row);
        byId.set(row.antibody.id, row);
        continue;
      }
      for (const field of ["relationships", "evidence", "sources", "match_fields"]) {
        existing[field] = [...new Set([...(existing[field] || []), ...(row[field] || [])])];
      }
      const interactions = new Map(
        (existing.interactions || []).map(interaction => [interaction.id, interaction]),
      );
      for (const interaction of row.interactions || [])
        interactions.set(interaction.id, interaction);
      existing.interactions = [...interactions.values()];
    }
  }

  function descendantTargets(target) {
    const output = [];
    const pending = [...(target.children || [])];
    const seen = new Set();
    while (pending.length) {
      const id = pending.shift();
      if (seen.has(id)) continue;
      seen.add(id);
      const child = state.targets.find(item => item.id === id);
      if (!child) continue;
      output.push(child);
      pending.push(...(child.children || []));
    }
    return output;
  }

  async function loadDescendantRows(category = activeTargetPageCategory()) {
    if (!state.selected || !state.includeDescendants) return;
    const descendants = descendantTargets(state.selected);
    const requests = [];
    for (const child of descendants) {
      const pageCount =
        category === "negative"
          ? child.negative_page_count || 0
          : category === "functional"
            ? child.functional_page_count || 0
            : child.page_count || 0;
      for (let page = 1; page <= pageCount; page += 1) {
        const key = `${child.id}:${category}:${page}`;
        if (state.descendantRowsLoaded.has(key)) continue;
        state.descendantRowsLoaded.add(key);
        requests.push(
          getJSON(
            category === "negative"
              ? negativeTargetPageUrl(child, page)
              : category === "functional"
                ? functionalTargetPageUrl(child, page)
                : targetPageUrl(child, page),
          ),
        );
      }
    }
    const pages = await Promise.all(requests);
    pages.forEach(mergeRowsIntoState);
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
    state.includeDescendants =
      !push && new URL(location.href).searchParams.get("scope") === "descendants";
    resetResultState();
    updateFilterAvailability();
    closeSuggestions();
    syncSearchInputs(target.name);
    els.main.classList.add("active");
    els.targetName.textContent = target.name;
    els.main.dataset.view = "target";
    document.title = `${target.name} antibodies — PAIRS`;
    renderTargetEntitySummary(target);
    updateTargetMeta();
    els.summary.innerHTML = "";
    els.results.innerHTML = loadingState("Loading target results…");
    els.loadMore.hidden = true;
    try {
      await loadTargetPage(1, "positive");
      await loadDescendantRows("positive");
      apply();
      rememberView("target", target.id, target.name);
      if (push) {
        const url = new URL(location.href);
        url.search = "";
        url.searchParams.set("target", target.id);
        if (state.includeDescendants) url.searchParams.set("scope", "descendants");
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
    let countLabel =
      positiveCount == null
        ? `${fmt(state.selected.result_count)} indexed antibodies (negative-only rows excluded)`
        : `${fmt(positiveCount)} antibodies with positive evidence`;
    if (hasAdvancedFilters() || state.includeDescendants)
      countLabel = `${fmt(state.filtered.length)} antibodies match active filters`;
    const scopeLabel = state.includeDescendants
      ? `including ${descendantTargets(state.selected).length} descendant targets`
      : "exact target only";
    els.targetMeta.textContent = `${countLabel} · ${scopeLabel} · ${fmt(state.selected.count)} source-level relationships · ${state.selected.sources.join(" · ")}`;
  }

  function renderTargetEntitySummary(target) {
    const entity = target.entity || {};
    const identifiers = entity.identifiers || {};
    const identifierLinks = [
      identifiers.uniprot &&
        `<a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(identifiers.uniprot)}/entry" target="_blank" rel="noopener">UniProt ${esc(identifiers.uniprot)}</a>`,
      identifiers.hgnc &&
        `<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${encodeURIComponent(identifiers.hgnc)}" target="_blank" rel="noopener">${esc(identifiers.hgnc)}</a>`,
      (identifiers.ncbi_taxonomy || entity.organism?.ncbi_taxonomy_id) &&
        `<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${encodeURIComponent(identifiers.ncbi_taxonomy || entity.organism.ncbi_taxonomy_id)}" target="_blank" rel="noopener">NCBI Taxonomy ${esc(identifiers.ncbi_taxonomy || entity.organism.ncbi_taxonomy_id)}</a>`,
    ].filter(Boolean);
    const parent = target.hierarchy?.parent_id
      ? state.targets.find(item => item.id === target.hierarchy.parent_id)
      : null;
    const children = (target.children || [])
      .map(id => state.targets.find(item => item.id === id))
      .filter(Boolean);
    els.entitySummary.hidden = false;
    els.entitySummary.innerHTML = `<div><span class="eyebrow">${esc(entity.kind || "locally resolved target")}</span><strong>${esc(target.name)}</strong><span>${esc(entity.organism?.name || "Organism not authoritatively mapped")}</span><span>${identifierLinks.join(" · ") || "No authoritative external accession assigned"}</span>${parent ? `<span>Parent: <a href="${esc(buildViewUrl({ target: parent.id }))}" data-target-link="${esc(parent.id)}">${esc(parent.name)}</a> · ${esc(target.hierarchy.relation || "related_to")}</span>` : ""}</div><div>${children.length ? `<span>Child targets: ${children.map(child => `<a href="${esc(buildViewUrl({ target: child.id }))}" data-target-link="${esc(child.id)}">${esc(child.name)}</a>`).join(" · ")}</span><label class="scope-toggle"><input type="checkbox" data-include-descendants${state.includeDescendants ? " checked" : ""}> Include subdomains/children</label>` : '<span class="meta">Exact target scope</span>'}</div>`;
  }

  async function setDescendantScope(enabled) {
    if (!state.selected) return;
    state.includeDescendants = enabled;
    resetResultState();
    renderTargetEntitySummary(state.selected);
    els.results.innerHTML = loadingState("Loading scoped target results…");
    await loadTargetPage(1, "positive");
    await loadDescendantRows("positive");
    const url = new URL(location.href);
    if (enabled) url.searchParams.set("scope", "descendants");
    else url.searchParams.delete("scope");
    history.replaceState({}, "", url);
    apply();
  }

  function renderStaleEntity(kind, id) {
    state.selected = null;
    state.rawResults = [];
    state.filtered = [];
    els.main.classList.add("active");
    els.main.dataset.view = kind;
    els.targetName.textContent = `${kind === "target" ? "Target" : "Antibody"} not present`;
    els.targetMeta.textContent = `The requested ID is not present in snapshot ${state.manifest?.snapshot || "unknown"}.`;
    els.entitySummary.hidden = true;
    els.results.innerHTML = `<div class="empty state-panel"><strong>Stale or unknown PAIRS ID</strong><span>${esc(id)} was not found in this snapshot. PAIRS has not substituted another entity.</span><div class="empty-actions"><a class="secondary" href="#browse">Browse current targets</a></div></div>`;
    els.loadMore.hidden = true;
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
      state.sequenceAlignmentQuery = null;
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
      els.main.dataset.view = "antibody";
      document.title = `${antibody.name || antibody.id} — PAIRS`;
      els.entitySummary.hidden = false;
      els.entitySummary.innerHTML = `<div><span class="eyebrow">Sequence entity</span><strong>${esc(antibody.name || antibody.id)}</strong><span>${esc(antibody.id)} · ${(antibody.sources || []).length} upstream source${(antibody.sources || []).length === 1 ? "" : "s"}</span></div><div><span>This page represents a PAIRS sequence entity. Source records and multispecific construct arms remain separate below.</span></div>`;
      els.targetMeta.textContent = `${fmt(antibody.target_count)} target annotations · ${(antibody.sources || []).join(" · ")} · ${antibody.heavy && antibody.light ? "paired VH + VL" : "sequence record"}`;
      apply();
      rememberView("antibody", antibody.id, antibody.name || antibody.id);
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

  function reportedOrigin(antibody) {
    const value = `${antibody.organism || ""} ${antibody.metadata?.genetics || ""}`.toLowerCase();
    if (!value.trim()) return "unknown";
    if (/human|humanized|homo sapiens/.test(value)) return "human";
    if (/mouse|murine|mus musculus/.test(value)) return "mouse";
    return "other";
  }

  function passesAdvancedFilters(result) {
    const antibody = result.antibody;
    const source = els.resultSourceFilter.value;
    const species = els.resultSpeciesFilter.value;
    const sequenceFilter = els.resultSequenceFilter.value;
    const evidence = els.resultEvidenceFilter.value;
    if (source && !result.sources.includes(source)) return false;
    if (species && reportedOrigin(antibody) !== species) return false;
    const quality = derivedSequenceQuality(antibody);
    if (sequenceFilter === "paired" && quality.pairing !== "paired") return false;
    if (sequenceFilter === "vhh" && !quality.explicit_vhh) return false;
    if (sequenceFilter === "heavy_only" && quality.pairing !== "heavy_only") return false;
    if (sequenceFilter === "light_only" && quality.pairing !== "light_only") return false;
    if (sequenceFilter === "cdr_only" && quality.pairing !== "no_full_chain") return false;
    if (evidence && !result.evidence.includes(evidence)) return false;
    return true;
  }

  function hasAdvancedFilters() {
    return Boolean(
      els.resultSourceFilter.value ||
      els.resultSpeciesFilter.value ||
      els.resultSequenceFilter.value ||
      els.resultEvidenceFilter.value,
    );
  }

  function passes(result) {
    const antibody = result.antibody;
    if (!passesAdvancedFilters(result)) return false;
    if (state.filter === "negative") return result.relationships.some(isNegative);
    if (state.filter === "functional") return result.relationships.some(isFunctionalPositive);
    if (state.mode === "target" && !result.relationships.some(isPrimaryPositive)) return false;
    if (state.filter === "paired") return Boolean(antibody.has_heavy && antibody.has_light);
    if (state.filter === "therapeutic") return Boolean(antibody.therapeutic_status);
    if (state.filter === "structure") return hasExactStructure(antibody);
    return true;
  }

  function rank(result) {
    if (result.similarity_score != null) return result.similarity_score * 10;
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
    if (
      state.mode === "target" &&
      state.selected?.stats &&
      !hasAdvancedFilters() &&
      !state.includeDescendants
    ) {
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
    const quality = derivedSequenceQuality(result.antibody);
    if (result.antibody.has_heavy && result.antibody.has_light) {
      output.push('<span class="badge paired">VH + VL</span>');
    } else if (quality.explicit_vhh) {
      output.push('<span class="badge paired">VHH · source-labelled</span>');
    } else if (quality.pairing === "heavy_only") {
      output.push('<span class="badge">VH only</span>');
    } else if (quality.pairing === "light_only") {
      output.push('<span class="badge">VL only</span>');
    } else {
      output.push('<span class="badge">CDR / metadata only</span>');
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
    for (const [field, similarity] of Object.entries(result.similarity || {})) {
      output.push(
        `<span class="badge functional">${esc(field.toUpperCase())} · ${similarity.identity.toFixed(1)}% ID</span>`,
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
        const measurements = (interaction.measurements || [])
          .map(item =>
            [
              item.metric,
              `${item.qualifier || ""}${(item.value ?? item.raw_value) || ""}`,
              item.unit,
            ]
              .filter(Boolean)
              .join(" "),
          )
          .join("; ");
        const details = [
          interaction.evidence,
          interaction.source_record_id &&
            `${provenance.scope === "record" ? "Exact source record" : "Source record ID"}: ${interaction.source_record_id}`,
          interaction.assertion_origin === "derived_hierarchy" && "PAIRS-derived parent target",
          interaction.assertion_origin === "source_epitope" && "Derived from source epitope field",
          interaction.epitope && `Epitope: ${interaction.epitope}`,
          interaction.assay && `Assay: ${interaction.assay}`,
          measurements && `Measurements: ${measurements}`,
          interaction.reference,
          interaction.note,
        ]
          .filter(Boolean)
          .join(" · ");
        return `<div class="evidence-row ${isNegative(interaction.relationship) ? "negative" : ""}"><strong>${esc(relationLabel(interaction.relationship))} · ${provenance.html}</strong><span>${esc(details)}</span></div>`;
      })
      .join("");
  }

  function evidenceLine(result) {
    if (state.mode === "sequence") {
      if (result.similarity) {
        const metrics = Object.entries(result.similarity)
          .map(
            ([field, similarity]) =>
              `${field.toUpperCase()} ${similarity.identity.toFixed(1)}% identity · ${similarity.coverage.toFixed(1)}% coverage`,
          )
          .join(" · ");
        return `<button class="evidence-trigger" data-open-evidence type="button">Aligned similarity candidate · ${esc(metrics)}</button>`;
      }
      const fields = (result.match_fields || []).map(field => field.toUpperCase()).join(" + ");
      return `<span class="evidence-trigger static">Exact ${esc(fields || "sequence")} match · PAIRS sequence index</span>`;
    }
    const eligible = (result.interactions || []).filter(interaction => {
      if (state.filter === "negative") return isNegative(interaction.relationship);
      if (state.filter === "functional") return isFunctionalPositive(interaction.relationship);
      return isPrimaryPositive(interaction.relationship);
    });
    const interaction = eligible.sort(
      (left, right) => (evidenceRank[right.evidence] ?? 0) - (evidenceRank[left.evidence] ?? 0),
    )[0];
    if (!interaction)
      return '<span class="evidence-trigger static">Sequence and provenance record</span>';
    const sourceName = state.manifest?.sources?.[interaction.source]?.name || interaction.source;
    const origin = interaction.record_url ? "exact source record" : "source record";
    return `<button class="evidence-trigger" data-open-evidence type="button">${esc(relationLabel(interaction.relationship))} · ${esc(sourceName)} · ${esc(origin)}</button>`;
  }

  function saveSelection() {
    try {
      localStorage.setItem(SELECTION_STORAGE_KEY, JSON.stringify([...state.selectedIds]));
    } catch {
      // Selection remains available for this session when storage is unavailable.
    }
    saveWorkspaceLocal();
  }

  function restoreSelection() {
    try {
      const ids = JSON.parse(localStorage.getItem(SELECTION_STORAGE_KEY) || "[]");
      state.selectedIds = new Set(ids.filter(id => typeof id === "string" && id.startsWith("ab_")));
    } catch {
      state.selectedIds = new Set();
    }
    updateSelectionBar();
  }

  function workspacePayload() {
    return {
      workspace_schema: "pairs-workspace",
      workspace_version: 1,
      pairs_snapshot: state.manifest?.snapshot || null,
      pairs_schema_version: state.manifest?.schema_version || SUPPORTED_SCHEMA,
      saved_at: new Date().toISOString(),
      selected_ids: [...state.selectedIds],
      recent_views: state.recentViews,
      current_view: location.search || null,
      settings: {
        filter: state.filter,
        sort: els.sort.value,
        dataset_mode: els.datasetMode.value,
        source: els.resultSourceFilter.value,
        species: els.resultSpeciesFilter.value,
        sequence: els.resultSequenceFilter.value,
        evidence: els.resultEvidenceFilter.value,
      },
    };
  }

  function saveWorkspaceLocal() {
    try {
      localStorage.setItem(WORKSPACE_STORAGE_KEY, JSON.stringify(workspacePayload()));
    } catch {
      // The workspace remains in memory when browser storage is unavailable.
    }
  }

  function renderRecentViews() {
    els.recentSearches.innerHTML = '<option value="">Recent views</option>';
    for (const [index, view] of state.recentViews.entries()) {
      els.recentSearches.insertAdjacentHTML(
        "beforeend",
        `<option value="${index}">${esc(view.label)}</option>`,
      );
    }
  }

  function restoreWorkspace() {
    restoreSelection();
    try {
      const payload = JSON.parse(localStorage.getItem(WORKSPACE_STORAGE_KEY) || "null");
      if (payload?.workspace_schema !== "pairs-workspace" || payload.workspace_version !== 1) {
        renderRecentViews();
        return;
      }
      if (Array.isArray(payload.selected_ids)) {
        state.selectedIds = new Set(
          payload.selected_ids.filter(id => /^ab_[0-9a-f]{20}$/.test(id)).slice(0, 5000),
        );
      }
      state.recentViews = Array.isArray(payload.recent_views)
        ? payload.recent_views
            .filter(
              view =>
                view &&
                typeof view.id === "string" &&
                view.id.length <= 300 &&
                ((view.kind === "target" && /^target:/.test(view.id)) ||
                  (view.kind === "antibody" && /^ab_[0-9a-f]{20}$/.test(view.id))) &&
                typeof view.label === "string" &&
                view.label.length <= 200,
            )
            .slice(0, 30)
        : [];
      const settings = payload.settings || {};
      state.pendingWorkspaceFilter = [
        "all",
        "paired",
        "therapeutic",
        "structure",
        "functional",
        "negative",
      ].includes(settings.filter)
        ? settings.filter
        : "all";
      if ([...els.sort.options].some(option => option.value === settings.sort))
        els.sort.value = settings.sort;
      if ([...els.datasetMode.options].some(option => option.value === settings.dataset_mode))
        els.datasetMode.value = settings.dataset_mode;
      for (const [element, value] of [
        [els.resultSourceFilter, settings.source],
        [els.resultSpeciesFilter, settings.species],
        [els.resultSequenceFilter, settings.sequence],
        [els.resultEvidenceFilter, settings.evidence],
      ]) {
        if ([...element.options].some(option => option.value === value)) element.value = value;
      }
      if (payload.pairs_snapshot && payload.pairs_snapshot !== state.manifest?.snapshot) {
        els.workspaceStatus.textContent =
          "Saved workspace uses another snapshot; IDs were retained.";
      }
    } catch {
      state.recentViews = [];
    }
    renderRecentViews();
    updateSelectionBar();
  }

  function rememberView(kind, id, label) {
    state.recentViews = [
      { kind, id, label, viewed_at: new Date().toISOString() },
      ...state.recentViews.filter(view => !(view.kind === kind && view.id === id)),
    ].slice(0, 12);
    renderRecentViews();
    saveWorkspaceLocal();
  }

  function exportWorkspace() {
    downloadBlob(
      "pairs-workspace.json",
      JSON.stringify(workspacePayload(), null, 2),
      "application/json",
    );
  }

  async function importWorkspace(file) {
    if (!file || file.size > 1_000_000) throw new Error("Workspace file is missing or too large.");
    const payload = JSON.parse(await file.text());
    if (
      !payload ||
      payload.workspace_schema !== "pairs-workspace" ||
      payload.workspace_version !== 1 ||
      !Array.isArray(payload.selected_ids)
    ) {
      throw new Error("This is not a supported PAIRS workspace.");
    }
    if (payload.selected_ids.length > 5000) throw new Error("Workspace selection is too large.");
    const invalidId = payload.selected_ids.find(id => !/^ab_[0-9a-f]{20}$/.test(id));
    if (invalidId) throw new Error("Workspace contains an invalid antibody ID.");
    state.selectedIds = new Set(payload.selected_ids);
    state.recentViews = Array.isArray(payload.recent_views)
      ? payload.recent_views
          .filter(
            view =>
              view &&
              typeof view.id === "string" &&
              view.id.length <= 300 &&
              ((view.kind === "target" && /^target:/.test(view.id)) ||
                (view.kind === "antibody" && /^ab_[0-9a-f]{20}$/.test(view.id))) &&
              typeof view.label === "string" &&
              view.label.length <= 200,
          )
          .slice(0, 30)
      : [];
    const settings = payload.settings || {};
    if ([...els.sort.options].some(option => option.value === settings.sort))
      els.sort.value = settings.sort;
    if ([...els.datasetMode.options].some(option => option.value === settings.dataset_mode))
      els.datasetMode.value = settings.dataset_mode;
    for (const [element, value] of [
      [els.resultSourceFilter, settings.source],
      [els.resultSpeciesFilter, settings.species],
      [els.resultSequenceFilter, settings.sequence],
      [els.resultEvidenceFilter, settings.evidence],
    ]) {
      if ([...element.options].some(option => option.value === value)) element.value = value;
    }
    const snapshotMismatch =
      payload.pairs_snapshot && payload.pairs_snapshot !== state.manifest?.snapshot;
    saveSelection();
    renderRecentViews();
    updateSelectionBar();
    const filter = ["all", "paired", "therapeutic", "structure", "functional", "negative"].includes(
      settings.filter,
    )
      ? settings.filter
      : "all";
    await handleFilter(filter);
    saveWorkspaceLocal();
    els.workspaceStatus.textContent = snapshotMismatch
      ? "Imported. Snapshot differs; retained IDs may be stale."
      : "Workspace imported.";
  }

  function updateSelectionBar() {
    const count = state.selectedIds.size;
    if (els.selectionBar) {
      // Keep the export controls out of the reading flow until they have a job.
      els.selectionBar.hidden = count === 0;
    }
    if (!count && els.selectionExportSheet) {
      els.selectionExportSheet.open = false;
      els.exportSelected?.setAttribute("aria-expanded", "false");
    }
    if (els.selectionCount) els.selectionCount.textContent = `${fmt(count)} selected`;
    if (els.clearSelection) els.clearSelection.disabled = count === 0;
    if (els.exportSelected) els.exportSelected.disabled = count === 0;
    if (els.compareSelected) els.compareSelected.disabled = count < 2 || count > 10;
    if (els.compareSelected)
      els.compareSelected.title =
        count > 10 ? "Comparison supports up to 10 antibodies" : "Compare 2–10 antibodies";
    $$("#results [data-select-ab]").forEach(input => {
      input.checked = state.selectedIds.has(input.dataset.selectAb);
    });
    updateDatasetPreview();
  }

  async function updateDatasetPreview() {
    const token = ++state.datasetPreviewToken;
    const count = state.selectedIds.size;
    if (!count) {
      els.datasetPreview.textContent = "";
      return;
    }
    if (els.datasetMode.value === "entities") {
      els.datasetPreview.textContent = `${fmt(count)} entities`;
      return;
    }
    els.datasetPreview.textContent = "Counting…";
    try {
      const { rows, dataset } = await prepareRows("selected");
      if (token !== state.datasetPreviewToken) return;
      els.datasetPreview.textContent = `${fmt(dataset.selected_entity_count)} → ${fmt(rows.length)} ${dataset.mode === "source_records" ? "source records" : "exact pairs"}`;
    } catch {
      if (token === state.datasetPreviewToken) els.datasetPreview.textContent = "Count unavailable";
    }
  }

  function setSelected(antibodyId, selected) {
    if (selected) state.selectedIds.add(antibodyId);
    else state.selectedIds.delete(antibodyId);
    saveSelection();
    updateSelectionBar();
    els.comparisonPanel.hidden = true;
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
        const checked = state.selectedIds.has(antibody.id) ? " checked" : "";
        return `<article class="card" data-ab="${esc(antibody.id)}" data-shard="${esc(antibody.shard)}"><div class="card-main"><label class="result-select" title="Add to export selection"><input type="checkbox" data-select-ab="${esc(antibody.id)}"${checked} /><span class="sr-only">Select ${esc(antibody.name || antibody.id)}</span></label><div><a class="name-link" href="${esc(directUrl)}" data-ab-link="${esc(antibody.id)}">${esc(antibody.name || antibody.id)}</a><div class="meta">${esc(alias || antibody.organism || antibody.format || "Public antibody record")}</div>${evidenceLine(result)}</div><div class="badges">${badges(result)}</div><div class="sources">${esc(result.sources.join(" · ") || "Public source")}<br>${result.evidence.length ? esc(result.evidence.join(" · ")) : esc((result.match_fields || []).length ? "Exact sequence match" : "Sequence/provenance record")}</div><button class="download-one" data-download-ab="${esc(antibody.id)}" type="button" aria-label="Download ${esc(antibody.name || antibody.id)} amino-acid FASTA">FASTA</button><button class="expand" type="button" aria-label="Show antibody details" aria-expanded="false">Details <span class="expand-chevron" aria-hidden="true">⌄</span></button></div><div class="detail"><div class="detail-grid"><div class="sequence-slot"><div class="meta">Sequence details load only when this record is expanded.</div></div><div><div class="section-title">Evidence for this view</div><div class="evidence-list">${evidenceRows(result.interactions)}</div><div class="full-record-slot"></div></div></div></div></article>`;
      })
      .join("");
    $$("#results .card").forEach((card, index) =>
      card.style.setProperty("--card-index", Math.min(index, 11)),
    );
    updateSelectionBar();
  }

  function annotatedSequenceBlock(label, value, cdr3, vGene, jGene) {
    if (!value) return "";
    const start = cdr3 ? value.indexOf(cdr3) : -1;
    const annotated =
      start >= 0
        ? `${coloredSequence(value.slice(0, start))}<mark class="cdr3-region">${coloredSequence(cdr3)}</mark>${coloredSequence(value.slice(start + cdr3.length))}`
        : coloredSequence(value);
    const calls = [vGene && `V: ${vGene}`, jGene && `J: ${jGene}`].filter(Boolean).join(" · ");
    const cdrStatus = cdr3
      ? start >= 0
        ? `Source-supplied CDR3 highlighted · ${fmt(cdr3.length)} aa`
        : "Source CDR3 supplied separately; boundary not found in stored variable region"
      : "CDR3 boundary not supplied by the source";
    return `<div class="seq annotated-seq"><div class="seq-head"><span>${esc(label)} · ${fmt(value.length)} aa${calls ? ` · ${esc(calls)}` : ""}</span><button class="copy" data-copy-seq="${esc(value)}" aria-label="Copy ${esc(label)} sequence">Copy chain</button></div><code class="sequence-residues" aria-label="${esc(value)}">${annotated}</code><div class="annotation-meta"><span>${esc(cdrStatus)}</span>${cdr3 ? `<button class="copy" data-copy-seq="${esc(cdr3)}" aria-label="Copy ${esc(label)} CDR3">Copy CDR3</button>` : ""}</div></div>`;
  }

  function sourceRegionCopies(label, regions) {
    const values = Object.entries(regions || {}).filter(([, value]) => value);
    if (!values.length) return "";
    return `<div class="annotation-meta"><span>${esc(label)} source regions</span>${values
      .map(
        ([region, value]) =>
          `<button class="copy" data-copy-seq="${esc(value)}" aria-label="Copy ${esc(label)} ${esc(region)}">Copy ${esc(region)}</button>`,
      )
      .join("")}</div>`;
  }

  function sequenceQualityHtml(antibody) {
    const quality = derivedSequenceQuality(antibody);
    const pairing = {
      paired: "Paired VH + VL",
      heavy_only: quality.explicit_vhh ? "Source-labelled VHH" : "Heavy chain only",
      light_only: "Light chain only",
      no_full_chain: "No full variable-domain chain",
    }[quality.pairing];
    const ambiguity = quality.ambiguous_residues?.length
      ? `Ambiguous/non-standard residues: ${quality.ambiguous_residues.join(", ")}`
      : "No ambiguous/non-standard residues detected";
    const completeness =
      quality.completeness === "unknown_not_inferred"
        ? "Variable-region completeness not asserted by the source"
        : `Source-reported completeness: ${quality.completeness}`;
    return `<div class="quality-grid"><span>${esc(pairing)}</span><span>VH ${fmt(quality.heavy_length)} aa · VL ${fmt(quality.light_length)} aa</span><span>${esc(ambiguity)}</span><span>${esc(completeness)}</span>${quality.source_format_quarantined ? '<span class="quality-warning">Source sequence formatting quarantined</span>' : ""}</div>`;
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
    const directTargets = antibody.direct_targets || [];
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
          const pdbId = /^[0-9][A-Za-z0-9]{3}$/.test(entry.id) ? entry.id.toUpperCase() : "";
          const external = `<a class="pdb-link" href="https://www.rcsb.org/structure/${encodeURIComponent(entry.id)}" target="_blank" rel="noopener">${esc(entry.id)}${suffix}</a>`;
          return pdbId
            ? `${external}<button class="pdb-link viewer-trigger" type="button" data-view-structure="${esc(pdbId)}" data-structure-tier="${esc(entry.tier)}">3D</button>`
            : external;
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
    const timelineRecords = [...(antibody.source_records || [])].sort((left, right) => {
      const leftTime = left.record_date ? Date.parse(left.record_date) : Number.POSITIVE_INFINITY;
      const rightTime = right.record_date
        ? Date.parse(right.record_date)
        : Number.POSITIVE_INFINITY;
      return (
        (Number.isNaN(leftTime) ? Number.POSITIVE_INFINITY : leftTime) -
          (Number.isNaN(rightTime) ? Number.POSITIVE_INFINITY : rightTime) ||
        String(left.source).localeCompare(String(right.source)) ||
        String(left.record_id).localeCompare(String(right.record_id))
      );
    });
    const timeline = timelineRecords
      .map((record, index) => {
        const sourceName = state.manifest?.sources?.[record.source]?.name || record.source;
        return `<div class="timeline-row"><span>${record.record_date ? esc(record.record_date) : "Date unavailable"}</span><strong>${esc(sourceName)} · ${esc(record.record_id)}</strong><span>${index === 0 && record.record_date ? "Earliest dated record indexed by PAIRS" : esc(record.record_date_field ? `Source field: ${record.record_date_field}` : "No source record date supplied")}</span></div>`;
      })
      .join("");
    const conflicts = (antibody.source_conflicts || [])
      .map(
        conflict =>
          `<div class="evidence-row negative"><strong>Conflicting source annotation · ${esc(conflict.field)}</strong><span>${esc(conflict.incoming_source || "Source")} · ${esc(conflict.incoming_source_record_id || "record ID unavailable")}</span></div>`,
      )
      .join("");
    const conflictHtml = conflicts
      ? `<div class="section-title" style="margin-top:16px">Source conflicts</div><div class="evidence-list">${conflicts}</div>`
      : "";
    const relatedAction =
      antibody.heavy || antibody.light
        ? `<button class="secondary" data-related-ab="${esc(antibody.id)}">Find related sequences</button>`
        : "";
    const familyScope =
      antibody.heavy && antibody.light
        ? "paired"
        : antibody.heavy
          ? "heavy"
          : antibody.light
            ? "light"
            : "";
    const familyAction = familyScope
      ? `<select class="sort" data-family-threshold aria-label="Sequence cluster identity"><option value="99">99%</option><option value="95" selected>95%</option><option value="90">90%</option></select><button class="secondary" data-family-ab="${esc(antibody.id)}" data-family-scope="${familyScope}">View sequence cluster</button>`
      : "";
    const otherEvidence = `${renderConstructContext(antibody)}${targets || '<div class="meta">No direct target evidence stored.</div>'}${functional}${negative}`;
    const structuresSection = structureHtml
      ? `<details class="record-section"><summary>Structures</summary><div class="record-section-body">${structureContext}${structureHtml}<div class="structure-viewer-slot"></div></div></details>`
      : '<div class="structure-viewer-slot"></div>';
    const literatureSection = literature
      ? `<details class="record-section"><summary>Literature</summary><div class="record-section-body">${literature}</div></details>`
      : "";
    return `<div class="record-disclosure"><details class="record-section"><summary>Other target evidence</summary><div class="record-section-body">${otherEvidence}</div></details>${structuresSection}${literatureSection}<details class="record-section"><summary>Provenance &amp; source conflicts</summary><div class="record-section-body">${conflictHtml}<div class="section-title">Record provenance</div><div class="evidence-list">${provenance || '<div class="meta">No source-record details stored.</div>'}</div><div class="section-title" style="margin-top:16px">Indexed provenance timeline</div><div class="timeline"><div class="meta">Ordered only by explicit source record dates; this is not a novelty or invention timeline.</div>${timeline}</div></div></details><details class="record-section"><summary>More actions</summary><div class="record-section-body"><div class="family-slot"></div><div class="detail-actions"><button class="secondary" data-copy-ab-url="${esc(antibody.id)}">Copy antibody URL</button>${relatedAction}${familyAction}<select class="sort" data-report-category aria-label="Correction category"><option value="wrong target">Wrong target</option><option value="wrong sequence">Wrong sequence</option><option value="wrong pairing">Wrong pairing</option><option value="broken source">Broken source</option><option value="duplicate">Duplicate</option></select><button class="secondary" data-report-record="${esc(antibody.id)}">Report this record</button></div></div></details></div>`;
  }

  async function showSequenceCluster(button) {
    const antibodyId = button.dataset.familyAb;
    const scope = button.dataset.familyScope;
    const threshold = button
      .closest(".detail-actions")
      .querySelector("[data-family-threshold]").value;
    const slot = button.closest(".detail").querySelector(".family-slot");
    slot.innerHTML = loadingState("Loading sequence cluster…");
    try {
      const lookup = await getJSON(
        `${DATA_ROOT}/clusters/${scope}/${threshold}/lookup/${antibodyShardFromId(antibodyId)}.json`,
      );
      const clusterId = lookup[antibodyId];
      if (!clusterId) {
        slot.innerHTML = `<div class="family-panel"><strong>No multi-member ${esc(threshold)}% ${esc(scope)} sequence cluster</strong><span>This record is a singleton in the indexed candidate graph.</span></div>`;
        return;
      }
      const index = await getJSON(`${DATA_ROOT}/clusters/${scope}/${threshold}/index.json`);
      const cluster = index.clusters.find(item => item.id === clusterId);
      if (!cluster) throw new Error("Cluster lookup is inconsistent with its index");
      const memberIds = cluster.members.slice(0, 100);
      const records = await fetchFullRecords(memberIds);
      const members = memberIds
        .map(
          id =>
            `<a href="${esc(buildViewUrl({ ab: id }))}" data-ab-link="${esc(id)}">${esc(records.get(id)?.name || id)}</a>`,
        )
        .join("");
      slot.innerHTML = `<div class="family-panel"><strong>${fmt(cluster.size)} members · ${esc(threshold)}% ${esc(scope)} sequence cluster</strong><span>Global edit identity, ≥90% length coverage, single-link clustering. Sequence similarity only—not a clonal lineage.</span><div class="target-pills">${members}</div>${cluster.size > 100 ? `<span>First 100 of ${fmt(cluster.size)} members shown.</span>` : ""}</div>`;
    } catch (error) {
      slot.innerHTML = errorState("Sequence cluster unavailable", error.message);
    }
  }

  async function loadCardDetails(card) {
    if (card.dataset.loaded === "true") return;
    const antibody = await fetchAntibody(card.dataset.ab, card.dataset.shard);
    if (!antibody) return;
    const sequenceSlot = card.querySelector(".sequence-slot");
    const sourceNtAvailable = Boolean(
      antibody.vh_nt_source ||
      antibody.vl_nt_source ||
      antibody.heavy_nt_source ||
      antibody.light_nt_source,
    );
    const query = state.sequenceAlignmentQuery;
    const alignment = query
      ? `<div class="section-title" style="margin-top:16px">Query alignment</div>${query.heavy ? alignmentHtml("Query", antibody.name || antibody.id, query.heavy, antibody.heavy, "VH / VHH") : ""}${query.light ? alignmentHtml("Query", antibody.name || antibody.id, query.light, antibody.light, "VL") : ""}`
      : "";
    sequenceSlot.innerHTML = `<div class="section-title">Sequence quality</div>${sequenceQualityHtml(antibody)}<div class="section-title" style="margin-top:16px">Source-defensible sequence annotation</div><div class="annotation-note">Only source-supplied regions and gene calls are shown. Boundaries and numbering are never inferred.</div>${annotatedSequenceBlock("VH / VHH", antibody.heavy, antibody.cdrh3, antibody.heavy_v, antibody.heavy_j)}${sourceRegionCopies("VH", { CDR1: antibody.cdrh1, CDR2: antibody.cdrh2, CDR3: antibody.cdrh3 })}${annotatedSequenceBlock("VL", antibody.light, antibody.cdrl3, antibody.light_v, antibody.light_j)}${sourceRegionCopies("VL", { CDR1: antibody.cdrl1, CDR2: antibody.cdrl2, CDR3: antibody.cdrl3 })}${!antibody.heavy && !antibody.light ? '<div class="meta">No variable-domain sequence is available in the imported record.</div>' : ""}<div class="nucleotide-state ${sourceNtAvailable ? "available" : ""}">${sourceNtAvailable ? "Source nucleotide sequence available" : "No source nucleotide sequence available"}</div>${alignment}`;
    card.querySelector(".full-record-slot").innerHTML = renderFullRecordSlot(antibody);
    card.dataset.loaded = "true";
  }

  const conservativeResidueGroups = [
    "STA",
    "NEQK",
    "NHQK",
    "NDEQ",
    "QHRK",
    "MILV",
    "MILF",
    "HY",
    "FYW",
  ];

  function proteinScore(left, right) {
    if (left === right) return 2;
    return conservativeResidueGroups.some(group => group.includes(left) && group.includes(right))
      ? 1
      : -1;
  }

  function alignProteins(reference, query) {
    const rows = reference.length + 1;
    const columns = query.length + 1;
    const scores = Array.from({ length: rows }, () => new Int16Array(columns));
    const trace = Array.from({ length: rows }, () => new Uint8Array(columns));
    for (let row = 1; row < rows; row += 1) {
      scores[row][0] = row * -2;
      trace[row][0] = 2;
    }
    for (let column = 1; column < columns; column += 1) {
      scores[0][column] = column * -2;
      trace[0][column] = 3;
    }
    for (let row = 1; row < rows; row += 1) {
      for (let column = 1; column < columns; column += 1) {
        const diagonal =
          scores[row - 1][column - 1] + proteinScore(reference[row - 1], query[column - 1]);
        const up = scores[row - 1][column] - 2;
        const left = scores[row][column - 1] - 2;
        const best = Math.max(diagonal, up, left);
        scores[row][column] = best;
        trace[row][column] = best === diagonal ? 1 : best === up ? 2 : 3;
      }
    }
    let row = reference.length;
    let column = query.length;
    let alignedReference = "";
    let alignedQuery = "";
    while (row || column) {
      const direction = trace[row][column];
      if (direction === 1) {
        alignedReference = reference[row - 1] + alignedReference;
        alignedQuery = query[column - 1] + alignedQuery;
        row -= 1;
        column -= 1;
      } else if (direction === 2) {
        alignedReference = reference[row - 1] + alignedReference;
        alignedQuery = `-${alignedQuery}`;
        row -= 1;
      } else {
        alignedReference = `-${alignedReference}`;
        alignedQuery = query[column - 1] + alignedQuery;
        column -= 1;
      }
    }
    let matches = 0;
    let aligned = 0;
    let coveredReference = 0;
    for (let index = 0; index < alignedReference.length; index += 1) {
      const left = alignedReference[index];
      const right = alignedQuery[index];
      if (left !== "-" && right !== "-") {
        aligned += 1;
        if (left === right) matches += 1;
      }
      if (left !== "-" && right !== "-") coveredReference += 1;
    }
    return {
      reference: alignedReference,
      query: alignedQuery,
      identity: aligned ? (matches / aligned) * 100 : 0,
      coverage: reference.length ? (coveredReference / reference.length) * 100 : 0,
      score: scores[reference.length][query.length],
    };
  }

  function alignmentHtml(referenceName, queryName, reference, query, chain) {
    if (!reference || !query) {
      return `<div class="alignment-missing">${esc(chain)} alignment unavailable: ${!reference ? "reference chain missing" : "comparison chain missing"}.</div>`;
    }
    const alignment = alignProteins(reference, query);
    const blocks = [];
    for (let offset = 0; offset < alignment.reference.length; offset += 70) {
      const left = alignment.reference.slice(offset, offset + 70);
      const right = alignment.query.slice(offset, offset + 70);
      const match = [...left]
        .map((residue, index) => {
          if (residue === right[index]) return "|";
          return residue !== "-" && right[index] !== "-" && proteinScore(residue, right[index]) > 0
            ? ":"
            : " ";
        })
        .join("");
      blocks.push(
        `<pre><span>${esc(referenceName.slice(0, 14).padEnd(14))}</span> ${esc(left)}\n<span>${" ".repeat(14)}</span> ${esc(match)}\n<span>${esc(queryName.slice(0, 14).padEnd(14))}</span> ${esc(right)}</pre>`,
      );
    }
    return `<section class="alignment"><div class="alignment-title"><strong>${esc(chain)}</strong><span>${alignment.identity.toFixed(1)}% identity · ${alignment.coverage.toFixed(1)}% reference coverage · global AA alignment (match 2, conservative 1, mismatch −1, gap −2)</span></div>${blocks.join("")}</section>`;
  }

  async function compareSelectedAntibodies() {
    const ids = [...state.selectedIds];
    if (ids.length < 2 || ids.length > 10) return;
    await withBusyButton(els.compareSelected, "Loading…", async () => {
      const records = await fetchFullRecords(ids);
      const antibodies = ids.map(id => records.get(id)).filter(Boolean);
      if (antibodies.length < 2) return;
      const reference = antibodies[0];
      const tableRows = antibodies
        .map(antibody => {
          const structures = structureCollections(antibody);
          return `<tr><th><a href="${esc(buildViewUrl({ ab: antibody.id }))}">${esc(antibody.name || antibody.id)}</a><small>${esc(antibody.id)}</small></th><td>${antibody.heavy && antibody.light ? "VH + VL" : antibody.heavy ? "VH/VHH only" : antibody.light ? "VL only" : "No full chain"}</td><td>${esc(directTargetNames(antibody).slice(0, 4).join(" · ") || "—")}</td><td><code>${esc(antibody.cdrh3 || "—")}</code></td><td>${fmt(antibody.source_record_count || (antibody.source_records || []).length)}</td><td>${structures.exact.length ? `${structures.exact.length} exact` : structures.homologous.length ? `${structures.homologous.length} homologous` : "—"}</td></tr>`;
        })
        .join("");
      const alignments = antibodies
        .slice(1)
        .map(
          antibody =>
            `<div class="comparison-pair"><h4>${esc(antibody.name || antibody.id)} vs ${esc(reference.name || reference.id)}</h4>${alignmentHtml(reference.name || reference.id, antibody.name || antibody.id, reference.heavy, antibody.heavy, "VH / VHH")}${alignmentHtml(reference.name || reference.id, antibody.name || antibody.id, reference.light, antibody.light, "VL")}</div>`,
        )
        .join("");
      els.comparisonContent.innerHTML = `<div class="comparison-table-wrap"><table class="comparison-table"><thead><tr><th>Antibody</th><th>Pairing</th><th>Direct targets</th><th>CDRH3</th><th>Source records</th><th>Structures</th></tr></thead><tbody>${tableRows}</tbody></table></div><div class="comparison-method">The first selected antibody is the alignment reference. Alignments are display comparisons, not evidence of shared specificity or lineage.</div>${alignments}`;
      els.comparisonPanel.hidden = false;
      els.comparisonPanel.scrollIntoView({ behavior: "smooth", block: "start" });
    });
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
    await loadDescendantRows(activeTargetPageCategory());
    apply();
  }

  function populateAdvancedFilters() {
    els.resultSourceFilter.innerHTML = '<option value="">Any source</option>';
    for (const [key, source] of Object.entries(state.manifest?.sources || {})) {
      els.resultSourceFilter.insertAdjacentHTML(
        "beforeend",
        `<option value="${esc(key)}">${esc(source.name || key)}</option>`,
      );
    }
    els.resultEvidenceFilter.innerHTML = '<option value="">Any class</option>';
    for (const evidence of Object.keys(evidenceRank)) {
      els.resultEvidenceFilter.insertAdjacentHTML(
        "beforeend",
        `<option value="${esc(evidence)}">${esc(evidence.replaceAll("_", " "))}</option>`,
      );
    }
  }

  async function handleAdvancedFilters() {
    state.shown = PAGE_RENDER_SIZE;
    if (state.mode === "target" && state.selected) {
      els.targetMeta.textContent = "Loading all result pages for global filters…";
      await loadAllTargetPages(activeTargetPageCategory());
    }
    apply();
  }

  async function clearAdvancedFilters() {
    els.resultSourceFilter.value = "";
    els.resultSpeciesFilter.value = "";
    els.resultSequenceFilter.value = "";
    els.resultEvidenceFilter.value = "";
    await handleAdvancedFilters();
    saveWorkspaceLocal();
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

  const exportFields = [
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
    "source_nucleotide_status",
    "sequence_pairing",
    "sequence_completeness",
    "ambiguous_residues",
    "source_record_source",
    "source_record_id",
  ];

  function exportRecord(row, antibody) {
    const quality = derivedSequenceQuality(antibody);
    return {
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
      legacy_associated_target_annotations: directTargetNames(antibody).length ? "" : "",
      relationships: row.relationships.join(";"),
      evidence: row.evidence.join(";"),
      source_nucleotide_status:
        antibody.vh_nt_source ||
        antibody.vl_nt_source ||
        antibody.heavy_nt_source ||
        antibody.light_nt_source
          ? "source nucleotide available"
          : "no source nucleotide sequence available",
      sequence_pairing: quality.pairing,
      sequence_completeness: quality.completeness,
      ambiguous_residues: (quality.ambiguous_residues || []).join(";"),
      source_record_source: row.exportSourceRecord?.source || "",
      source_record_id: row.exportSourceRecord?.record_id || "",
    };
  }

  function buildCSV(rows, fullRecords) {
    const output = [exportFields.join(",")];
    for (const row of rows) {
      const record = exportRecord(row, fullRecords.get(row.antibody.id) || {});
      output.push(exportFields.map(field => csvEscape(record[field])).join(","));
    }
    return output.join("\n");
  }

  const fastaValue = value =>
    String(value || "unknown")
      .trim()
      .replace(/[^A-Za-z0-9_.-]+/g, "_")
      .slice(0, 80);

  const wrapSequence = sequence =>
    String(sequence || "")
      .match(/.{1,80}/g)
      ?.join("\n") || "";

  function sourceNucleotide(antibody, chain) {
    return chain === "VH"
      ? antibody.vh_nt_source || antibody.heavy_nt_source || ""
      : antibody.vl_nt_source || antibody.light_nt_source || "";
  }

  const codingPresets = {
    human_common_v1: {
      label: "Human common-codon v1",
      table: {
        A: "GCC",
        C: "TGC",
        D: "GAC",
        E: "GAG",
        F: "TTC",
        G: "GGC",
        H: "CAC",
        I: "ATC",
        K: "AAG",
        L: "CTG",
        M: "ATG",
        N: "AAC",
        P: "CCC",
        Q: "CAG",
        R: "AGA",
        S: "AGC",
        T: "ACC",
        V: "GTG",
        W: "TGG",
        Y: "TAC",
      },
    },
    ecoli_common_v1: {
      label: "E. coli common-codon v1",
      table: {
        A: "GCG",
        C: "TGC",
        D: "GAT",
        E: "GAA",
        F: "TTT",
        G: "GGC",
        H: "CAT",
        I: "ATT",
        K: "AAA",
        L: "CTG",
        M: "ATG",
        N: "AAC",
        P: "CCG",
        Q: "CAG",
        R: "CGT",
        S: "AGC",
        T: "ACC",
        V: "GTG",
        W: "TGG",
        Y: "TAT",
      },
    },
  };

  function backTranslate(sequence, presetKey) {
    const preset = codingPresets[presetKey];
    if (!preset) throw new Error("Unknown coding-DNA preset.");
    const invalid = [...new Set(String(sequence || "").match(/[^ACDEFGHIKLMNPQRSTVWY]/g) || [])];
    if (invalid.length)
      throw new Error(
        `Cannot back-translate ambiguous/non-codable residues: ${invalid.join(", ")}`,
      );
    return [...sequence].map(residue => preset.table[residue]).join("");
  }

  async function generatedCodingFiles(rows, fullRecords, presetKey) {
    const preset = codingPresets[presetKey];
    const fasta = [];
    const records = [];
    const errors = [];
    for (const row of rows) {
      const antibody = fullRecords.get(row.antibody.id) || {};
      for (const [chain, aminoAcids] of [
        ["VH", antibody.heavy],
        ["VL", antibody.light],
      ]) {
        if (!aminoAcids) continue;
        try {
          const codingDna = backTranslate(aminoAcids, presetKey);
          const inputDigest = await sha256Hex(aminoAcids);
          fasta.push(
            `>pairs_id=${fastaValue(row.antibody.id)}|chain=${chain}|sequence=GENERATED_FROM_AA|preset=${presetKey}|fragment=variable_domain|aa_sha256=${inputDigest}\n${wrapSequence(codingDna)}`,
          );
          records.push({
            pairs_id: row.antibody.id,
            chain,
            input_aa_sha256: inputDigest,
            input_aa_length: aminoAcids.length,
            output_nt_length: codingDna.length,
          });
        } catch (error) {
          errors.push({ pairs_id: row.antibody.id, chain, error: error.message });
        }
      }
    }
    return {
      "GENERATED_FROM_AA.fasta": fasta.join("\n"),
      "generation-manifest.json": JSON.stringify(
        {
          format: "PAIRS-synthetic-coding-DNA",
          format_version: 1,
          generated_at: new Date().toISOString(),
          preset: {
            id: presetKey,
            label: preset.label,
            algorithm: "one deterministic common codon per canonical amino acid",
          },
          semantics:
            "Synthetic variable-domain fragment generated from amino-acid sequence; not source or experimental nucleotide data and not a complete expression construct.",
          pairs_snapshot: state.manifest?.snapshot || null,
          records,
          errors,
        },
        null,
        2,
      ),
      "README.txt":
        "GENERATED_FROM_AA.fasta is synthetic back-translation of variable-domain amino-acid sequences. It is not source DNA, is not expression optimization, and lacks signal peptide, constant region, regulatory elements, and cloning context.\n",
    };
  }

  function fastaFor(rows, fullRecords, mode = "aa") {
    const output = [];
    for (const row of rows) {
      const antibody = fullRecords.get(row.antibody.id) || {};
      const target = directTargetNames(antibody)[0] || "unspecified";
      const header = chain =>
        `pairs_id=${fastaValue(row.antibody.id)}|chain=${chain}|name=${fastaValue(antibody.name || row.antibody.name)}|target=${fastaValue(target)}${row.exportSourceRecord ? `|source=${fastaValue(row.exportSourceRecord.source)}|source_record=${fastaValue(row.exportSourceRecord.record_id)}` : ""}`;
      const entries =
        mode === "heavy"
          ? [["VH", antibody.heavy]]
          : mode === "light"
            ? [["VL", antibody.light]]
            : mode === "cdr"
              ? [
                  ["CDRH3", antibody.cdrh3],
                  ["CDRL3", antibody.cdrl3],
                ]
              : mode === "nt"
                ? [
                    ["VH_SOURCE_NT", sourceNucleotide(antibody, "VH")],
                    ["VL_SOURCE_NT", sourceNucleotide(antibody, "VL")],
                  ]
                : [
                    ["VH", antibody.heavy],
                    ["VL", antibody.light],
                  ];
      for (const [chain, sequence] of entries) {
        if (!sequence) continue;
        const kind = mode === "nt" ? "|sequence=source_nucleotide" : "|sequence=amino_acid";
        output.push(`>${header(chain)}${kind}\n${wrapSequence(sequence)}`);
      }
    }
    return output.join("\n");
  }

  function exportManifest(rows, fullRecords, scope, dataset = null) {
    return {
      format: "PAIRS-export",
      format_version: 1,
      exported_at: new Date().toISOString(),
      pairs_snapshot: state.manifest?.snapshot || null,
      pairs_schema_version: state.manifest?.schema_version || SUPPORTED_SCHEMA,
      query:
        scope === "selected"
          ? {
              mode: "cross-view selection",
              current_target_id: state.selected?.id || null,
              current_target_name: state.selected?.name || null,
            }
          : state.selected
            ? { target_id: state.selected.id, target_name: state.selected.name }
            : { mode: state.mode, sequence_query: state.sequenceQueryLabel || null },
      filter: state.filter,
      scope,
      record_count: rows.length,
      unique_antibody_count: new Set(rows.map(row => row.antibody.id)).size,
      antibody_ids: rows.map(row => row.antibody.id),
      redundancy: dataset || { mode: "current filtered sequence entities" },
      sources: Object.fromEntries(
        Object.entries(state.manifest?.sources || {}).map(([key, source]) => [
          key,
          {
            name: source.name,
            last_modified: source.last_modified || null,
            sha256: source.sha256 || null,
          },
        ]),
      ),
      nucleotide: {
        semantics: "Only source-reported nucleotide sequences; never back-translated DNA.",
        records_available: rows.filter(row => {
          const antibody = fullRecords.get(row.antibody.id) || {};
          return sourceNucleotide(antibody, "VH") || sourceNucleotide(antibody, "VL");
        }).length,
      },
    };
  }

  function provenancePayload(rows, fullRecords) {
    return rows.map(row => {
      const antibody = fullRecords.get(row.antibody.id) || {};
      return {
        pairs_id: row.antibody.id,
        name: antibody.name || row.antibody.name,
        source_records: antibody.source_records || [],
        direct_targets: antibody.direct_targets || [],
        functional_targets: antibody.functional_targets || [],
        negative_evidence: antibody.negative_evidence || [],
        literature_mentions: antibody.literature_mentions || [],
        nucleotide_provenance: antibody.nucleotide_provenance || {},
        source_conflicts: antibody.source_conflicts || [],
        export_source_record: row.exportSourceRecord || null,
      };
    });
  }

  function citationRecords(rows, fullRecords) {
    const records = new Map();
    for (const row of rows) {
      const antibody = fullRecords.get(row.antibody.id) || {};
      for (const record of antibody.source_records || []) {
        if (!record.reference) continue;
        const url = record.record_url || record.source_url || "";
        const identity = [record.source, record.record_id, record.reference, url].join("\u001f");
        if (!records.has(identity)) {
          records.set(identity, {
            source: record.source || "",
            source_record_id: record.record_id || "",
            reference: String(record.reference),
            url,
            pairs_ids: new Set(),
          });
        }
        records.get(identity).pairs_ids.add(row.antibody.id);
      }
    }
    return [...records.values()].map(record => ({
      ...record,
      pairs_ids: [...record.pairs_ids].sort(),
    }));
  }

  function citationBib(rows, fullRecords) {
    const bibValue = value =>
      String(value || "")
        .replaceAll("\\", "\\textbackslash{}")
        .replace(/([{}%&#_$])/g, "\\$1")
        .replace(/[\r\n]+/g, " ");
    return citationRecords(rows, fullRecords)
      .map((record, index) => {
        const key =
          `PAIRS_${record.source || "source"}_${record.source_record_id || index + 1}`.replace(
            /[^A-Za-z0-9_:.-]/g,
            "_",
          );
        return `@misc{${key},\n  note = {${bibValue(record.reference)}},\n  howpublished = {${bibValue(record.source || "Public antibody source")}},\n  url = {${bibValue(record.url)}},\n  annote = {${bibValue(`PAIRS IDs: ${record.pairs_ids.join(", ")}`)}}\n}`;
      })
      .join("\n\n");
  }

  function citationRis(rows, fullRecords) {
    return citationRecords(rows, fullRecords)
      .map(
        record =>
          `TY  - GEN\nTI  - ${record.reference.replace(/[\r\n]+/g, " ")}\nDB  - ${record.source || "PAIRS upstream source"}\n${record.source_record_id ? `ID  - ${record.source_record_id}\n` : ""}${record.url ? `UR  - ${record.url}\n` : ""}N1  - PAIRS IDs: ${record.pairs_ids.join(", ")}; snapshot: ${state.manifest?.snapshot || "unknown"}\nER  -`,
      )
      .join("\n\n");
  }

  function citationCSV(rows, fullRecords) {
    const fields = ["source", "source_record_id", "reference", "url", "pairs_ids"];
    return [
      fields.join(","),
      ...citationRecords(rows, fullRecords).map(record =>
        fields
          .map(field =>
            csvEscape(field === "pairs_ids" ? record.pairs_ids.join(";") : record[field]),
          )
          .join(","),
      ),
    ].join("\n");
  }

  let crcTable;
  function crc32(bytes) {
    if (!crcTable) {
      crcTable = Array.from({ length: 256 }, (_, value) => {
        let crc = value;
        for (let bit = 0; bit < 8; bit += 1) crc = crc & 1 ? 0xedb88320 ^ (crc >>> 1) : crc >>> 1;
        return crc >>> 0;
      });
    }
    let crc = 0xffffffff;
    for (const byte of bytes) crc = crcTable[(crc ^ byte) & 0xff] ^ (crc >>> 8);
    return (crc ^ 0xffffffff) >>> 0;
  }

  function zipBytes(files) {
    const encoder = new TextEncoder();
    const chunks = [];
    const directory = [];
    let offset = 0;
    const u16 = value => new Uint8Array([value & 255, (value >>> 8) & 255]);
    const u32 = value =>
      new Uint8Array([
        value & 255,
        (value >>> 8) & 255,
        (value >>> 16) & 255,
        (value >>> 24) & 255,
      ]);
    const join = parts => {
      const output = new Uint8Array(parts.reduce((sum, part) => sum + part.length, 0));
      let cursor = 0;
      for (const part of parts) {
        output.set(part, cursor);
        cursor += part.length;
      }
      return output;
    };
    for (const [filename, content] of Object.entries(files)) {
      const name = encoder.encode(filename);
      const data = encoder.encode(content);
      const crc = crc32(data);
      const local = join([
        u32(0x04034b50),
        u16(20),
        u16(0x0800),
        u16(0),
        u16(0),
        u16(0),
        u32(crc),
        u32(data.length),
        u32(data.length),
        u16(name.length),
        u16(0),
        name,
        data,
      ]);
      chunks.push(local);
      directory.push(
        join([
          u32(0x02014b50),
          u16(20),
          u16(20),
          u16(0x0800),
          u16(0),
          u16(0),
          u16(0),
          u32(crc),
          u32(data.length),
          u32(data.length),
          u16(name.length),
          u16(0),
          u16(0),
          u16(0),
          u16(0),
          u32(0),
          u32(offset),
          name,
        ]),
      );
      offset += local.length;
    }
    const directoryBytes = join(directory);
    const body = join(chunks);
    const end = join([
      u32(0x06054b50),
      u16(0),
      u16(0),
      u16(directory.length),
      u16(directory.length),
      u32(directoryBytes.length),
      u32(body.length),
      u16(0),
    ]);
    return join([body, directoryBytes, end]);
  }

  async function prepareRows(scope) {
    if (scope === "selected") {
      const ids = [...state.selectedIds];
      const fullRecords = await fetchFullRecords(ids);
      const rows = ids
        .filter(id => fullRecords.has(id))
        .map(id => {
          const antibody = fullRecords.get(id);
          return {
            antibody: summaryFromFull(antibody),
            relationships: [],
            evidence: [],
            sources: antibody.sources || [],
            interactions: [],
          };
        });
      const mode = els.datasetMode.value;
      if (mode === "source_records") {
        const expanded = rows.flatMap(row => {
          const records = fullRecords.get(row.antibody.id)?.source_records || [];
          return records.length
            ? records.map(record => ({ ...row, exportSourceRecord: record }))
            : [row];
        });
        return {
          rows: expanded,
          fullRecords,
          dataset: {
            mode: "source_records",
            selected_entity_count: rows.length,
            exported_record_count: expanded.length,
          },
        };
      }
      if (mode === "unique_pairs") {
        const seen = new Map();
        const unique = [];
        const discarded = [];
        for (const row of rows) {
          const antibody = fullRecords.get(row.antibody.id) || {};
          const key =
            antibody.heavy && antibody.light
              ? `pair:${antibody.heavy}|${antibody.light}`
              : `source_scoped:${row.antibody.id}`;
          if (seen.has(key)) {
            discarded.push({ antibody_id: row.antibody.id, representative_id: seen.get(key) });
          } else {
            seen.set(key, row.antibody.id);
            unique.push(row);
          }
        }
        return {
          rows: unique,
          fullRecords,
          dataset: {
            mode: "unique_exact_vh_vl_pairs",
            algorithm: "exact normalized VH+VL equality; incomplete chains remain source-scoped",
            selected_entity_count: rows.length,
            exported_entity_count: unique.length,
            discarded,
          },
        };
      }
      if (/^cluster_(99|95|90)$/.test(mode)) {
        const threshold = mode.slice(-2);
        const clustered = await Promise.all(
          rows.map(async row => {
            const antibody = fullRecords.get(row.antibody.id) || {};
            const scope =
              antibody.heavy && antibody.light
                ? "paired"
                : antibody.heavy
                  ? "heavy"
                  : antibody.light
                    ? "light"
                    : "";
            if (!scope) return { row, key: `singleton:${row.antibody.id}`, scope: "none" };
            const lookup = await getJSON(
              `${DATA_ROOT}/clusters/${scope}/${threshold}/lookup/${antibodyShardFromId(row.antibody.id)}.json`,
            );
            return {
              row,
              key: lookup[row.antibody.id] || `singleton:${row.antibody.id}`,
              scope,
            };
          }),
        );
        const seen = new Map();
        const unique = [];
        const discarded = [];
        for (const item of clustered) {
          if (seen.has(item.key)) {
            discarded.push({
              antibody_id: item.row.antibody.id,
              representative_id: seen.get(item.key),
              cluster_id: item.key,
            });
          } else {
            seen.set(item.key, item.row.antibody.id);
            unique.push(item.row);
          }
        }
        return {
          rows: unique,
          fullRecords,
          dataset: {
            mode: `sequence_clusters_${threshold}`,
            algorithm:
              "single-link global edit identity; >=90% length coverage; paired records must pass VH and VL independently; sequence clusters are not lineages",
            selected_entity_count: rows.length,
            exported_entity_count: unique.length,
            discarded,
          },
        };
      }
      return {
        rows,
        fullRecords,
        dataset: { mode: "pairs_sequence_entities", selected_entity_count: rows.length },
      };
    }
    const rows = await allRowsForCurrentView();
    return { rows, fullRecords: await fetchFullRecords(rows.map(row => row.antibody.id)) };
  }

  async function exportCSV() {
    await withBusyButton(els.csv, "Preparing…", async () => {
      const { rows, fullRecords } = await prepareRows("filtered");
      downloadBlob(
        `${exportBaseName()}.csv`,
        buildCSV(rows, fullRecords),
        "text/csv;charset=utf-8",
      );
    });
  }

  async function exportFASTA() {
    await withBusyButton(els.fasta, "Preparing…", async () => {
      const { rows, fullRecords } = await prepareRows("filtered");
      downloadBlob(
        `${exportBaseName()}.fasta`,
        fastaFor(rows, fullRecords),
        "text/plain;charset=utf-8",
      );
    });
  }

  async function exportSelected() {
    const mode = els.selectionExportType.value;
    await withBusyButton(
      els.confirmExportSelected || els.exportSelected,
      "Preparing…",
      async () => {
        const { rows, fullRecords, dataset } = await prepareRows("selected");
        const base = "pairs-selected";
        if (mode === "csv") {
          downloadBlob(`${base}.csv`, buildCSV(rows, fullRecords), "text/csv;charset=utf-8");
        } else if (mode === "json") {
          const payload = {
            manifest: exportManifest(rows, fullRecords, "selected", dataset),
            records: provenancePayload(rows, fullRecords),
          };
          downloadBlob(`${base}.json`, JSON.stringify(payload, null, 2), "application/json");
        } else if (mode === "bundle") {
          const nt = fastaFor(rows, fullRecords, "nt");
          const files = {
            "sequences.fasta": fastaFor(rows, fullRecords),
            "heavy.fasta": fastaFor(rows, fullRecords, "heavy"),
            "light.fasta": fastaFor(rows, fullRecords, "light"),
            "cdr3.fasta": fastaFor(rows, fullRecords, "cdr"),
            "metadata.csv": buildCSV(rows, fullRecords),
            "provenance.json": JSON.stringify(provenancePayload(rows, fullRecords), null, 2),
            "references.bib": citationBib(rows, fullRecords),
            "references.ris": citationRis(rows, fullRecords),
            "references.csv": citationCSV(rows, fullRecords),
            "manifest.json": JSON.stringify(
              exportManifest(rows, fullRecords, "selected", dataset),
              null,
              2,
            ),
            "README.txt": nt
              ? "source-nucleotide.fasta contains only nucleotide sequences reported by an upstream source.\n"
              : "No source nucleotide sequences are available for these records. No DNA was generated from amino-acid sequences.\n",
          };
          if (nt) files["source-nucleotide.fasta"] = nt;
          downloadBlob(`${base}.zip`, zipBytes(files), "application/zip");
        } else if (["bib", "ris", "citations_csv"].includes(mode)) {
          const content =
            mode === "bib"
              ? citationBib(rows, fullRecords)
              : mode === "ris"
                ? citationRis(rows, fullRecords)
                : citationCSV(rows, fullRecords);
          const extension = mode === "citations_csv" ? "csv" : mode;
          downloadBlob(
            `${base}-references.${extension}`,
            content,
            mode === "citations_csv" ? "text/csv;charset=utf-8" : "text/plain;charset=utf-8",
          );
        } else if (mode === "generated_nt") {
          const files = await generatedCodingFiles(rows, fullRecords, els.codingPreset.value);
          if (!files["GENERATED_FROM_AA.fasta"])
            window.alert("No selected canonical amino-acid sequences could be generated.");
          else downloadBlob(`${base}-generated-dna.zip`, zipBytes(files), "application/zip");
        } else {
          const content = fastaFor(rows, fullRecords, mode);
          if (mode === "nt" && !content) {
            window.alert("No source nucleotide sequence is available for the selected records.");
            return;
          }
          downloadBlob(`${base}-${mode}.fasta`, content, "text/plain;charset=utf-8");
        }
      },
    );
    if (els.selectionExportSheet) {
      els.selectionExportSheet.open = false;
      els.exportSelected?.setAttribute("aria-expanded", "false");
    }
  }

  async function downloadAntibodyFASTA(antibodyId, button) {
    await withBusyButton(button, "…", async () => {
      const antibody = await fetchAntibody(antibodyId);
      if (!antibody) return;
      const row = {
        antibody: summaryFromFull(antibody),
        relationships: [],
        evidence: [],
        sources: antibody.sources || [],
        interactions: [],
      };
      downloadBlob(
        `pairs-${fastaValue(antibody.name || antibodyId).toLowerCase()}.fasta`,
        fastaFor([row], new Map([[antibodyId, antibody]])),
        "text/plain;charset=utf-8",
      );
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

  function feedbackRecordIssueUrl(antibodyId, category) {
    if (!location.hostname.endsWith(".github.io")) return "";
    const owner = location.hostname.split(".")[0];
    const repo = location.pathname.split("/").filter(Boolean)[0];
    if (!owner || !repo) return "";
    const title = `PAIRS record report: ${category} — ${antibodyId}`;
    const body = `PAIRS ID: ${antibodyId}\nCategory: ${category}\nSnapshot: ${state.manifest?.snapshot || "unknown"}\n\nPlease describe the issue and cite the upstream source record where possible. This report is not scientific evidence and does not itself correct the indexed record.\n`;
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
    // Do not draw arbitrary edges: the index has no target relationship graph.
    // Browse remains available through the ranked target grid.
    return;
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
    state.browseView = "grid";
    const field = false;
    els.browseGrid?.classList.toggle("active", !field);
    els.browseField?.classList.toggle("active", field);
    els.browseGrid?.setAttribute("aria-pressed", String(!field));
    els.browseField?.setAttribute("aria-pressed", String(field));
    els.targetGrid.hidden = field;
    if (els.bindingField) els.bindingField.hidden = true;
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

  function parseDelimitedText(text, delimiter) {
    const rows = [];
    let values = [];
    let value = "";
    let quoted = false;
    const input = String(text ?? "");
    for (let index = 0; index < input.length; index += 1) {
      const character = input[index];
      if (character === '"') {
        if (quoted && input[index + 1] === '"') {
          value += '"';
          index += 1;
        } else quoted = !quoted;
      } else if (character === delimiter && !quoted) {
        values.push(value.trim());
        value = "";
      } else if ((character === "\n" || character === "\r") && !quoted) {
        if (character === "\r" && input[index + 1] === "\n") index += 1;
        values.push(value.trim());
        if (values.some(Boolean)) rows.push(values);
        values = [];
        value = "";
      } else value += character;
    }
    values.push(value.trim());
    if (values.some(Boolean)) rows.push(values);
    if (quoted) throw new Error("Unclosed quoted field in batch CSV.");
    return rows;
  }

  function parseBatchCsv(value) {
    const firstLine = String(value ?? "").split(/\r?\n/, 1)[0];
    const delimiter = firstLine.includes("\t") ? "\t" : ",";
    const table = parseDelimitedText(value, delimiter);
    if (table.length < 2) return null;
    const headers = table[0].map(header =>
      header
        .trim()
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "_"),
    );
    const findColumn = names => names.map(name => headers.indexOf(name)).find(index => index >= 0);
    const nameColumn = findColumn(["name", "query_id", "id", "clone"]);
    const heavyColumn = findColumn(["vh", "heavy", "heavy_sequence", "vh_sequence"]);
    const lightColumn = findColumn(["vl", "light", "light_sequence", "vl_sequence"]);
    if (heavyColumn === undefined && lightColumn === undefined) return null;
    return table.slice(1).map((values, index) => {
      const query_id = values[nameColumn] || `pair_${index + 1}`;
      try {
        const query_heavy =
          heavyColumn === undefined || !values[heavyColumn]
            ? ""
            : normalizePastedSequence(values[heavyColumn]);
        const query_light =
          lightColumn === undefined || !values[lightColumn]
            ? ""
            : normalizePastedSequence(values[lightColumn]);
        if (!query_heavy && !query_light) throw new Error("No VH or VL sequence.");
        return {
          query_id,
          query_heavy,
          query_light,
          normalized_sequence: query_heavy || query_light,
          length: query_heavy.length + query_light.length,
          query_type: query_heavy && query_light ? "paired_vh_vl" : query_heavy ? "heavy" : "light",
          search_key: `pair:${query_heavy}|${query_light}`,
          validation_status: "valid",
        };
      } catch (error) {
        return {
          query_id,
          query_heavy: "",
          query_light: "",
          normalized_sequence: "",
          length: 0,
          query_type: "paired_vh_vl",
          search_key: `invalid:${index}`,
          validation_status: error.message,
        };
      }
    });
  }

  function parseBatchInput(value) {
    let csvRows;
    try {
      csvRows = parseBatchCsv(value);
    } catch (error) {
      return [
        {
          query_id: "batch_input",
          normalized_sequence: "",
          length: 0,
          query_type: "unknown",
          search_key: "invalid:batch_input",
          validation_status: error.message,
        },
      ];
    }
    if (csvRows) return csvRows;
    const raw = String(value ?? "");
    const records = raw.includes(">")
      ? parseBatchFasta(raw)
      : raw
          .split(/\r?\n/)
          .filter(line => line.trim())
          .flatMap((line, index) => {
            try {
              const normalized_sequence = normalizePastedSequence(line);
              return [
                {
                  query_id: `sequence_${index + 1}`,
                  normalized_sequence,
                  length: normalized_sequence.length,
                  validation_status: "valid",
                },
              ];
            } catch (error) {
              return [
                {
                  query_id: `sequence_${index + 1}`,
                  normalized_sequence: "",
                  length: 0,
                  validation_status: error.message,
                },
              ];
            }
          });
    return records.map(row => ({
      ...row,
      query_type: els.sequenceType.value,
      search_key: `sequence:${row.normalized_sequence}`,
    }));
  }

  function renderBatchRows() {
    els.batchResults.innerHTML = state.batchRows
      .map(
        row =>
          `<div class="batch-row"><span>${esc(row.query_id)}${row.query_heavy && row.query_light ? " · VH+VL" : ""}</span><strong>${esc(row.match_type)}</strong><span>${row.antibody_id ? `<a href="?ab=${esc(row.antibody_id)}">${esc(row.antibody_name)}</a>` : "—"}</span><span>${esc(row.heavy_identity ? `VH ${row.heavy_identity.toFixed(1)}% / ${row.heavy_coverage.toFixed(1)}% coverage${row.light_identity ? ` · VL ${row.light_identity.toFixed(1)}% / ${row.light_coverage.toFixed(1)}% coverage` : ""}` : row.matched_sequence || row.validation_status || "—")}</span></div>`,
      )
      .join("");
  }
  async function runBatchSearch() {
    const type = els.sequenceType.value,
      records = parseBatchInput(els.batchQuery.value).map(row => {
        const requested =
          type === "paired_similarity"
            ? [row.query_heavy || row.normalized_sequence, row.query_light]
            : type.endsWith("_similarity")
              ? [row.query_heavy || row.query_light || row.normalized_sequence]
              : [];
        const similarityError =
          type === "paired_similarity" && !requested[1]
            ? "Paired similarity requires explicit VH and VL columns."
            : type === "heavy_similarity" && row.query_light && !row.query_heavy
              ? "VH similarity requires a VH sequence."
              : type === "light_similarity" && row.query_heavy && !row.query_light
                ? "VL similarity requires a VL sequence."
                : requested.some(sequence => sequence && sequence.length < 40)
                  ? "Full-chain similarity requires at least 40 amino acids per requested chain."
                  : "";
        return {
          ...row,
          query_type: type.endsWith("_similarity") ? type : row.query_type,
          validation_status:
            row.validation_status === "valid" && similarityError
              ? similarityError
              : row.validation_status,
        };
      }),
      unique = new Map(
        records.filter(row => row.validation_status === "valid").map(row => [row.search_key, row]),
      );
    els.batchSearch.disabled = true;
    els.batchProgress.textContent = `Parsed ${records.length} queries; searching ${unique.size} unique sequence inputs…`;
    try {
      const found = new Map();
      await Promise.all(
        [...unique.values()].map(async row => {
          let hits;
          if (row.query_heavy || row.query_light) {
            if (type === "paired_similarity" && row.query_heavy && row.query_light) {
              const [heavyMatches, lightMatches] = await Promise.all([
                alignedSimilarityCandidates(row.query_heavy, "heavy"),
                alignedSimilarityCandidates(row.query_light, "light"),
              ]);
              const lightById = new Map(lightMatches.map(match => [match.id, match]));
              hits = heavyMatches
                .filter(match => lightById.has(match.id))
                .map(match => {
                  const lightMatch = lightById.get(match.id);
                  return {
                    antibody_ids: [match.id],
                    sequence: [match.antibody.heavy, match.antibody.light].join(" | "),
                    distance: 0,
                    heavy_identity: match.identity,
                    heavy_coverage: match.coverage,
                    light_identity: lightMatch.identity,
                    light_coverage: lightMatch.coverage,
                  };
                });
            } else if (
              (type === "heavy_similarity" && row.query_heavy) ||
              (type === "light_similarity" && row.query_light)
            ) {
              const field = type === "heavy_similarity" ? "heavy" : "light";
              const matches = await alignedSimilarityCandidates(
                field === "heavy" ? row.query_heavy : row.query_light,
                field,
              );
              hits = matches.map(match => ({
                antibody_ids: [match.id],
                sequence: match.antibody[field],
                distance: 0,
                [`${field}_identity`]: match.identity,
                [`${field}_coverage`]: match.coverage,
              }));
            } else {
              const heavyHits = row.query_heavy
                ? (await lookupSequence(row.query_heavy)).filter(hit => hit.field === "heavy")
                : [];
              const lightHits = row.query_light
                ? (await lookupSequence(row.query_light)).filter(hit => hit.field === "light")
                : [];
              const combined =
                row.query_heavy && row.query_light
                  ? combineSequenceMatches(heavyHits, lightHits)
                  : combineSequenceMatches(row.query_heavy ? heavyHits : lightHits);
              hits = combined.map(hit => ({
                antibody_ids: [hit.id],
                sequence: [row.query_heavy, row.query_light].filter(Boolean).join(" | "),
                distance: 0,
              }));
            }
          } else {
            if (["heavy_similarity", "light_similarity"].includes(type)) {
              const field = type === "heavy_similarity" ? "heavy" : "light";
              hits = (await alignedSimilarityCandidates(row.normalized_sequence, field)).map(
                match => ({
                  antibody_ids: [match.id],
                  sequence: match.antibody[field],
                  distance: 0,
                  [`${field}_identity`]: match.identity,
                  [`${field}_coverage`]: match.coverage,
                }),
              );
            } else {
              hits = ["cdrh3", "cdrl3"].includes(type)
                ? await lookupCdr(row.normalized_sequence, type, els.nearMatches.checked)
                : (await lookupSequence(row.normalized_sequence))
                    .filter(hit => type === "auto" || hit.field === type)
                    .map(hit => ({
                      antibody_ids: [hit.id],
                      sequence: row.normalized_sequence,
                      distance: 0,
                    }));
            }
          }
          found.set(
            row.search_key,
            hits.flatMap(hit => hit.antibody_ids.map(id => ({ ...hit, id }))),
          );
        }),
      );
      const full = await fetchFullRecords([
        ...new Set([...found.values()].flatMap(hits => hits.map(hit => hit.id))),
      ]);
      state.batchRows = records.flatMap(row => {
        if (row.validation_status !== "valid") return [{ ...row, match_type: "INVALID" }];
        const hits = found.get(row.search_key) || [];
        return hits.length
          ? hits.map(hit => ({
              ...row,
              match_type:
                hit.heavy_identity || hit.light_identity
                  ? "SIMILAR"
                  : hit.distance
                    ? "NEAR"
                    : "EXACT",
              matched_sequence: hit.sequence,
              antibody_id: hit.id,
              antibody_name: full.get(hit.id)?.name || hit.id,
              edit_distance: hit.distance,
              heavy_identity: hit.heavy_identity,
              heavy_coverage: hit.heavy_coverage,
              light_identity: hit.light_identity,
              light_coverage: hit.light_coverage,
              direct_targets: directTargetNames(full.get(hit.id) || {}).join("; "),
              legacy_associated_target_annotations: directTargetNames(full.get(hit.id) || {}).length
                ? ""
                : "",
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
      "query_vh",
      "query_vl",
      "query_pairing",
      "query_length",
      "query_type",
      "match_type",
      "matched_sequence",
      "edit_distance",
      "heavy_identity",
      "heavy_coverage",
      "light_identity",
      "light_coverage",
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
                query_vh: row.query_heavy,
                query_vl: row.query_light,
                query_pairing: row.query_heavy && row.query_light ? "paired" : "single_chain",
                query_length: row.length,
                query_type: row.query_type || els.sequenceType.value,
                match_type: row.match_type,
                matched_sequence: row.matched_sequence,
                edit_distance: row.edit_distance,
                heavy_identity: row.heavy_identity,
                heavy_coverage: row.heavy_coverage,
                light_identity: row.light_identity,
                light_coverage: row.light_coverage,
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

  const modalFocusableSelector =
    'button:not([disabled]), [href], input:not([disabled]), select:not([disabled]), textarea:not([disabled]), [tabindex]:not([tabindex="-1"])';

  function modalFocusable() {
    return els.modal ? [...els.modal.querySelectorAll(modalFocusableSelector)] : [];
  }

  function openSourcesModal(opener = els.sourceStatusBtn) {
    if (!els.modal) return;
    state.modalOpener = opener || document.activeElement;
    els.modal.classList.add("open");
    els.modal.setAttribute("aria-hidden", "false");
    const focusTarget = els.closeModal || modalFocusable()[0];
    requestAnimationFrame(() => focusTarget?.focus());
  }

  function closeSourcesModal(restoreFocus = true) {
    if (!els.modal) return;
    els.modal.classList.remove("open");
    els.modal.setAttribute("aria-hidden", "true");
    const opener = state.modalOpener;
    state.modalOpener = null;
    if (restoreFocus && opener && typeof opener.focus === "function" && opener.isConnected)
      requestAnimationFrame(() => opener.focus());
  }

  function trapSourcesModalFocus(event) {
    if (!els.modal?.classList.contains("open") || event.key !== "Tab") return;
    const focusables = modalFocusable();
    if (!focusables.length) {
      event.preventDefault();
      return;
    }
    const first = focusables[0];
    const last = focusables[focusables.length - 1];
    if (
      event.shiftKey &&
      (document.activeElement === first || !els.modal.contains(document.activeElement))
    ) {
      event.preventDefault();
      last.focus();
    } else if (
      !event.shiftKey &&
      (document.activeElement === last || !els.modal.contains(document.activeElement))
    ) {
      event.preventDefault();
      first.focus();
    }
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
      $("#warningSourceBtn")?.addEventListener("click", event =>
        openSourcesModal(event.currentTarget),
      );
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
      populateAdvancedFilters();
      restoreWorkspace();
      renderStatus();
      renderSourceStatus();
      renderHeroField();
      // Binding Field was a positional illustration, not a measured
      // relationship view. Keep the catalogue in its defensible grid mode.
      if (els.browseField) {
        els.browseField.hidden = true;
        els.browseField.disabled = true;
        els.browseField.setAttribute("aria-hidden", "true");
      }
      if (els.bindingField) els.bindingField.hidden = true;
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
        else renderStaleEntity("target", targetId);
      } else if (query) {
        syncSearchInputs(query);
        await search();
      }
      if (state.mode === "target" && state.selected && state.pendingWorkspaceFilter !== "all") {
        await handleFilter(state.pendingWorkspaceFilter);
        saveWorkspaceLocal();
      }
    } catch (error) {
      els.status.innerHTML = `<span class="bad">Dataset failed to load: ${esc(error.message)} Serve the repository over HTTP rather than opening index.html directly.</span>`;
      console.error(error);
    }
  }

  els.q.addEventListener("input", () => {
    syncSearchInputs(els.q.value);
    showSuggestions();
  });
  els.headerSearchInput?.addEventListener("input", () => {
    syncSearchInputs(els.headerSearchInput.value);
  });
  els.headerSearchForm?.addEventListener("submit", async event => {
    event.preventDefault();
    syncSearchInputs(els.headerSearchInput.value);
    await search();
  });
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
  els.sequenceIntent?.addEventListener("change", syncSequenceType);
  els.sequenceChain?.addEventListener("change", syncSequenceType);
  [els.textMode, els.sequenceMode].filter(Boolean).forEach(tab =>
    tab.addEventListener("keydown", event => {
      if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
      const tabs = [els.textMode, els.sequenceMode].filter(Boolean);
      if (!tabs.length) return;
      event.preventDefault();
      const current = Math.max(0, tabs.indexOf(event.currentTarget));
      const next =
        event.key === "Home"
          ? 0
          : event.key === "End"
            ? tabs.length - 1
            : (current + (event.key === "ArrowRight" ? 1 : -1) + tabs.length) % tabs.length;
      const nextTab = tabs[next];
      setSearchMode(nextTab === els.sequenceMode ? "sequence" : "text", { focusInput: false });
      nextTab.focus();
    }),
  );
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
      syncSearchInputs(chip.dataset.q);
      search();
    }),
  );

  els.filters.addEventListener("click", async event => {
    const filter = event.target.closest(".filter");
    if (filter && !filter.hidden) {
      await handleFilter(filter.dataset.filter);
      saveWorkspaceLocal();
    }
  });
  els.entitySummary.addEventListener("click", event => {
    const targetLink = event.target.closest("[data-target-link]");
    if (!targetLink) return;
    event.preventDefault();
    const target = state.targets.find(item => item.id === targetLink.dataset.targetLink);
    if (target) selectTarget(target, true);
  });
  els.entitySummary.addEventListener("change", event => {
    const scope = event.target.closest("[data-include-descendants]");
    if (scope) setDescendantScope(scope.checked);
  });
  for (const filter of [
    els.resultSourceFilter,
    els.resultSpeciesFilter,
    els.resultSequenceFilter,
    els.resultEvidenceFilter,
  ]) {
    filter.addEventListener("change", async () => {
      await handleAdvancedFilters();
      saveWorkspaceLocal();
    });
  }
  els.clearAdvancedFilters.addEventListener("click", clearAdvancedFilters);
  els.sort.addEventListener("change", async () => {
    await handleSort();
    saveWorkspaceLocal();
  });
  els.loadMore.addEventListener("click", handleLoadMore);

  els.results.addEventListener("click", async event => {
    const selection = event.target.closest("[data-select-ab]");
    if (selection) {
      setSelected(selection.dataset.selectAb, selection.checked);
      return;
    }

    const evidence = event.target.closest("[data-open-evidence]");
    if (evidence) {
      const card = evidence.closest(".card");
      card.classList.add("open");
      const expandButton = card.querySelector(".expand");
      setExpandButtonState(expandButton, true);
      await loadCardDetails(card);
      card
        .querySelector(".evidence-list")
        ?.scrollIntoView({ behavior: "smooth", block: "nearest" });
      return;
    }

    const download = event.target.closest("[data-download-ab]");
    if (download) {
      await downloadAntibodyFASTA(download.dataset.downloadAb, download);
      return;
    }

    const related = event.target.closest("[data-related-ab]");
    if (related) {
      await findRelatedSequences(related.dataset.relatedAb, related);
      return;
    }

    const structure = event.target.closest("[data-view-structure]");
    if (structure) {
      const pdbId = structure.dataset.viewStructure;
      if (!/^[0-9][A-Z0-9]{3}$/.test(pdbId)) return;
      const slot = structure.closest(".detail").querySelector(".structure-viewer-slot");
      const tier =
        structure.dataset.structureTier === "exact"
          ? "exact sequence"
          : structure.dataset.structureTier === "homologous"
            ? "homologous sequence"
            : "identity tier unavailable";
      slot.innerHTML = `<div class="structure-viewer"><div><strong>${esc(pdbId)} · ${esc(tier)}</strong><button class="secondary compact" type="button" data-close-structure>Close</button></div><iframe loading="lazy" title="3D structure ${esc(pdbId)}" src="https://molstar.org/viewer/?pdb=${encodeURIComponent(pdbId)}" referrerpolicy="no-referrer" allow="fullscreen"></iframe><span>Sequence-level structure only; PAIRS does not assert that this entry is an antibody–target complex. If embedding is blocked, <a href="https://www.rcsb.org/structure/${encodeURIComponent(pdbId)}" target="_blank" rel="noopener">open RCSB PDB</a>.</span></div>`;
      slot.scrollIntoView({ behavior: "smooth", block: "nearest" });
      return;
    }

    const closeStructure = event.target.closest("[data-close-structure]");
    if (closeStructure) {
      closeStructure.closest(".structure-viewer-slot").innerHTML = "";
      return;
    }

    const family = event.target.closest("[data-family-ab]");
    if (family) {
      await showSequenceCluster(family);
      return;
    }

    const expand = event.target.closest(".expand");
    if (expand) {
      const card = expand.closest(".card");
      card.classList.toggle("open");
      setExpandButtonState(expand, card.classList.contains("open"));
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
      return;
    }

    const reportRecord = event.target.closest("[data-report-record]");
    if (reportRecord) {
      const category =
        reportRecord.closest(".detail-actions")?.querySelector("[data-report-category]")?.value ||
        "wrong target";
      const issue = feedbackRecordIssueUrl(reportRecord.dataset.reportRecord, category);
      if (issue) window.open(issue, "_blank", "noopener");
      else
        copyText(
          `PAIRS ID: ${reportRecord.dataset.reportRecord}\nCategory: ${category}\nSnapshot: ${state.manifest?.snapshot || "unknown"}`,
          reportRecord,
        );
    }
  });

  els.copyUrl.addEventListener("click", () => copyText(location.href, els.copyUrl));
  els.csv.addEventListener("click", exportCSV);
  els.fasta.addEventListener("click", exportFASTA);
  els.selectVisible.addEventListener("click", () => {
    for (const row of state.filtered.slice(0, state.shown)) state.selectedIds.add(row.antibody.id);
    saveSelection();
    updateSelectionBar();
    els.comparisonPanel.hidden = true;
  });
  els.selectFiltered.addEventListener("click", () =>
    withBusyButton(els.selectFiltered, "Selecting…", async () => {
      const rows = await allRowsForCurrentView();
      for (const row of rows) state.selectedIds.add(row.antibody.id);
      saveSelection();
      updateSelectionBar();
      els.comparisonPanel.hidden = true;
    }),
  );
  els.clearSelection.addEventListener("click", () => {
    state.selectedIds.clear();
    saveSelection();
    updateSelectionBar();
    els.comparisonPanel.hidden = true;
  });
  els.exportSelected.addEventListener("click", () => {
    if (!els.selectionExportSheet) return;
    els.selectionExportSheet.open = !els.selectionExportSheet.open;
    els.exportSelected.setAttribute("aria-expanded", String(els.selectionExportSheet.open));
    if (els.selectionExportSheet.open) els.selectionExportType?.focus();
  });
  els.confirmExportSelected?.addEventListener("click", exportSelected);
  els.selectionExportType.addEventListener("change", () => {
    if (els.codingPresetControl)
      els.codingPresetControl.hidden = els.selectionExportType.value !== "generated_nt";
  });
  els.datasetMode.addEventListener("change", () => {
    updateDatasetPreview();
    saveWorkspaceLocal();
  });
  els.compareSelected.addEventListener("click", compareSelectedAntibodies);
  els.closeComparison.addEventListener("click", () => {
    els.comparisonPanel.hidden = true;
  });
  els.recentSearches.addEventListener("change", () => {
    const selected = els.recentSearches.value;
    els.recentSearches.value = "";
    if (selected === "") return;
    const view = state.recentViews[Number(selected)];
    if (!view) return;
    if (view.kind === "target") {
      const target = state.targets.find(item => item.id === view.id);
      if (target) selectTarget(target, true);
    } else if (view.kind === "antibody") {
      openStandaloneAntibody(view.id, true);
    }
  });
  els.exportWorkspace.addEventListener("click", exportWorkspace);
  els.importWorkspace.addEventListener("change", async () => {
    els.workspaceStatus.textContent = "Importing…";
    try {
      await importWorkspace(els.importWorkspace.files?.[0]);
    } catch (error) {
      els.workspaceStatus.textContent = error.message;
    } finally {
      els.importWorkspace.value = "";
    }
  });
  document.addEventListener("keydown", event => {
    if (els.modal?.classList.contains("open")) {
      if (event.key === "Escape") {
        event.preventDefault();
        closeSourcesModal();
      } else {
        trapSourcesModalFocus(event);
      }
      return;
    }
    if (event.key === "Escape" && !els.comparisonPanel.hidden) els.comparisonPanel.hidden = true;
  });

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

  els.sourceStatusBtn?.addEventListener("click", event => openSourcesModal(event.currentTarget));
  els.closeModal?.addEventListener("click", () => closeSourcesModal());
  els.modal.addEventListener("click", event => {
    if (event.target === els.modal) closeSourcesModal();
  });

  window.addEventListener("popstate", async () => {
    const parameters = new URLSearchParams(location.search);
    const antibodyId = parameters.get("ab");
    const targetId = parameters.get("target");
    if (antibodyId) await openStandaloneAntibody(antibodyId, false);
    else if (targetId) {
      const target = state.targets.find(item => item.id === targetId);
      if (target) await selectTarget(target, false);
      else renderStaleEntity("target", targetId);
    }
  });

  requestAnimationFrame(() => {
    document.body.classList.add("page-ready");
    setSearchMode("text");
  });
  if (els.modal) els.modal.setAttribute("aria-hidden", "true");
  syncSequenceType();
  setupHeaderSearch();
  init();
})();
