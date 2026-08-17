(() => {
  "use strict";

  const SUPPORTED_SCHEMA = 2;
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
    browseMore: $("#browseMore"),
    filters: $("#filters"),
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
    targetPageCache: new Map(),
    loadedTargetPages: new Set(),
    browseShown: 60,
    sequenceQueryLabel: "",
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
    })[relationship] || relationship.replaceAll("_", " ");

  const esc = value =>
    String(value ?? "").replace(
      /[&<>"']/g,
      character =>
        ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[
          character
        ],
    );

  const norm = value =>
    String(value ?? "")
      .toLowerCase()
      .normalize("NFKD")
      .replace(/[^a-z0-9]+/g, " ")
      .trim();

  const compact = value => norm(value).replaceAll(" ", "");
  const fmt = value => new Intl.NumberFormat().format(value || 0);
  const isNegative = relationship =>
    relationship.includes("does_not") || relationship.includes("not_");

  async function getJSON(url) {
    const response = await fetch(url, { cache: "no-cache" });
    if (!response.ok) {
      if (response.status === 404 && url.startsWith(DATA_ROOT)) {
        throw new Error(
          "A versioned dataset file is unavailable. The site may have been updated while this tab was open; reload the page.",
        );
      }
      throw new Error(`${url}: HTTP ${response.status}`);
    }
    return response.json();
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
        for (const alias of antibody.aliases || []) score = Math.max(score, textScore(alias, query));
        return [score, antibody];
      })
      .filter(([score]) => score >= 0)
      .sort((left, right) => right[0] - left[0])
      .slice(0, limit)
      .map(([, antibody]) => antibody);
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
      html.push(`<div class="suggestion" role="option" data-index="${index}"><div class="s-main"><strong>${esc(target.name)}</strong><span>${esc(aliases || target.sources.join(" · "))}</span></div><span class="count">${fmt(target.result_count)} antibodies</span></div>`);
    }
    for (const antibody of antibodies) {
      const index = state.suggestionItems.push({ kind: "antibody", antibody }) - 1;
      html.push(`<div class="suggestion" role="option" data-index="${index}"><div class="s-main"><strong>${esc(antibody.name)}</strong><span>Antibody · ${esc((antibody.targets || []).slice(0, 4).map(target => target.name).join(" · ") || antibody.sources.join(" · "))}</span></div><span class="count">${antibody.paired ? "VH + VL" : "sequence record"}</span></div>`);
    }

    els.suggestions.innerHTML = html.join("") || renderNoSuggestion(query);
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
    return withoutHeaders.toUpperCase().replace(/[^A-Z*]/g, "");
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
      therapeutic_status: antibody.therapeutic_status || "",
      shard: antibodyShardFromId(antibody.id),
    };
  }

  async function runSequenceSearch(firstRaw, secondRaw = "") {
    const first = normalizePastedSequence(firstRaw);
    const second = normalizePastedSequence(secondRaw);
    if (!first) {
      showSearchError("Paste at least one amino-acid sequence.");
      return;
    }

    els.sequenceSearch.disabled = true;
    const original = els.sequenceSearch.textContent;
    els.sequenceSearch.textContent = "Searching…";
    try {
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
    els.results.innerHTML = `<div class="empty"><strong>No exact public sequence match</strong>PAIRS found no exact ${paired ? "paired " : ""}match for the normalized ${length}-aa query in this snapshot.<div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a></div></div>`;
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
    els.results.innerHTML = `<div class="empty"><strong>Could not run sequence search</strong>${esc(message)}</div>`;
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
    resetFilterButtons();
  }

  function resetFilterButtons() {
    $$(".filter").forEach(button =>
      button.classList.toggle("active", button.dataset.filter === "all"),
    );
  }

  function updateFilterAvailability() {
    const negative = $('.filter[data-filter="negative"]');
    if (negative) negative.hidden = state.mode !== "target";
  }

  function targetPageUrl(target, pageNumber) {
    return `${DATA_ROOT}/targets/${target.dir}/page-${String(pageNumber).padStart(3, "0")}.json`;
  }

  async function loadTargetPage(pageNumber) {
    if (!state.selected || pageNumber < 1 || pageNumber > state.selected.page_count) return [];
    const key = `${state.selected.id}:${pageNumber}`;
    if (state.targetPageCache.has(key)) return state.targetPageCache.get(key);
    const rows = await getJSON(targetPageUrl(state.selected, pageNumber));
    state.targetPageCache.set(key, rows);
    if (!state.loadedTargetPages.has(pageNumber)) {
      state.loadedTargetPages.add(pageNumber);
      state.rawResults.push(...rows);
    }
    return rows;
  }

  async function loadNextTargetPage() {
    if (!state.selected) return false;
    for (let page = 1; page <= state.selected.page_count; page += 1) {
      if (!state.loadedTargetPages.has(page)) {
        await loadTargetPage(page);
        return true;
      }
    }
    return false;
  }

  async function loadAllTargetPages() {
    if (!state.selected) return;
    const missing = [];
    for (let page = 1; page <= state.selected.page_count; page += 1) {
      if (!state.loadedTargetPages.has(page)) missing.push(page);
    }
    const rows = await Promise.all(missing.map(page => getJSON(targetPageUrl(state.selected, page))));
    rows.forEach((pageRows, index) => {
      const page = missing[index];
      const key = `${state.selected.id}:${page}`;
      state.targetPageCache.set(key, pageRows);
      state.loadedTargetPages.add(page);
      state.rawResults.push(...pageRows);
    });
  }

  async function selectTarget(target, push = true) {
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
    els.results.innerHTML = '<div class="empty">Loading first result page…</div>';
    els.loadMore.hidden = true;
    try {
      await loadTargetPage(1);
      apply();
      if (push) {
        const url = new URL(location.href);
        url.search = "";
        url.searchParams.set("target", target.id);
        history.pushState({}, "", url);
      }
      scrollToResults();
    } catch (error) {
      els.results.innerHTML = `<div class="empty"><strong>Could not load target results</strong>${esc(error.message)}</div>`;
    }
  }

  function updateTargetMeta() {
    if (!state.selected) return;
    els.targetMeta.textContent = `${fmt(state.selected.result_count)} unique antibodies · ${fmt(state.selected.count)} source-level relationships · ${state.selected.sources.join(" · ")} · loaded ${state.loadedTargetPages.size}/${state.selected.page_count} pages`;
  }

  async function openStandaloneAntibody(antibodyId, push = true) {
    closeSuggestions();
    els.main.classList.add("active");
    els.targetName.textContent = "Antibody record";
    els.targetMeta.textContent = "Loading public sequence record…";
    els.results.innerHTML = '<div class="empty">Loading antibody shard…</div>';
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
      els.results.innerHTML = `<div class="empty"><strong>Could not open antibody</strong>${esc(error.message)}</div>`;
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

    await showSuggestions();
    if (state.suggestionItems.length) {
      activateSuggestion(0);
      return;
    }
    renderTextNoMatch(query);
  }

  function renderTextNoMatch(query) {
    state.mode = "target";
    state.selected = null;
    state.rawResults = [];
    state.filtered = [];
    els.main.classList.add("active");
    els.targetName.textContent = "No local match";
    els.targetMeta.textContent = `No target or antibody matched “${query}” in this snapshot.`;
    els.summary.innerHTML = "";
    const issue = feedbackIssueUrl(query);
    els.results.innerHTML = `<div class="empty"><strong>No indexed match</strong>Try a synonym, browse the target catalogue, or report the missing query so it can be reviewed for a future alias/source update.<div class="empty-actions"><a class="secondary" href="#browse">Browse targets</a>${issue ? `<a class="secondary" href="${esc(issue)}" target="_blank" rel="noopener">Report missing target</a>` : `<button class="secondary" data-copy-missing="${esc(query)}">Copy missing-query report</button>`}</div></div>`;
    els.loadMore.hidden = true;
    scrollToResults();
  }

  function passes(result) {
    const antibody = result.antibody;
    if (state.filter === "paired") return Boolean(antibody.has_heavy && antibody.has_light);
    if (state.filter === "therapeutic") return Boolean(antibody.therapeutic_status);
    if (state.filter === "structure") return Boolean(antibody.structures?.length);
    if (state.filter === "negative") return result.relationships.some(isNegative);
    return true;
  }

  function rank(result) {
    const evidence = Math.max(...result.evidence.map(item => evidenceRank[item] ?? 0), 0);
    return (
      evidence * 100 +
      result.sources.length * 10 +
      (result.antibody.structures?.length ? 5 : 0) +
      (result.antibody.therapeutic_status ? 3 : 0)
    );
  }

  function apply() {
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
        (left, right) => rank(right) - rank(left) || left.antibody.name.localeCompare(right.antibody.name),
      );
    }
    renderSummary();
    renderResults();
    updateTargetMeta();
  }

  function renderSummary() {
    if (state.mode === "target" && state.selected?.stats) {
      const stats = state.selected.stats;
      els.summary.innerHTML = [
        [stats.unique_results, "Unique antibodies"],
        [stats.paired, "VH + VL pairs"],
        [stats.structure, "With structures"],
        [stats.negative, "Negative evidence"],
      ]
        .map(
          ([value, label]) =>
            `<div class="metric"><b>${fmt(value)}</b><span>${label}</span></div>`,
        )
        .join("");
      return;
    }

    const results = state.filtered;
    const paired = results.filter(result => result.antibody.has_heavy && result.antibody.has_light).length;
    const structures = results.filter(result => result.antibody.structures?.length).length;
    const therapeutics = results.filter(result => result.antibody.therapeutic_status).length;
    els.summary.innerHTML = [
      [results.length, "Results"],
      [paired, "VH + VL pairs"],
      [structures, "With structures"],
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
      output.push('<span class="badge strong">VH + VL</span>');
    }
    if (result.antibody.therapeutic_status) {
      output.push('<span class="badge strong">THERAPEUTIC</span>');
    }
    if (result.antibody.structures?.length) output.push('<span class="badge">PDB</span>');
    for (const field of result.match_fields || []) {
      output.push(`<span class="badge sequence-match">EXACT ${esc(field.toUpperCase())}</span>`);
    }
    for (const relationship of result.relationships.slice(0, 4)) {
      output.push(
        `<span class="badge ${isNegative(relationship) ? "negative" : ""}">${esc(relationLabel(relationship).toUpperCase())}</span>`,
      );
    }
    return output.join("");
  }

  function evidenceRows(interactions) {
    if (!interactions.length) return '<div class="meta">Open the target links below to inspect target-specific evidence.</div>';
    return interactions
      .map(interaction => {
        const source = state.manifest?.sources?.[interaction.source];
        const sourceName = source?.name || interaction.source;
        const sourceText = source?.homepage
          ? `<a href="${esc(source.homepage)}" target="_blank" rel="noopener">${esc(sourceName)}</a>`
          : esc(sourceName);
        const details = [
          interaction.evidence,
          interaction.epitope && `Epitope: ${interaction.epitope}`,
          interaction.assay && `Assay: ${interaction.assay}`,
          interaction.reference,
          interaction.note,
        ]
          .filter(Boolean)
          .join(" · ");
        return `<div class="evidence-row ${isNegative(interaction.relationship) ? "negative" : ""}"><strong>${esc(relationLabel(interaction.relationship))} · ${sourceText}</strong><span>${esc(details)}</span></div>`;
      })
      .join("");
  }

  function renderResults() {
    const subset = state.filtered.slice(0, state.shown);
    const targetHasMorePages =
      state.mode === "target" &&
      state.selected &&
      state.loadedTargetPages.size < state.selected.page_count;
    els.loadMore.hidden = state.shown >= state.filtered.length && !targetHasMorePages;

    if (!subset.length) {
      els.results.innerHTML = '<div class="empty">No loaded results match the current filter.</div>';
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
  }

  function sequenceBlock(label, value) {
    if (!value) return "";
    return `<div class="seq"><div class="seq-head"><span>${esc(label)} · ${fmt(value.length)} aa</span><button class="copy" data-copy-seq="${esc(value)}">Copy</button></div><code>${esc(value)}</code></div>`;
  }

  function renderFullRecordSlot(antibody) {
    const targets = (antibody.targets || [])
      .slice(0, 40)
      .map(
        target =>
          `<a class="target-pill" href="${esc(buildViewUrl({ target: target.id }))}" data-target-link="${esc(target.id)}">${esc(target.name)}</a>`,
      )
      .join("");
    const structures = (antibody.structures || [])
      .map(
        pdb =>
          `<a class="pdb-link" href="https://www.rcsb.org/structure/${encodeURIComponent(pdb)}" target="_blank" rel="noopener">PDB ${esc(pdb)}</a>`,
      )
      .join("");
    const provenance = (antibody.source_records || [])
      .map(record => {
        const label = state.manifest?.sources?.[record.source]?.name || record.source;
        const linked = record.source_url
          ? `<a href="${esc(record.source_url)}" target="_blank" rel="noopener">${esc(label)}</a>`
          : esc(label);
        return `<div class="evidence-row"><strong>${linked} · ${esc(record.record_id)}</strong><span>${esc(record.reference || "Source record")}</span></div>`;
      })
      .join("");
    return `<div class="section-title" style="margin-top:16px">Targets</div><div class="target-pills">${targets || '<span class="meta">No normalized target annotation.</span>'}</div>${structures ? `<div class="section-title" style="margin-top:16px">Structures</div><div class="target-pills">${structures}</div>` : ""}<div class="section-title" style="margin-top:16px">Record provenance</div><div class="evidence-list">${provenance || '<div class="meta">No source-record details stored.</div>'}</div><div class="detail-actions"><button class="secondary" data-copy-ab-url="${esc(antibody.id)}">Copy antibody URL</button></div>`;
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
    while (
      state.filtered.length < minimum &&
      state.loadedTargetPages.size < state.selected.page_count
    ) {
      const loaded = await loadNextTargetPage();
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
      state.loadedTargetPages.size < state.selected.page_count
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
      const before = state.filtered.length;
      while (state.loadedTargetPages.size < state.selected.page_count) {
        const loaded = await loadNextTargetPage();
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
      await loadAllTargetPages();
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
    const label = state.selected?.name || (state.mode === "sequence" ? "sequence-search" : "antibody");
    return `pairs-${label}`.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
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
        "structures",
        "sources",
        "targets",
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
          structures: (antibody.structures || []).join(";"),
          sources: (antibody.sources || row.sources || []).join(";"),
          targets: (antibody.targets || []).map(target => target.name).join(";"),
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
        const safeName = String(antibody.name || row.antibody.name || row.antibody.id).replace(/\s+/g, "_");
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
        button.textContent = "Copied";
        setTimeout(() => {
          button.textContent = original;
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

  function renderBrowse() {
    const query = norm(els.browseQuery.value);
    let targets = state.targets.filter(target => {
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
    els.browseMore.hidden = visible.length >= targets.length;
  }

  function renderSourceStatus() {
    if (!state.manifest) return;
    els.sourceList.innerHTML = Object.entries(state.manifest.sources)
      .map(([key, source]) => {
        const details = source.ok
          ? `${fmt(source.records)} records · ${fmt(source.interactions)} evidence rows${source.bytes ? ` · ${formatBytes(source.bytes)}` : ""}${source.discovered ? " · current link discovered at build time" : ""}`
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

  function renderStatus() {
    const stats = state.manifest.stats;
    const sourceOkay = state.manifest.sources_ok ?? Object.values(state.manifest.sources).filter(source => source.ok).length;
    const sourceExpected = state.manifest.sources_expected ?? Object.keys(state.manifest.sources).length;
    const partial = sourceOkay < sourceExpected;
    els.status.innerHTML = `<span><span class="status-dot ${partial ? "warning" : ""}"></span>${fmt(stats.antibodies)} antibodies · ${fmt(stats.interactions)} evidence records · ${fmt(stats.targets)} targets</span><span>${sourceOkay}/${sourceExpected} public sources included</span>`;
    const snapshot = new Date(state.manifest.snapshot);
    els.footer.textContent = `Snapshot ${Number.isNaN(snapshot.getTime()) ? state.manifest.snapshot : snapshot.toLocaleDateString()} · PAIRS ${state.manifest.app_version} · schema v${state.manifest.schema_version}`;

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
      setActiveSuggestion(state.activeSuggestion <= 0 ? state.suggestionItems.length - 1 : state.activeSuggestion - 1);
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
    if (copyMissing) copyText(`Missing PAIRS query: ${copyMissing.dataset.copyMissing}`, copyMissing);
  });

  document.addEventListener("click", event => {
    if (!event.target.closest(".search-wrap")) closeSuggestions();
    const copyMissing = event.target.closest("[data-copy-missing]");
    if (copyMissing) copyText(`Missing PAIRS query: ${copyMissing.dataset.copyMissing}`, copyMissing);
  });

  els.textMode.addEventListener("click", () => setSearchMode("text"));
  els.sequenceMode.addEventListener("click", () => setSearchMode("sequence"));
  els.sequenceSearch.addEventListener("click", () =>
    runSequenceSearch(els.heavySequence.value, els.lightSequence.value),
  );

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
  els.browseMore.addEventListener("click", () => {
    state.browseShown += 60;
    renderBrowse();
  });
  els.targetGrid.addEventListener("click", event => {
    const card = event.target.closest("[data-browse-target]");
    if (!card) return;
    const target = state.targets.find(item => item.id === card.dataset.browseTarget);
    if (target) selectTarget(target, true);
  });

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

  init();
})();
