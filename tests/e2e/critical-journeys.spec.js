const { test, expect } = require("@playwright/test");
const fs = require("fs");

test("unknown text stays unchanged and does not auto-select a fuzzy result", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  const query = "definitely-not-a-real-antibody?";
  await page.locator("#query").fill(query);
  await expect(page.locator("#suggestions")).toContainText("No suggestions yet");
  await expect(page.locator("#suggestions")).not.toContainText("No exact match");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).toContainText("No exact match");
  await expect(page.locator("#query")).toHaveValue(query);
  expect(page.url()).not.toContain(encodeURIComponent(query));
});

test("trastuzumab resolves through the generated exact antibody index", async ({ page }) => {
  await page.goto("/");
  await page.locator("#query").fill("trastuzumab");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).not.toContainText("No local match");
  await expect(page.locator("#results")).not.toContainText("No exact match");
  await expect(page.locator("#targetName")).toContainText("Trastuzumab");
  await expect(page.locator("#entitySummary")).toContainText("ABH71318");
});

test("autocomplete uses lossless antibody names and preserves the matched name", async ({
  page,
}) => {
  await page.goto("/");
  await page.locator("#query").fill("trastu");
  const suggestion = page
    .locator("#suggestions .suggestion[data-index]")
    .filter({ hasText: "Trastuzumab" })
    .first();
  await expect(suggestion).toBeVisible();
  await expect(page.locator("#suggestions")).not.toContainText("No exact match");
  await suggestion.click();
  await expect(page.locator("#targetName")).toContainText("Trastuzumab");
  await expect(page.locator("#entitySummary")).toContainText("ABH71318");
});

test("CD19 autocomplete stays visible above the Browse targets section", async ({
  page,
}, testInfo) => {
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  await page.locator("#query").fill("CD19");
  const suggestions = page.locator("#suggestions");
  await expect(suggestions).toHaveClass(/open/);
  await expect(suggestions.locator(".suggestion[data-index]").first()).toBeVisible();
  await page.waitForTimeout(220);

  const overlay = await page.evaluate(() => {
    const dropdown = document.querySelector("#suggestions");
    const browse = document.querySelector("#browse");
    const searchPanel = document.querySelector("#textSearchPanel");
    const searchWrap = document.querySelector(".search-wrap");
    const dropdownRect = dropdown.getBoundingClientRect();
    const browseRect = browse.getBoundingClientRect();
    const overlapTop = Math.max(dropdownRect.top, browseRect.top);
    const overlapBottom = Math.min(dropdownRect.bottom, browseRect.bottom);
    const point = {
      x: dropdownRect.left + dropdownRect.width / 2,
      y: overlapTop + Math.min(8, Math.max(1, overlapBottom - overlapTop - 1)),
    };
    const layers = document.elementsFromPoint(point.x, point.y);

    return {
      dropdown: {
        top: dropdownRect.top,
        bottom: dropdownRect.bottom,
        left: dropdownRect.left,
        right: dropdownRect.right,
      },
      browse: { top: browseRect.top, left: browseRect.left, right: browseRect.right },
      overlap: { top: overlapTop, bottom: overlapBottom },
      position: getComputedStyle(dropdown).position,
      dropdownZ: Number(getComputedStyle(dropdown).zIndex),
      searchPanelZ: Number(getComputedStyle(searchPanel).zIndex),
      searchWrapZ: Number(getComputedStyle(searchWrap).zIndex),
      topLayerWithinSuggestions: layers[0]?.closest("#suggestions") === dropdown,
      browseCoversPoint: layers[0]?.closest("#browse") === browse,
    };
  });

  expect(overlay.position).toBe("absolute");
  expect(overlay.dropdown.top).toBeLessThan(overlay.browse.top);
  expect(overlay.dropdown.bottom).toBeGreaterThan(overlay.browse.top);
  expect(overlay.overlap.bottom).toBeGreaterThan(overlay.overlap.top);
  expect(overlay.dropdown.left).toBeGreaterThanOrEqual(overlay.browse.left);
  expect(overlay.dropdown.right).toBeLessThanOrEqual(overlay.browse.right);
  expect(overlay.searchPanelZ).toBeGreaterThanOrEqual(100);
  expect(overlay.searchWrapZ).toBeGreaterThanOrEqual(100);
  expect(overlay.dropdownZ).toBeGreaterThan(overlay.searchWrapZ);
  expect(overlay.topLayerWithinSuggestions).toBe(true);
  expect(overlay.browseCoversPoint).toBe(false);

  await page.screenshot({
    path: testInfo.outputPath("cd19-autocomplete-overlay.png"),
    fullPage: false,
  });
});

test("antibody typo remains an unconfirmed no-exact result", async ({ page }) => {
  await page.goto("/");
  const query = "trastuzimab";
  await page.locator("#query").fill(query);
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).toContainText("No exact match");
  await expect(page.locator("#query")).toHaveValue(query);
  expect(page.url()).not.toContain("ab=");
  expect(page.url()).not.toContain("target=");
});

test("a shared exact antibody name opens a selection state", async ({ page }) => {
  await page.goto("/");
  await page.locator("#query").fill("scfv16");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).toContainText("Exact matches");
  expect(await page.locator("[data-exact-ab]").count()).toBeGreaterThan(1);
  expect(page.url()).not.toContain("ab=");
  await page.locator("[data-exact-ab]").first().click();
  await expect.poll(() => page.url()).toContain("ab=");
});

test("a close target suggestion opens only after an explicit click", async ({ page }) => {
  await page.goto("/");
  await page.locator("#query").fill("ERBB");
  const suggestion = page
    .locator("#suggestions .suggestion[data-index]")
    .filter({ hasText: "ERBB2" })
    .first();
  await expect(suggestion).toBeVisible();
  expect(page.url()).not.toContain("target=");
  await suggestion.click();
  await expect(page.locator("#targetName")).toContainText("ERBB2");
  await expect(page.locator("#targetMeta")).toContainText("direct positive evidence");
  await expect(page.locator("#targetMeta")).not.toContainText("source-level relationships");
  await expect(page.locator("#results .card .sources")).toHaveCount(0);
  await expect.poll(() => page.url()).toContain("target=");
});

test("FASTA export uses the active target context", async ({ page }) => {
  await page.goto("/");
  await page.locator("#query").fill("ERBB");
  await page
    .locator("#suggestions .suggestion[data-index]")
    .filter({ hasText: "ERBB2" })
    .first()
    .click();
  await expect(page.locator("#targetName")).toContainText("ERBB2");
  await page.locator("details.collection-export summary").click();
  const [download] = await Promise.all([
    page.waitForEvent("download"),
    page.locator("#fastaBtn").click(),
  ]);
  const content = fs.readFileSync(await download.path(), "utf8");
  expect(content).toContain("target=ERBB2");
});

test("exact-index outages are not reported as antibody absence", async ({ page }) => {
  await page.route("**/data/v4/antibody-exact/tr.json", route => route.abort());
  await page.goto("/");
  await page.locator("#query").fill("trastuzumab");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).toContainText("Search index unavailable");
  await expect(page.locator("#results")).toContainText(
    "PAIRS cannot determine whether this antibody is present",
  );
});

test("a pasted sequence never enters the URL", async ({ page }) => {
  await page.goto("/");
  await page.locator("#sequenceModeBtn").click();
  const sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWV";
  await page.locator("#heavySequenceQuery").fill(sequence);
  await page.locator("#sequenceIntent").selectOption("similarity");
  await page.locator("#sequenceChain").selectOption("heavy");
  await page.locator("#sequenceSearchBtn").click();
  await page.waitForTimeout(1000);
  expect(page.url()).not.toContain(sequence);
});

test("exact sequence region selection does not cross-match CDR fields", async ({ page }) => {
  await page.goto("/");
  await page.locator("#sequenceModeBtn").click();
  const cdrh3 = "ARSRDLLLFPHHALSP";
  const cdrl3 = "SSRDSSGNHWV";

  await page.locator("#heavySequenceQuery").fill(cdrh3);
  await page.locator("#sequenceChain").selectOption("heavy");
  await page.locator("#sequenceSearchBtn").click();
  await expect(page.locator("#results")).toContainText("No exact public sequence match");

  await page.locator("#sequenceChain").selectOption("cdrh3");
  await page.locator("#sequenceSearchBtn").click();
  await expect(page.locator("#results")).toContainText("ABL52827");

  await page.locator("#sequenceChain").selectOption("paired");
  await page.locator("#lightSequenceQuery").fill(cdrl3);
  await page.locator("#sequenceSearchBtn").click();
  await expect(page.locator("#results")).toContainText("No exact public sequence match");
});

test("stable antibody deep links load record details", async ({ page }) => {
  await page.goto("/?ab=ab_b4adc696667f200be5ba");
  await expect(page.locator("#targetName")).toContainText("7D11");
  await expect(page.locator("#targetMeta")).not.toContainText("target annotations");
  await expect(page.locator("#results .card")).toHaveCount(1);
});

test("search modes expose complete tab relationships", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#textModeBtn")).toHaveAttribute("aria-controls", "textSearchPanel");
  await expect(page.locator("#sequenceModeBtn")).toHaveAttribute(
    "aria-controls",
    "sequenceSearchPanel",
  );
  await expect(page.locator("#textSearchPanel")).toHaveAttribute("role", "tabpanel");
  await expect(page.locator("#sequenceSearchPanel")).toHaveAttribute("role", "tabpanel");
  await expect(page.locator("#filters")).toContainText("Evidence");
  await expect(page.locator("#filters")).toContainText("Direct positive");
  await expect(page.locator("#filters")).toContainText("Filter");
  await expect(page.locator("#browseSort option").first()).toHaveText(
    "Most direct-positive antibodies",
  );
});

test("search remains available after the hero leaves view", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  await page.evaluate(() => window.scrollTo(0, document.body.scrollHeight));
  await expect(page.locator("#headerSearchForm")).toBeVisible();
  await page.locator("#headerSearchInput").fill("ERBB2");
  await page.locator("#headerSearchInput").press("Enter");
  await expect(page.locator("#targetName")).toContainText("ERBB2");
});

test("selection and export controls appear only when needed", async ({ page }) => {
  await page.goto("/?ab=ab_b4adc696667f200be5ba");
  await expect(page.locator("#selectionBar")).toBeHidden();
  await page.locator("[data-select-ab]").check();
  await expect(page.locator("#selectionBar")).toBeVisible();
  // Chromium's mobile emulation may auto-scroll fixed controls underneath the
  // page during actionability checks even though hit-testing places the tray first.
  await page.locator("#exportSelectedBtn").click({ force: true });
  await expect(page.locator("#selectionExportSheet")).toHaveAttribute("open", "");
  await expect(page.locator("#codingPresetControl")).toBeHidden();
  await page.locator("#selectionExportType").selectOption("generated_nt");
  await expect(page.locator("#codingPresetControl")).toBeVisible();
});

test("sources dialog closes with Escape and restores focus", async ({ page }) => {
  await page.goto("/");
  const mobile = (page.viewportSize()?.width || 0) <= 620;
  if (mobile) await page.locator(".mobile-nav > summary").click();
  const opener = page.locator(mobile ? "#mobileSourceStatusBtn" : "#sourceStatusBtn");
  await opener.click();
  await expect(page.locator("#sourceModal")).toHaveClass(/open/);
  await expect(page.locator("#closeModal")).toBeFocused();
  await page.keyboard.press("Escape");
  await expect(page.locator("#sourceModal")).not.toHaveClass(/open/);
  await expect(opener).toBeFocused();
});
