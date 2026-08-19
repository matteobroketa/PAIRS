const { test, expect } = require("@playwright/test");

test("unknown text stays unchanged and does not auto-select a fuzzy result", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  const query = "definitely-not-a-real-antibody-target";
  await page.locator("#query").fill(query);
  await page.locator("#searchBtn").click();
  await expect(page.locator("#suggestions")).toContainText("No exact match");
  await expect(page.locator("#query")).toHaveValue(query);
  expect(page.url()).not.toContain(encodeURIComponent(query));
});

test("trastuzumab resolves through the generated exact antibody index", async ({ page }) => {
  await page.goto("/");
  await page.locator("#query").fill("trastuzumab");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).not.toContainText("No local match");
  await expect(page.locator("#results")).not.toContainText("No exact match");
});

test("antibody typo remains an unconfirmed no-exact result", async ({ page }) => {
  await page.goto("/");
  const query = "trastuzimab";
  await page.locator("#query").fill(query);
  await page.locator("#searchBtn").click();
  await expect(page.locator("#suggestions")).toContainText("No exact match");
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
