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
  await expect.poll(() => page.url()).toContain("target=");
});

test("a pasted sequence never enters the URL", async ({ page }) => {
  await page.goto("/");
  await page.locator("#sequenceModeBtn").click();
  const sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWV";
  await page.locator("#heavySequenceQuery").fill(sequence);
  await page.locator("#sequenceType").selectOption("heavy_similarity");
  await page.locator("#sequenceSearchBtn").click();
  await page.waitForTimeout(1000);
  expect(page.url()).not.toContain(sequence);
});

test("stable antibody deep links load record details", async ({ page }) => {
  await page.goto("/?ab=ab_b4adc696667f200be5ba");
  await expect(page.locator("#targetName")).toContainText("7D11");
  await expect(page.locator("#results .card")).toHaveCount(1);
});
