const { test, expect } = require("@playwright/test");
const fs = require("fs");
const path = require("path");

// Skip animation timing so the tests exercise overlapping network requests.
test.use({ reducedMotion: "reduce" });

async function browseTarget(page, name) {
  await page.locator("#browseQuery").fill(name);
  await page.locator('[data-browse-target="target:' + name.toLowerCase() + '"]').click();
}

function deferred() {
  let resolve;
  const promise = new Promise(done => {
    resolve = done;
  });
  return { promise, resolve };
}

async function delayDataFile(page, relativePath) {
  const received = deferred();
  const released = deferred();
  const finished = deferred();
  await page.route("**/" + relativePath, async route => {
    received.resolve();
    await released.promise;
    await route.fulfill({
      contentType: "application/json",
      body: fs.readFileSync(path.join(__dirname, "../..", relativePath), "utf8"),
    });
    finished.resolve();
  });
  return {
    received: received.promise,
    async release() {
      released.resolve();
      await finished.promise;
      // Let the completed fetch and its render callbacks drain before asserting.
      await page.waitForTimeout(150);
    },
  };
}

test("Back restores the home view and Forward restores the target", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  const homeTitle = await page.title();
  await browseTarget(page, "CD19");
  await expect(page).toHaveURL(/target=target%3Acd19/);
  await expect(page.locator("#results .card").first()).toBeVisible();
  await page.goBack();
  await expect(page.locator("#main")).not.toHaveClass(/active/);
  await expect(page.locator("#results .card")).toHaveCount(0);
  await expect(page.locator("#query")).toHaveValue("");
  await expect(page).toHaveTitle(homeTitle);
  await page.goForward();
  await expect(page.locator("#targetName")).toHaveText("CD19");
  await expect(page.locator("#results .card").first()).toBeVisible();
});

test("a failed name search removes the previous entity summary", async ({ page }) => {
  await page.goto("/?ab=ab_b4adc696667f200be5ba");
  await expect(page.locator("#entitySummary")).toBeVisible();
  await page.locator("#query").fill("definitely-not-a-real-antibody?");
  await page.locator("#searchBtn").click();
  await expect(page.locator("#targetName")).toContainText("No exact match");
  await expect(page.locator("#entitySummary")).toBeHidden();
  await expect(page).toHaveTitle("No exact match — PAIRS");
});

test("late target pages cannot mix evidence into a newer target", async ({ page }) => {
  const delayed = await delayDataFile(page, "data/v4/targets/d56051e9ea1733df81b2/page-001.json");
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  await browseTarget(page, "ERBB2");
  await delayed.received;
  await browseTarget(page, "CD19");
  await expect(page).toHaveURL(/target=target%3Acd19/);
  await expect(page.locator("#results .card").first()).toBeVisible();
  const cards = await page
    .locator("#results .card")
    .evaluateAll(elements => elements.map(element => element.dataset.ab));
  const summary = await page.locator("#summaryGrid").innerText();
  await delayed.release();
  await expect(page).toHaveURL(/target=target%3Acd19/);
  await expect(page.locator("#targetName")).toHaveText("CD19");
  expect(
    await page
      .locator("#results .card")
      .evaluateAll(elements => elements.map(element => element.dataset.ab)),
  ).toEqual(cards);
  await expect(page.locator("#summaryGrid")).toHaveText(summary, { useInnerText: true });
});

test("late exact-name searches cannot replace a newer target", async ({ page }) => {
  const delayed = await delayDataFile(page, "data/v4/antibody-exact/tr.json");
  await page.goto("/");
  await expect(page.locator("#footerSnapshot")).not.toBeEmpty();
  await page.locator("#query").fill("trastuzumab");
  await page.locator("#searchBtn").click();
  await delayed.received;
  await browseTarget(page, "CD19");
  await expect(page).toHaveURL(/target=target%3Acd19/);
  await delayed.release();
  await expect(page.locator("#targetName")).toHaveText("CD19");
  await expect(page.locator("#query")).toHaveValue("CD19");
  await expect(page).toHaveURL(/target=target%3Acd19/);
});

test("typing while the catalogue loads still produces target suggestions", async ({ page }) => {
  const delayed = await delayDataFile(page, "data/v4/targets.json");
  await page.goto("/");
  await delayed.received;
  await page.locator("#query").fill("ERBB");
  await delayed.release();
  await expect(
    page.locator("#suggestions .suggestion[data-index]").filter({ hasText: "ERBB2" }).first(),
  ).toBeVisible();
});

for (const previousView of ["target=target%3Acd19", "ab=ab_b4adc696667f200be5ba"]) {
  test(`changing scope during target loading replaces the previous ${previousView.split("=")[0]} URL`, async ({ page }) => {
    await page.goto("/?" + previousView);
    await expect(page.locator("#results .card").first()).toBeVisible();
    const delayed = await delayDataFile(page, "data/v4/targets/95952a55b092a8e39824/page-001.json");
    await browseTarget(page, "HIV");
    await delayed.received;
    await page.locator("[data-include-descendants]").check();
    await delayed.release();
    await expect(page).toHaveURL(/target=target%3Ahiv/);
    expect(new URL(page.url()).searchParams.has("ab")).toBe(false);
    expect(new URL(page.url()).searchParams.get("scope")).toBe("descendants");
    await expect(page.locator("#targetName")).toHaveText("HIV");
    await page.reload();
    await expect(page.locator("#targetName")).toHaveText("HIV");
    await expect(page.locator("[data-include-descendants]")).toBeChecked();
    await expect(page.locator("#results .card").first()).toBeVisible();
  });
}
