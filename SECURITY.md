# Security

PAIRS is a static site with no login, server-side database, analytics collector, sequence-upload endpoint, or runtime call to an antibody source. Search terms and pasted sequences remain in the browser.

The build pipeline downloads public source files during CI. Source URLs are discovered from configured upstream homepages, downloads use retry/backoff, and required CSV columns are validated before data is published.

For security issues in this repository, open a private GitHub security advisory if enabled. For scientific data errors, use a normal issue and identify the source record where possible.

Do not place credentials or private dataset URLs in `config/sources.json`; GitHub Pages artifacts are public.
