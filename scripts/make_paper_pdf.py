#!/usr/bin/env python3
"""
Render paper/recovering_ks1970.html to PDF, one file per language.

The paper is a single bilingual page whose language switch is a runtime toggle, so
a PDF has to pick one. The page honours ?lang=bg / ?lang=en in its URL, which is
what this script uses to produce the two editions; everything else -- page size,
margins, printing on white, keeping figure fills, not breaking a figure or table
across a page -- lives in the paper's own `@media print` block, so pressing Ctrl+P
in a browser gives the same result as running this.

Headless Chrome is used because the figures are inline SVG with CSS custom
properties and a light/dark token set: a real browser is the only thing that
resolves those the way the published page does. The theme is pinned to light for
print, since the viewer's OS setting should not decide what a printed page looks
like.

Usage:
    python3 make_paper_pdf.py                 # both languages
    python3 make_paper_pdf.py --lang en       # just one
    python3 make_paper_pdf.py --browser /path/to/chrome
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper" / "recovering_ks1970.html"

LANGS = ["bg", "en"]

# Chrome's binary is called several different things depending on the platform and
# packaging; take the first one that exists.
BROWSER_CANDIDATES = [
    "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
    "chrome", "microsoft-edge",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
]


def find_browser(explicit: str | None) -> str:
    if explicit:
        if shutil.which(explicit) or Path(explicit).exists():
            return explicit
        sys.exit(f"browser not found: {explicit}")
    for name in BROWSER_CANDIDATES:
        found = shutil.which(name) or (name if Path(name).exists() else None)
        if found:
            return found
    sys.exit(
        "No Chrome/Chromium found. Install one, or pass --browser /path/to/chrome.\n"
        "Any Chromium-based browser with --headless --print-to-pdf will do."
    )


def render(browser: str, lang: str, out_path: Path, timeout: int = 120) -> None:
    """Print one language edition to `out_path`.

    The page is copied into a temp directory with the theme pinned and the language
    stamped in, rather than passing a query string: a file:// URL with a query is
    handled inconsistently across Chromium versions, and copying also keeps the
    source paper untouched."""
    html = PAPER.read_text()
    # Pin the light theme for print and force the language, ahead of the page's own
    # script so its localStorage lookup cannot override either.
    stamp = (
        '<script>document.documentElement.setAttribute("data-theme","light");'
        f'try{{localStorage.setItem("ks1970-lang","{lang}");}}catch(e){{}}</script>'
    )
    html = html.replace("<div class=\"sheet\">", stamp + "\n<div class=\"sheet\">", 1)

    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / f"paper_{lang}.html"
        src.write_text(html)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            browser,
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--no-pdf-header-footer",
            "--generate-pdf-document-outline",
            f"--print-to-pdf={out_path}",
            src.as_uri() + f"?lang={lang}",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if not out_path.exists() or out_path.stat().st_size == 0:
            sys.exit(f"Chrome produced no PDF for {lang}.\n"
                     f"stdout: {proc.stdout[-2000:]}\nstderr: {proc.stderr[-2000:]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lang", choices=LANGS, help="only this edition (default: both)")
    ap.add_argument("--browser", help="path to a Chrome/Chromium binary")
    ap.add_argument("--outdir", type=Path, default=ROOT / "paper")
    args = ap.parse_args()

    if not PAPER.exists():
        sys.exit(f"paper not found: {PAPER}")
    browser = find_browser(args.browser)
    print(f"Using {browser}")

    for lang in ([args.lang] if args.lang else LANGS):
        out = args.outdir / f"recovering_ks1970_{lang}.pdf"
        render(browser, lang, out)
        print(f"  wrote {out.relative_to(ROOT)}  ({out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
