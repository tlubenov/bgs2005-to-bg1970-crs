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
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .paths import paper_dir

LANGS = ("bg", "en")

# Chrome's binary is called several different things depending on the platform and
# packaging; take the first one that exists.
BROWSER_CANDIDATES = [
    "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
    "chrome", "microsoft-edge",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
]


def find_browser(explicit: str | None = None) -> str:
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


def render(browser: str, lang: str, out_path: Path, paper: Path, timeout: int = 120) -> None:
    """Print one language edition to `out_path`.

    The page is copied into a temp directory with the theme pinned and the language
    stamped in, rather than passing a query string: a file:// URL with a query is
    handled inconsistently across Chromium versions, and copying also keeps the
    source paper untouched."""
    html = Path(paper).read_text()
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


def build_pdfs(langs=LANGS, root=None, browser: str | None = None,
               verbose: bool = True) -> list[Path]:
    """Render the paper to PDF, one file per language. Returns the paths written."""
    say = print if verbose else (lambda *a, **k: None)
    base = paper_dir(root)
    paper = base / "recovering_ks1970.html"
    if not paper.exists():
        raise FileNotFoundError(paper)
    exe = find_browser(browser)
    say(f"Using {exe}")
    written = []
    for lang in langs:
        out = base / f"recovering_ks1970_{lang}.pdf"
        render(exe, lang, out, paper)
        written.append(out)
        say(f"  wrote {out.name}  ({out.stat().st_size/1024:.0f} KB)")
    return written
