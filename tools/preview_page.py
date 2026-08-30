"""Preview docs/index.md locally, before pushing it to GitHub Pages.

    pip install grip
    python tools/preview_page.py            # then open the printed URL
    python tools/preview_page.py --port 8080

Renders through GitHub's own markdown API, so headings, tables and fenced code
come out byte-identical to what GitHub will show. Two honest caveats:

  * The live site is styled by the Primer *Jekyll theme*; this is styled like a
    GitHub README. The page chrome differs slightly - the content does not.
  * The YAML front matter is stripped before rendering, because Jekyll consumes
    it and GitHub's markdown API would print it as literal text. That is the
    only edit made to the file.

Relative image paths are served from docs/, so the light/dark <picture> swap is
testable here: switch your OS theme, or use the browser devtools
"Emulate CSS prefers-color-scheme" setting, and reload.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE = REPO_ROOT / "docs" / "index.md"
# Leading underscore so Jekyll ignores it too, if it is ever committed by accident.
PREVIEW = REPO_ROOT / "docs" / "_preview.md"


def strip_front_matter(text: str) -> str:
    """Drop a leading YAML front-matter block, leaving the body untouched."""
    if not text.startswith("---"):
        return text
    lines = text.splitlines(keepends=True)
    for index in range(1, len(lines)):
        if lines[index].rstrip() in ("---", "..."):
            return "".join(lines[index + 1:]).lstrip("\n")
    return text  # unterminated block - render as-is rather than silently eat the file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=6419)
    args = parser.parse_args()

    if not SOURCE.exists():
        print(f"not found: {SOURCE.relative_to(REPO_ROOT)}")
        return 1

    PREVIEW.write_text(strip_front_matter(SOURCE.read_text(encoding="utf-8")),
                       encoding="utf-8")
    print(f"rendering {SOURCE.relative_to(REPO_ROOT)} (front matter stripped)")
    print(f"  http://localhost:{args.port}/    - Ctrl+C to stop\n")

    try:
        # grip takes the listen address as a positional ("<host>:<port>" or just
        # "<port>"); there is no --port flag. --title overrides the default, which
        # would otherwise be the scratch filename.
        return subprocess.call([sys.executable, "-m", "grip", str(PREVIEW),
                                str(args.port), "--browser", "--title=StructScan"])
    except KeyboardInterrupt:
        return 0
    finally:
        PREVIEW.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
