"""Update the workspace and finalize any curated notes for tagpr."""

from datetime import UTC, datetime
import os
from pathlib import Path
import subprocess


def main():
    tag = os.environ["TAGPR_NEXT_VERSION"]
    manifest = Path("Cargo.toml")
    original = manifest.read_bytes()
    subprocess.run(
        ["cargo", "set-version", "--workspace", tag.removeprefix("v")],
        check=True,
    )
    # Let tagpr apply the root version change after Cargo updates the dependencies and lockfile.
    manifest.write_bytes(original)

    changelog = Path("CHANGELOG.md")
    before, heading, after = changelog.read_text(encoding="utf-8").partition(
        "\n## Unreleased\n"
    )
    if heading and after.split("\n## ", 1)[0].strip():
        today = datetime.now(UTC).date()
        changelog.write_text(f"{before}\n## {tag} - {today}\n{after}", encoding="utf-8")


if __name__ == "__main__":
    main()
