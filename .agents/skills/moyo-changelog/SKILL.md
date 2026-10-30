---
name: moyo-changelog
description: Write or review moyo CHANGELOG.md entries for a pull request, or audit and finalize the changelog for a requested release.
---

# Moyo changelog

Update the repository-root `CHANGELOG.md` from the actual changes. Use its current
layout and preserve existing entries.

## Entries for a pull request

Read the requested PR's current diff, base, and linked issue before writing.
Use the implementation and regression tests to establish what changed; an issue
may describe symptoms that require several separate fixes. Attribute only the
behavior addressed by the requested PR.

Add a bullet under `## Unreleased`. If that section is absent, insert it before
the latest release. Use the applicable crate heading, following the existing
order: `moyo`, `moyopy`, `moyoc`, `moyo-wasm`. Create only headings that need entries.

- Put core symmetry-search changes under `moyo`. Do not repeat them for each
  binding merely because the bindings use the core.
- Put binding-specific API, packaging, or compatibility changes under that
  binding. State minimum supported Python versions and ABI changes explicitly.
- Describe the user-visible change and the affected conditions. Include the
  mechanism only when it explains the result. Keep distinct fixes separate when
  they address different failure modes.
- Use the existing issue/PR reference style, such as `(#123, #124)`. Verify the
  numbers and avoid adding duplicate entries already present on the base branch.
- Do not claim that one fix solves every symptom of a broader issue. Distinguish
  incorrect results, rejected inputs, and retries with adjusted tolerances.

## Release changelog

For a requested version, identify the preceding release tag and audit all commits
and the cumulative diff through the proposed release head. Check the existing
Unreleased entries against that history, including merged sibling PRs. Do not
assume that every user-visible change already has an entry.

Include bug fixes, public API changes, and support changes. Omit routine internal
refactors, test cleanup, and dependency updates unless they affect users.

Finalize `## Unreleased` as `## vX.Y.Z - YYYY-MM-DD` using the requested version
and release date. If no date is specified, use the current date for preparation.
Follow the existing convention of not leaving an empty Unreleased section.
A changelog request alone does not authorize version bumps, tags, or publishing;
perform those only when included in the surrounding task.

## Verification

Read the final diff alongside the source changes and check each behavioral claim,
crate placement, reference, and version/date. Preserve newer entries when syncing
with the base branch. Run the repository's Markdown hooks and any checks required
by its contribution instructions. Report the changed file or PR and the checks
actually run.
