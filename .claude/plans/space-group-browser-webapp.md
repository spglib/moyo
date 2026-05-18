# Space Group Browser Web App (moyo-wasm)

A static web app that lets users browse crystallographic space groups (and layer / magnetic space groups) using the `moyo-wasm` bindings, deployed to GitHub Pages.

## 1. Goals & non-goals

**Goals**

- Browse all 230 space groups, 80 layer groups, and 1651 magnetic space groups offline in the browser.
- For each entry, show its identifying data (ITA number, HM symbols, Hall symbol, crystal system / family, centering, arithmetic crystal class) and the list of symmetry operations.
- Switch between settings (Standard / Spglib / specific Hall number) for the conventional space group view.
- Toggle conventional vs. primitive operations.
- Deep-linkable URLs per entry (e.g. `#/sg/227`) for sharing.
- Fully static, hosted on GitHub Pages from this monorepo.

**Non-goals (v1)**

- Crystal structure analysis (`analyze_cell`): out of scope for the "browser" view. Can be added later as a separate page.
- 3-D visualization of operations.
- Wyckoff position browser (not exposed in current wasm API).

## 2. Stack

- Build / dev: Vite (Svelte preset) + TypeScript.
- UI framework: Svelte 5 (runes API), single-page app with hash-based router.
- Styling: Tailwind CSS v4 (Vite plugin), dark mode via `class` strategy.
- WASM: `@spglib/moyo-wasm` (consumed locally during dev via `file:` link or workspace, switched to the published npm package in CI builds).
- Routing: `svelte-spa-router` (hash router; works on GitHub Pages without server config).
- Tests: Vitest + `@testing-library/svelte`; Playwright optional for one smoke test of the deployed bundle.
- Lint/format: ESLint, Prettier (with `prettier-plugin-svelte`, `prettier-plugin-tailwindcss`).

## 3. Repository layout

Add a new crate-sibling directory at the workspace root:

```
moyo/
  apps/
    web/                 # this app
      package.json
      vite.config.ts
      svelte.config.js
      tsconfig.json
      tailwind.config.ts (only if needed; Tailwind v4 uses CSS-first config)
      index.html
      src/
        app.css          # Tailwind entry
        main.ts          # mount + wasm init
        App.svelte
        lib/
          wasm.ts        # init + typed re-exports of moyo-wasm
          catalog.ts     # static arrays of valid numbers (1..230 etc.)
          format.ts      # matrix/translation pretty-print, Seitz / Jones-faithful symbols
        routes/
          Home.svelte
          SpaceGroup.svelte
          LayerGroup.svelte
          MagneticSpaceGroup.svelte
          NotFound.svelte
        components/
          SettingPicker.svelte
          PrimitiveToggle.svelte
          OperationsTable.svelte
          GroupListSidebar.svelte
          InfoGrid.svelte
        router.ts
      public/
        favicon.svg
      tests/
        format.test.ts
        wasm-smoke.test.ts
```

`apps/` is not a cargo workspace member; add it to root `.gitignore` only for `node_modules`, `dist`.

## 4. WASM integration

Vite + `wasm-pack --target web` output works directly. Key points:

- In `apps/web/package.json` depend on `"@spglib/moyo-wasm": "file:../../moyo-wasm/pkg"` for local dev; CI swaps in the published version (see deploy section).
- `lib/wasm.ts`:
  ```ts
  import init, * as moyo from '@spglib/moyo-wasm'
  import wasmUrl from '@spglib/moyo-wasm/moyo_wasm_bg.wasm?url'

  let ready: Promise<typeof moyo> | null = null
  export function getMoyo() {
    if (!ready) ready = init(wasmUrl).then(() => moyo)
    return ready
  }
  ```
- All route components `await getMoyo()` in a top-level `{#await}` block so the wasm payload is fetched once and cached by the browser.
- Catch and surface `Result<_, JsValue>` errors (the bindings return `JsValue` strings on failure) in a small `<ErrorCard>`.

## 5. Pages

### 5.1 Home (`/`)

- Three large cards: "Space groups (1-230)", "Layer groups (1-80)", "Magnetic space groups (1-1651)".
- Short intro paragraph + link to the moyo repo + license.

### 5.2 Space group detail (`#/sg/:number`)

- Header: `#227 Fd-3m` (HM short) + HM full + Hall symbol.
- Info grid (uses `space_group_type` + `hall_symbol_entry` + `arithmetic_crystal_class`):
  - ITA number, Hall number, HM short, HM full, Hall symbol.
  - Crystal family / system, lattice system, Bravais class, centering.
  - Arithmetic crystal class number + symbol, geometric crystal class, point group.
- Controls:
  - Setting picker: `Standard` (default) / `Spglib` / `HallNumber(n)` (dropdown of valid Hall numbers for the SG; computed lazily by scanning `hall_symbol_entry` 1..530 and grouping by `number`).
  - Primitive toggle (default off).
- Operations table: index, rotation (3x3 integer matrix), translation (fraction display), Seitz-like symbol if computable from matrix.
- Sidebar: collapsible list of all 230 groups, current one highlighted, filterable by crystal system.

### 5.3 Layer group detail (`#/lg/:number`)

- Same shape, using `layer_group_type`, `layer_hall_symbol_entry`, `layer_arithmetic_crystal_class`, `operations_from_layer_number`.

### 5.4 Magnetic space group detail (`#/msg/:uni_number`)

- Header: `UNI #1234` + BNS number + OG number + `magnetic_hall_symbol`.
- Info grid: uni / litvin / BNS / OG numbers, parent space-group number, construct type (1-4 with human label: `Type1`=Fedorov, `Type2`=grey, `Type3`=black-white I, `Type4`=black-white II).
- Operations table includes a `time_reversal` column.

## 6. State & routing

- `svelte-spa-router` with `replace: true` on selection so back/forward works.
- Each detail route binds `:number` (or `:uni_number`) param to a `$state` rune; route guards clamp to valid ranges and otherwise redirect to `NotFound`.
- Setting + primitive are stored in URL query (`?setting=Spglib&primitive=1`) so links are reproducible.

## 7. Styling

- Tailwind v4 CSS-first config in `app.css`:
  ```css
  @import "tailwindcss";
  @theme { --color-accent: oklch(0.62 0.18 250); }
  ```
- Use semantic Tailwind utility groups; no custom CSS unless necessary for matrix layout (a small `font-feature-settings: "tnum"` for numeric columns).
- Dark mode toggle persisted in `localStorage`.

## 8. Build & deployment (GitHub Pages)

Create `.github/workflows/deploy-web.yml`:

1. Trigger on `push` to `main` (paths: `moyo-wasm/**`, `moyo/**`, `apps/web/**`, the workflow itself) and `workflow_dispatch`.
1. Job `build`:
   - `actions/checkout`.
   - Install Rust toolchain + `wasm-pack`.
   - Build wasm: `wasm-pack build moyo-wasm --target web --release --scope spglib`.
   - Setup Node 22 + pnpm; `pnpm install` inside `apps/web` (with the `file:../../moyo-wasm/pkg` link resolving against the just-built `pkg/`).
   - `pnpm --filter @moyo/web build` with `--base=/moyo/` (or whatever the Pages base path is).
   - Upload `apps/web/dist` as the Pages artifact.
1. Job `deploy`: `actions/deploy-pages` reading the artifact. Use the `github-pages` environment + `pages: write, id-token: write` permissions.

Add `apps/web/vite.config.ts` reading `VITE_BASE` env (or compute from `process.env.GITHUB_REPOSITORY`) so local dev keeps `/` and CI emits the right base.

Add a `404.html` that is a copy of `index.html` so hash routing survives direct URL visits (Pages serves 404 then the SPA bootstraps).

## 9. Testing

- `tests/format.test.ts`: pure-TS tests for matrix / translation formatting (no wasm needed).
- `tests/wasm-smoke.test.ts`: vitest in `environment: 'jsdom'` (or `node` with `--experimental-wasm-modules`); `init()` the wasm, call `space_group_type(227)`, assert HM short is `Fd-3m`.
- CI runs `pnpm test` after build but before deploy.

## 10. Out-of-scope follow-ups

- Wyckoff position table (needs new wasm exports).
- Structure analyzer page powered by `analyze_cell` with a paste-in JSON or POSCAR drop zone.
- Search box that fuzzy-matches HM symbol -> SG number.
- Diagrammatic unit-cell view of generators.

## 11. Step-by-step execution

1. Scaffold `apps/web` with `pnpm create vite@latest apps/web -- --template svelte-ts`; commit baseline. ([structural])
1. Install Tailwind v4 + the Svelte/Vite plugins; verify `pnpm dev` renders a Tailwind-styled page. ([structural])
1. Add `svelte-spa-router` + create `routes/` skeleton with placeholder components and the `Home` card grid. ([structural])
1. Wire `lib/wasm.ts`; consume `@spglib/moyo-wasm` via `file:../../moyo-wasm/pkg`; build the wasm once locally with `just` and confirm `space_group_type(1)` returns. ([behavioral])
1. Implement `SpaceGroup.svelte` end-to-end (info grid + setting picker + primitive toggle + operations table). ([behavioral])
1. Duplicate for `LayerGroup.svelte` and `MagneticSpaceGroup.svelte`. ([behavioral])
1. Add `GroupListSidebar` with crystal-system filter; deep-link query params. ([behavioral])
1. Add unit + smoke tests; run `pnpm test`. ([behavioral])
1. Add `.github/workflows/deploy-web.yml`; enable GitHub Pages (Source: GitHub Actions) for the repo. ([structural])
1. Push, verify the deployed URL, fix base-path / 404 fallback issues. ([behavioral])
1. Update root `README.md` with a "Web demo" link. ([docs])

Each step lands as its own commit (structural changes first per Tidy First).
