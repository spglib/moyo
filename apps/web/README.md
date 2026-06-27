# moyo Crystallographic Group Viewer

Static Svelte + Vite app that browses crystallographic space groups, layer groups,
and magnetic space groups via the `moyo-wasm` WebAssembly bindings.

All commands are wrapped as `just` recipes from the repo root; the raw `npm`
commands below are equivalent and run from this directory.

## Development

```bash
# from repo root: rebuild moyo-wasm/pkg, then start the dev server
just web-install   # one-time: install npm deps
just web-dev       # runs js-build (wasm-pack) and `npm run dev`
```

Open `http://localhost:5173/`.

## Production build

```bash
just web-build     # outputs to apps/web/dist
just web-preview   # serves dist locally
```

For GitHub Pages the viewer is served under a subpath, so CI sets
`VITE_BASE=/<repo>/viewer/` before `npm run build` (the landing hub in
`apps/landing` occupies the site root). The workflow at
`.github/workflows/deploy-web.yml` does this automatically; routing is hash-based
(`svelte-spa-router`), so no server-side SPA fallback is needed.

## Tests

```bash
just web-check     # svelte-check (typecheck)
just web-test      # vitest unit tests
```
