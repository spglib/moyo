# moyo Space Group Browser

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

For GitHub Pages, set `VITE_BASE=/<repo>/` before `npm run build`. The CI workflow
at `.github/workflows/deploy-web.yml` does this automatically.

## Tests

```bash
just web-check     # svelte-check (typecheck)
just web-test      # vitest unit tests
```
