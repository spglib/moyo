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

For GitHub Pages, set `VITE_BASE=/<repo>/` before `npm run build`. The CI workflow
at `.github/workflows/deploy-web.yml` does this automatically.

### One-time GitHub setup for the deploy workflow

The first time `actions/deploy-pages` runs, GitHub auto-creates a `github-pages`
environment whose "Deployment branches and tags" rule may not match `main`. If
the deploy job reports

> Branch "main" is not allowed to deploy to github-pages due to environment
> protection rules.

then, in the repo settings:

1. **Settings -> Pages**: set **Source** to **GitHub Actions**.
1. **Settings -> Environments -> github-pages -> Deployment branches and tags**:
   either pick **No restriction**, or **Selected branches and tags** and add
   `main`.

Then re-run the failed workflow or push another commit to `main`.

## Tests

```bash
just web-check     # svelte-check (typecheck)
just web-test      # vitest unit tests
```
