# moyo Space Group Browser

Static Svelte + Vite app that browses crystallographic space groups, layer groups,
and magnetic space groups via the `moyo-wasm` WebAssembly bindings.

## Development

```bash
# 1. Build moyo-wasm into ../../moyo-wasm/pkg
cd ../../moyo-wasm
npm run build

# 2. Install + dev server
cd ../apps/web
npm install
npm run dev
```

Open `http://localhost:5173/`.

## Production build

```bash
npm run build       # outputs to dist/
npm run preview     # serves dist/ locally
```

For GitHub Pages, set `VITE_BASE=/<repo>/` before `npm run build`. The CI workflow at
`.github/workflows/deploy-web.yml` does this automatically.

## Tests

```bash
npm test
```
