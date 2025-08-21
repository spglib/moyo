# moyo-wasm

WASM bindings for the Rust crystal symmetry library `moyo`, for use in web apps.

## Usage (Vite/Svelte)

Install from a local path or registry:

```bash
pnpm add moyo-wasm
# or from this repo during development
pnpm add file:/path/to/moyo/moyo-wasm/pkg
```

Initialize and analyze a structure:

```ts
import init, { analyze_simple_cell } from 'moyo-wasm'
import wasm_url from 'moyo-wasm/moyo_wasm_bg.wasm?url'

await init(wasm_url)

const result = analyze_simple_cell(basis, positions, numbers, 1e-4, 'Standard')
console.log(result.number)
```

- basis: Float64Array length 9 (row-major 3x3)
- positions: Float64Array length 3N (fractional coords)
- numbers: Int32Array length N (atomic numbers)

Alternatively pass a JSON cell to `analyze_cell(json, symprec, setting)`.

## Building

```bash
wasm-pack build moyo-wasm --target web --release
```

The output package is in `moyo-wasm/pkg`.

## Publish

This crate is designed to be published as an npm package with CI when a Git tag is created in the Rust monorepo.
