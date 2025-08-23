# `moyo-wasm`

WASM bindings for the Rust crystal symmetry library `moyo`, for use in web apps.

## Usage

Install from a registry or local path:

```bash
pnpm add moyo-wasm
# or from cloned repo during development
pnpm add file:/path/to/moyo/moyo-wasm/pkg
```

Initialize and analyze a structure:

```ts
import init, { analyze_cell, type MoyoDataset, type MoyoCell } from 'moyo-wasm'
import wasm_url from 'moyo-wasm/moyo_wasm_bg.wasm?url'

await init(wasm_url)

// Build a JSON cell (row-major lattice matrix; fractional positions; atomic numbers)
const cell: MoyoCell = {
  lattice: { basis: [m00, m01, m02, m10, m11, m12, m20, m21, m22] },
  positions: [[fx, fy, fz], ...],
  numbers: [int, ...],
}
const result: MoyoDataset = analyze_cell(JSON.stringify(cell), 1e-4, 'Standard')
console.log(`Space group: ${result.number} (${result.hm_symbol})`)
console.log(`Hall number: ${result.hall_number}`)
console.log(`Pearson: ${result.pearson_symbol}`)
console.log(`# operations: ${result.operations.length}`)
console.log(`Wyckoffs: ${result.wyckoffs.join(', ')}`)
```

The package exports TypeScript types generated from Rust (e.g. `MoyoDataset`).

## Building

```bash
wasm-pack build moyo-wasm --target web --release
```

The package code ready for publishing is in `moyo-wasm/pkg`.

## Publish

This crate is published as an [NPM package](https://www.npmjs.com/package/moyo-wasm) with CI when a new `git` tag is pushed to the Rust monorepo.
