# moyo-wasm

[![image](https://img.shields.io/npm/v/%40spglib%2Fmoyo-wasm)](https://www.npmjs.com/package/@spglib/moyo-wasm)

JavaScript and WebAssembly interface of [moyo](https://github.com/spglib/moyo), a fast and robust crystal symmetry finder.

- npm: <https://www.npmjs.com/package/@spglib/moyo-wasm>

## Installation

```shell
npm install @spglib/moyo-wasm
# or from a cloned repo during development
npm install file:/path/to/moyo/moyo-wasm/pkg
```

## Usage

Initialize and analyze a structure:

```ts
import init, { analyze_cell, type MoyoDataset, type MoyoCell } from '@spglib/moyo-wasm'
import wasm_url from '@spglib/moyo-wasm/moyo_wasm_bg.wasm?url'

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

## Development

Run from the repo root with `just` (or the equivalent npm commands from this directory):

```shell
just js-install   # npm install
just js-build     # wasm-pack build --target web --release --scope spglib
just js-test      # npm test
```

The package code ready for publishing is generated in `moyo-wasm/pkg`. It is published to npm by CI when a new git tag is pushed to the monorepo.

## How to cite moyo-wasm

See the citation information in [the root README](https://github.com/spglib/moyo/blob/main/README.md)
