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

// Cartesian basis vectors a, b, c; fractional positions; atomic numbers
const cell: MoyoCell = {
  lattice: { basis: [ax, ay, az, bx, by, bz, cx, cy, cz] },
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

## Standardized cells and handedness

`analyze_cell` accepts either input handedness. A passive basis change, with
fractional positions transformed consistently, preserves the physical
space-group type. Physically reflecting a crystal instead exchanges the members
of an enantiomorphic pair.

Both `std_cell` and `prim_std_cell` have right-handed bases. Their coordinate
transformations, `std_linear` and `prim_std_linear`, have negative determinants
for left-handed input and positive determinants for right-handed input.
`std_rotation_matrix` is a proper Cartesian rotation with determinant +1.

All flat matrix fields use column-major storage: element `(i, j)` is at index
`i + 3*j`. In particular, `lattice.basis` concatenates the three Cartesian basis
vectors. With these column-wise matrices, the selected coordinates satisfy
`A_selected = A_input * std_linear` and
`x_selected = inverse(std_linear) * (x_input - std_origin_shift)`.
Lattice and position refinement follow this coordinate selection, then the
Cartesian rotation is applied. See the
[returned-cell specification](https://spglib.github.io/moyo/standardization/#returned-cell-specification)
for the symmetric refinement and primitive/conventional relation.

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
