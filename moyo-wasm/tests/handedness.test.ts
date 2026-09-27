import { expect, it } from "vitest";
import { analyze_cell } from "../pkg/moyo_wasm.js";

// WASM matrices are flattened column-major; basis arrays concatenate a, b, c.
function matvec(matrix: number[], vector: number[]): number[] {
  return [0, 1, 2].map((i) =>
    vector.reduce((sum, x, j) => sum + matrix[i + 3 * j] * x, 0)
  );
}

function multiply(a: number[], b: number[]): number[] {
  return [0, 1, 2].flatMap((j) => matvec(a, b.slice(3 * j, 3 * j + 3)));
}

function determinant(a: number[]): number {
  return a[0] * (a[4] * a[8] - a[7] * a[5])
    - a[3] * (a[1] * a[8] - a[7] * a[2])
    + a[6] * (a[1] * a[5] - a[4] * a[2]);
}

it.each([
  { handedness: 1, mirrored: false },
  { handedness: -1, mirrored: false },
  { handedness: 1, mirrored: true },
  { handedness: -1, mirrored: true },
])("preserves physical handedness for %j", ({ handedness, mirrored }) => {
  // P4_1 with two species, proper Cartesian rotation and a passive skew rebase
  // P = [[1,1,0],[0,1,1],[0,0,handedness]], plus a nonzero origin shift.
  const basis = [4, 0, 0, 4, 2.4, 3.2, 0, 2.4 - 4.8 * handedness, 3.2 + 3.6 * handedness];
  if (mirrored) {
    for (const i of [0, 3, 6]) basis[i] *= -1;
  }
  const positions: number[][] = [];
  const numbers: number[] = [];
  for (const [species, seed] of [[0.137, 0.271, 0.389], [0.219, 0.413, 0.157]].entries()) {
    let [x, y] = seed;
    for (let n = 0; n < 4; n++) {
      const z = (seed[2] + n / 4 - 0.19) / handedness;
      const rebasedY = y - 0.27 - z;
      positions.push([x - 0.13 - rebasedY, rebasedY, z]);
      numbers.push(species + 1);
      [x, y] = [-y, x];
    }
  }
  const result = analyze_cell(JSON.stringify({ lattice: { basis }, positions, numbers }), 1e-5, "Standard");
  expect(result.number).toBe(mirrored ? 78 : 76);
  expect(result.hall_number).toBe(mirrored ? 352 : 350);
  expect(result.operations).toHaveLength(4);
  const rotation = result.std_rotation_matrix;
  expect(determinant(rotation)).toBeCloseTo(1, 8);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const dot = [0, 1, 2].reduce((sum, k) => sum + rotation[k + 3 * i] * rotation[k + 3 * j], 0);
      expect(dot).toBeCloseTo(i === j ? 1 : 0, 8);
    }
  }
  for (const [cell, linear, shift, mapping] of [
    [result.std_cell, result.std_linear, result.std_origin_shift, undefined],
    [result.prim_std_cell, result.prim_std_linear, result.prim_std_origin_shift, result.mapping_std_prim],
  ] as const) {
    expect(cell.positions).toHaveLength(8);
    expect(determinant(cell.lattice.basis)).toBeGreaterThan(0);
    expect(determinant(linear) * determinant(basis)).toBeGreaterThan(0);
    const expectedBasis = multiply(rotation, multiply(basis, linear));
    cell.lattice.basis.forEach((x, i) => expect(x).toBeCloseTo(expectedBasis[i], 8));
    positions.forEach((position, i) => {
      const matches = cell.positions.flatMap((site, j) => {
        const transformed = matvec(linear, site);
        const match = cell.numbers[j] === numbers[i] && transformed.every((x, k) => {
          const delta = x + shift[k] - position[k];
          return Math.abs(delta - Math.round(delta)) < 1e-8;
        });
        return match ? [j] : [];
      });
      expect(matches).toHaveLength(1);
      if (mapping) expect(matches).toEqual([mapping[i]]);
    });
  }
});
