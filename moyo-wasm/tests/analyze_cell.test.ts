import { describe, expect, it } from "vitest";
import { analyze_cell } from "../pkg/moyo_wasm.js";

describe("analyze_cell", () => {
  // Error handling tests
  it.each([
    ["invalid JSON", "not json", /Input is not valid JSON:/],
    ["non-object", '"string"', /Input must be a JSON object/],
    [
      "missing lattice",
      '{"positions":[],"numbers":[]}',
      /Missing required field 'lattice'/,
    ],
    [
      "missing positions",
      '{"lattice":{"basis":[1,0,0,0,1,0,0,0,1]},"numbers":[]}',
      /Missing required field 'positions'/,
    ],
    [
      "missing numbers",
      '{"lattice":{"basis":[1,0,0,0,1,0,0,0,1]},"positions":[]}',
      /Missing required field 'numbers'/,
    ],
    [
      "lattice not object",
      '{"lattice":"invalid","positions":[],"numbers":[]}',
      /Field 'lattice' must be an object/,
    ],
    [
      "missing lattice.basis",
      '{"lattice":{},"positions":[],"numbers":[]}',
      /Missing required field 'lattice.basis'/,
    ],
    [
      "lattice.basis wrong length",
      '{"lattice":{"basis":[1,2,3]},"positions":[[0,0,0]],"numbers":[1]}',
      /JSON does not match expected Cell structure:/,
    ],
    [
      "positions wrong dimension",
      '{"lattice":{"basis":[1,0,0,0,1,0,0,0,1]},"positions":[[0,0]],"numbers":[1]}',
      /JSON does not match expected Cell structure:/,
    ],
    [
      "positions not 3D",
      '{"lattice":{"basis":[1,0,0,0,1,0,0,0,1]},"positions":[[0,0,0,0]],"numbers":[1]}',
      /JSON does not match expected Cell structure:/,
    ],
    [
      "positions with non-numeric values",
      '{"lattice":{"basis":[1,0,0,0,1,0,0,0,1]},"positions":[[0,0,0],["invalid",0.5,0.5]],"numbers":[1,1]}',
      /JSON does not match expected Cell structure:/,
    ],
  ])("should handle %s with helpful error", (_, input, expected_error) => {
    expect(() => analyze_cell(input, 1e-4, "Standard")).toThrow(expected_error);
  });

  // Crystal structure tests
  it.each([
    [
      "simple cubic",
      [1, 0, 0, 0, 1, 0, 0, 0, 1],
      [[0, 0, 0]],
      [1],
      "P m -3 m",
      221,
      1,
    ],
    [
      "body-centered cubic",
      [1, 0, 0, 0, 1, 0, 0, 0, 1],
      [[0, 0, 0], [0.5, 0.5, 0.5]],
      [1, 1],
      "I m -3 m",
      229,
      2,
    ],
    [
      "face-centered cubic",
      [1, 0, 0, 0, 1, 0, 0, 0, 1],
      [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
      [1, 1, 1, 1],
      "F m -3 m",
      225,
      4,
    ],
    [
      "hexagonal HCP",
      [1, 0, 0, -0.5, 0.866025, 0, 0, 0, 1.633],
      [[0, 0, 0], [1 / 3, 2 / 3, 0.5]],
      [1, 1],
      "P 6_3/m m c",
      194,
      2,
    ],
  ])(
    "should analyze %s structure",
    (
      _,
      basis,
      positions,
      numbers,
      expected_symbol,
      expected_number,
      expected_atoms,
    ) => {
      const result = analyze_cell(
        JSON.stringify({ lattice: { basis }, positions, numbers }),
        1e-4,
        "Standard",
      );
      expect(result.hm_symbol).toBe(expected_symbol);
      expect(result.number).toBe(expected_number);
      expect(result.wyckoffs.length).toBe(expected_atoms);
      expect(result.operations.length).toBeGreaterThan(0);
    },
  );

  // Settings tests
  it.each([
    ["Standard", "Standard"],
    ["Spglib", "Spglib"],
    ["spglib", "spglib"],
    ["Unknown", "Unknown"],
  ])("should work with %s setting", (setting) => {
    const cell = {
      lattice: { basis: [1, 0, 0, 0, 1, 0, 0, 0, 1] },
      positions: [[0, 0, 0]],
      numbers: [1],
    };
    const result = analyze_cell(JSON.stringify(cell), 1e-4, setting);
    expect(result.hm_symbol).toBe("P m -3 m");
  });

  // Symprec tests
  it.each([1e-3, 1e-4, 1e-5, 1e-6])(
    "should work with symprec %p",
    (symprec) => {
      const cell = {
        lattice: { basis: [1, 0, 0, 0, 1, 0, 0, 0, 1] },
        positions: [[0, 0, 0]],
        numbers: [1],
      };
      const result = analyze_cell(JSON.stringify(cell), symprec, "Standard");
      expect(result.hm_symbol).toBe("P m -3 m");
      expect(result.symprec).toBe(symprec);
    },
  );

  // Edge cases
  it.each([
    [1e-10, "very small"],
    [1e-1, "large"],
    [1e-2, "distorted"],
  ])("should handle %p symprec (%s)", (symprec) => {
    const basis = symprec === 1e-2
      ? [1.01, 0, 0, 0, 1.01, 0, 0, 0, 1.01]
      : [1, 0, 0, 0, 1, 0, 0, 0, 1];
    const cell = { lattice: { basis }, positions: [[0, 0, 0]], numbers: [1] };
    const result = analyze_cell(JSON.stringify(cell), symprec, "Standard");
    expect(result.hm_symbol).toBe("P m -3 m");
  });

  // Data integrity
  it("should return complete result structure", () => {
    const cell = {
      lattice: { basis: [2, 0, 0, 0, 2, 0, 0, 0, 2] },
      positions: [[0, 0, 0], [0.5, 0.5, 0.5]],
      numbers: [14, 14],
    };
    const result = analyze_cell(JSON.stringify(cell), 1e-4, "Standard");

    expect(result).toMatchObject({
      number: expect.any(Number),
      hall_number: expect.any(Number),
      hm_symbol: expect.any(String),
      pearson_symbol: expect.any(String),
      operations: expect.any(Array),
      wyckoffs: expect.any(Array),
      symprec: 1e-4,
    });

    expect(result.std_linear).toHaveLength(9);
    expect(result.std_origin_shift).toHaveLength(3);
    expect(result.std_rotation_matrix).toHaveLength(9);
  });
});
