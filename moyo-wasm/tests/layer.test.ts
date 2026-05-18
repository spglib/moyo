import { describe, expect, it } from "vitest";
import {
  layer_arithmetic_crystal_class,
  layer_group_type,
  layer_hall_symbol_entry,
  operations_from_layer_number,
} from "../pkg/moyo_wasm.js";

describe("operations_from_layer_number", () => {
  it("LG 1 (p1, Standard) returns 1 operation", () => {
    const ops = operations_from_layer_number(1, { type: "Standard" }, false);
    expect(ops.length).toBe(1);
  });

  it("LG 80 (p6/mmm, Standard) returns 24 operations", () => {
    const ops = operations_from_layer_number(80, { type: "Standard" }, false);
    expect(ops.length).toBe(24);
  });

  it("LG 10 (c211, Standard) is centered: 4 conventional / 2 primitive", () => {
    const ops_conv = operations_from_layer_number(
      10,
      { type: "Standard" },
      false,
    );
    expect(ops_conv.length).toBe(4);
    const ops_prim = operations_from_layer_number(
      10,
      { type: "Standard" },
      true,
    );
    expect(ops_prim.length).toBe(2);
  });

  it("throws on invalid layer number", () => {
    expect(() =>
      operations_from_layer_number(0, { type: "Standard" }, false),
    ).toThrow();
  });
});

describe("layer_hall_symbol_entry", () => {
  it("hall_number 1 -> LG 1 (p1)", () => {
    const entry = layer_hall_symbol_entry(1);
    expect(entry.hall_number).toBe(1);
    expect(entry.number).toBe(1);
    expect(entry.arithmetic_number).toBe(1);
    expect(entry.setting).toBe("");
    expect(entry.hall_symbol).toBe("p 1");
    expect(entry.hm_short).toBe("p 1");
    expect(entry.hm_full).toBe("p 1");
    expect(entry.centering).toBe("P");
  });

  it("hall_number 116 -> LG 80 (p6/mmm)", () => {
    const entry = layer_hall_symbol_entry(116);
    expect(entry.hall_number).toBe(116);
    expect(entry.number).toBe(80);
    expect(entry.arithmetic_number).toBe(43);
    expect(entry.hm_short).toBe("p 6/m m m");
  });

  it("throws on invalid hall_number", () => {
    expect(() => layer_hall_symbol_entry(0)).toThrow();
  });
});

describe("layer_arithmetic_crystal_class", () => {
  it("arithmetic 1 -> p1 (oblique primitive)", () => {
    const acc = layer_arithmetic_crystal_class(1);
    expect(acc.arithmetic_number).toBe(1);
    expect(acc.symbol).toBe("p1");
    expect(acc.geometric_crystal_class).toBe("1");
    expect(acc.bravais_class).toBe("mp");
    expect(acc.lattice_system).toBe("Oblique");
  });

  it("arithmetic 7 -> c211 (centered rectangular)", () => {
    const acc = layer_arithmetic_crystal_class(7);
    expect(acc.symbol).toBe("c211");
    expect(acc.bravais_class).toBe("oc");
    expect(acc.lattice_system).toBe("Rectangular");
  });

  it("arithmetic 43 -> p6/mmm (hexagonal primitive)", () => {
    const acc = layer_arithmetic_crystal_class(43);
    expect(acc.symbol).toBe("p6/mmm");
    expect(acc.geometric_crystal_class).toBe("6/mmm");
    expect(acc.bravais_class).toBe("hp");
    expect(acc.lattice_system).toBe("Hexagonal");
  });

  it("throws on invalid arithmetic_number", () => {
    expect(() => layer_arithmetic_crystal_class(0)).toThrow();
  });
});

describe("layer_group_type", () => {
  it("LG 1 (p1)", () => {
    const lgt = layer_group_type(1);
    expect(lgt.number).toBe(1);
    expect(lgt.hm_short).toBe("p 1");
    expect(lgt.hm_full).toBe("p 1");
    expect(lgt.arithmetic_number).toBe(1);
    expect(lgt.arithmetic_symbol).toBe("p1");
    expect(lgt.geometric_crystal_class).toBe("1");
    expect(lgt.bravais_class).toBe("mp");
    expect(lgt.lattice_system).toBe("Oblique");
  });

  it("LG 10 (c211, centered rectangular)", () => {
    const lgt = layer_group_type(10);
    expect(lgt.number).toBe(10);
    expect(lgt.arithmetic_symbol).toBe("c211");
    expect(lgt.geometric_crystal_class).toBe("2");
    expect(lgt.bravais_class).toBe("oc");
    expect(lgt.lattice_system).toBe("Rectangular");
  });

  it("LG 80 (p6/mmm, hexagonal)", () => {
    const lgt = layer_group_type(80);
    expect(lgt.number).toBe(80);
    expect(lgt.arithmetic_number).toBe(43);
    expect(lgt.arithmetic_symbol).toBe("p6/mmm");
    expect(lgt.geometric_crystal_class).toBe("6/mmm");
    expect(lgt.bravais_class).toBe("hp");
    expect(lgt.lattice_system).toBe("Hexagonal");
  });

  it("throws on invalid number", () => {
    expect(() => layer_group_type(0)).toThrow();
  });
});
