import { describe, expect, it } from "vitest";
import {
  arithmetic_crystal_class,
  hall_symbol_entry,
  operations_from_number,
  space_group_type,
} from "../pkg/moyo_wasm.js";

describe("operations_from_number", () => {
  it("Ia-3d (#230, Standard) returns 96 conventional operations", () => {
    const ops = operations_from_number(230, { type: "Standard" }, false);
    expect(ops.length).toBe(96);
  });

  it("Ia-3d (#230, Standard, primitive) returns 48 operations", () => {
    const ops = operations_from_number(230, { type: "Standard" }, true);
    expect(ops.length).toBe(48);
  });

  it("Pm-3m (#221, Standard) returns 48 operations", () => {
    const ops = operations_from_number(221, { type: "Standard" }, false);
    expect(ops.length).toBe(48);
  });

  it("P1 (#1, Standard) returns 1 operation", () => {
    const ops = operations_from_number(1, { type: "Standard" }, false);
    expect(ops.length).toBe(1);
  });

  it("Spglib setting also works (#15 C2/c)", () => {
    const ops = operations_from_number(15, { type: "Spglib" }, false);
    expect(ops.length).toBe(8);
  });

  it("operations have 3x3 rotation and length-3 translation", () => {
    const ops = operations_from_number(1, { type: "Standard" }, false);
    expect(ops[0].rotation.length).toBe(9);
    expect(ops[0].translation.length).toBe(3);
  });

  it("throws on invalid number", () => {
    expect(() =>
      operations_from_number(0, { type: "Standard" }, false),
    ).toThrow();
  });
});

describe("hall_symbol_entry", () => {
  it("returns Fd-3c entry for hall_number 528 (#228 origin choice 2)", () => {
    const entry = hall_symbol_entry(528);
    expect(entry.hall_number).toBe(528);
    expect(entry.number).toBe(228);
    expect(entry.setting).toBe("2");
    expect(entry.hm_short).toBe("F d -3 c");
    expect(entry.centering).toBe("F");
  });

  it("throws on invalid hall_number", () => {
    expect(() => hall_symbol_entry(0)).toThrow();
  });
});

describe("arithmetic_crystal_class", () => {
  it("returns 3m1P for arithmetic_number 45", () => {
    const c = arithmetic_crystal_class(45);
    expect(c.arithmetic_number).toBe(45);
    expect(c.symbol).toBe("3m1P");
    expect(c.geometric_crystal_class).toBe("3m");
    expect(c.bravais_class).toBe("hP");
  });

  it("throws on invalid arithmetic_number", () => {
    expect(() => arithmetic_crystal_class(0)).toThrow();
  });
});

describe("space_group_type", () => {
  it("Pm-3m (#221) reports cubic classifications", () => {
    const t = space_group_type(221);
    expect(t.number).toBe(221);
    expect(t.hm_short).toBe("P m -3 m");
    expect(t.arithmetic_number).toBe(71);
    expect(t.arithmetic_symbol).toBe("m-3mP");
    expect(t.geometric_crystal_class).toBe("m-3m");
    expect(t.crystal_system).toBe("Cubic");
    expect(t.bravais_class).toBe("cP");
    expect(t.lattice_system).toBe("Cubic");
    expect(t.crystal_family).toBe("Cubic");
  });

  it("R-3c (#167) reports trigonal/rhombohedral classifications", () => {
    const t = space_group_type(167);
    expect(t.number).toBe(167);
    expect(t.hm_short).toBe("R -3 c");
    expect(t.arithmetic_number).toBe(50);
    expect(t.arithmetic_symbol).toBe("-3mR");
    expect(t.geometric_crystal_class).toBe("-3m");
    expect(t.crystal_system).toBe("Trigonal");
    expect(t.bravais_class).toBe("hR");
    expect(t.lattice_system).toBe("Rhombohedral");
    expect(t.crystal_family).toBe("Hexagonal");
  });

  it("throws on invalid number", () => {
    expect(() => space_group_type(0)).toThrow();
  });
});
