import { describe, expect, it } from "vitest";
import {
  magnetic_hall_symbol_entry,
  magnetic_operations_from_uni_number,
  magnetic_space_group_type,
} from "../pkg/moyo_wasm.js";

describe("magnetic_operations_from_uni_number", () => {
  it("UNI 1242 (R_I3, BNS 146.12) returns 18 conventional / 6 primitive", () => {
    const ops_conv = magnetic_operations_from_uni_number(1242, false);
    expect(ops_conv.length).toBe(18);
    const ops_prim = magnetic_operations_from_uni_number(1242, true);
    expect(ops_prim.length).toBe(6);
  });

  it("UNI 1 (P 1, BNS 1.1) returns a single identity operation", () => {
    const ops = magnetic_operations_from_uni_number(1, false);
    expect(ops.length).toBe(1);
    expect(ops[0].time_reversal).toBe(false);
  });

  it("UNI 2 (P 1 1', BNS 1.2) returns 2 operations (identity + time reversal)", () => {
    const ops = magnetic_operations_from_uni_number(2, false);
    expect(ops.length).toBe(2);
    const timeReversals = ops.map((op) => op.time_reversal).sort();
    expect(timeReversals).toEqual([false, true]);
  });

  it("throws on invalid uni_number", () => {
    expect(() => magnetic_operations_from_uni_number(0, false)).toThrow();
    expect(() => magnetic_operations_from_uni_number(1652, false)).toThrow();
  });
});

describe("magnetic_hall_symbol_entry", () => {
  it("uni_number 1 -> 'P 1'", () => {
    const entry = magnetic_hall_symbol_entry(1);
    expect(entry.uni_number).toBe(1);
    expect(entry.magnetic_hall_symbol).toBe("P 1");
  });

  it("uni_number 2 -> 'P 1 1''", () => {
    const entry = magnetic_hall_symbol_entry(2);
    expect(entry.uni_number).toBe(2);
    expect(entry.magnetic_hall_symbol).toBe("P 1 1'");
  });

  it("throws on invalid uni_number", () => {
    expect(() => magnetic_hall_symbol_entry(0)).toThrow();
    expect(() => magnetic_hall_symbol_entry(1652)).toThrow();
  });
});

describe("magnetic_space_group_type", () => {
  it("uni_number 1 -> BNS 1.1, ITA 1, Type 1", () => {
    const msgt = magnetic_space_group_type(1);
    expect(msgt.uni_number).toBe(1);
    expect(msgt.litvin_number).toBe(1);
    expect(msgt.bns_number).toBe("1.1");
    expect(msgt.og_number).toBe("1.1.1");
    expect(msgt.number).toBe(1);
    expect(msgt.construct_type).toBe(1);
  });

  it("uni_number 1262 -> BNS 151.32, Type 4", () => {
    const msgt = magnetic_space_group_type(1262);
    expect(msgt.bns_number).toBe("151.32");
    expect(msgt.og_number).toBe("153.4.1270");
    expect(msgt.construct_type).toBe(4);
  });

  it("throws on invalid uni_number", () => {
    expect(() => magnetic_space_group_type(0)).toThrow();
    expect(() => magnetic_space_group_type(1652)).toThrow();
  });
});
