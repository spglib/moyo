export const SPACE_GROUP_COUNT = 230
export const LAYER_GROUP_COUNT = 80
export const MAGNETIC_SG_COUNT = 1651

export const SPACE_GROUP_NUMBERS = range(1, SPACE_GROUP_COUNT)
export const LAYER_GROUP_NUMBERS = range(1, LAYER_GROUP_COUNT)
export const MAGNETIC_UNI_NUMBERS = range(1, MAGNETIC_SG_COUNT)

export const CRYSTAL_SYSTEMS = [
  'Triclinic',
  'Monoclinic',
  'Orthorhombic',
  'Tetragonal',
  'Trigonal',
  'Hexagonal',
  'Cubic',
] as const

export type CrystalSystem = (typeof CRYSTAL_SYSTEMS)[number]

export function clampInt(value: number, lo: number, hi: number): number {
  const n = Math.trunc(value)
  if (!Number.isFinite(n)) return lo
  if (n < lo) return lo
  if (n > hi) return hi
  return n
}

function range(lo: number, hi: number): number[] {
  const out = new Array(hi - lo + 1)
  for (let i = 0; i < out.length; i++) out[i] = lo + i
  return out
}
