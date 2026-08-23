export const SPACE_GROUP_COUNT = 230
export const LAYER_GROUP_COUNT = 80
export const MAGNETIC_SG_COUNT = 1651

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
