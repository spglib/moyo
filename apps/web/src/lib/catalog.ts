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

/** Crystal system of each geometric crystal class symbol (as emitted by
 *  moyo, see moyo/src/data/classification.rs CrystalSystem). */
const CRYSTAL_SYSTEM_OF_POINT_GROUP: Record<string, CrystalSystem> = {
  '1': 'Triclinic',
  '-1': 'Triclinic',
  '2': 'Monoclinic',
  m: 'Monoclinic',
  '2/m': 'Monoclinic',
  '222': 'Orthorhombic',
  mm2: 'Orthorhombic',
  mmm: 'Orthorhombic',
  '4': 'Tetragonal',
  '-4': 'Tetragonal',
  '4/m': 'Tetragonal',
  '422': 'Tetragonal',
  '4mm': 'Tetragonal',
  '-42m': 'Tetragonal',
  '4/mmm': 'Tetragonal',
  '3': 'Trigonal',
  '-3': 'Trigonal',
  '32': 'Trigonal',
  '3m': 'Trigonal',
  '-3m': 'Trigonal',
  '6': 'Hexagonal',
  '-6': 'Hexagonal',
  '6/m': 'Hexagonal',
  '622': 'Hexagonal',
  '6mm': 'Hexagonal',
  '-6m2': 'Hexagonal',
  '6/mmm': 'Hexagonal',
  '23': 'Cubic',
  'm-3': 'Cubic',
  '432': 'Cubic',
  '-43m': 'Cubic',
  'm-3m': 'Cubic',
}

export function crystalSystemOfPointGroup(pointGroup: string): CrystalSystem | null {
  return CRYSTAL_SYSTEM_OF_POINT_GROUP[pointGroup] ?? null
}

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
