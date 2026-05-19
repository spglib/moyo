import type {
  MoyoSpaceGroupType,
  MoyoLayerGroupType,
  MoyoMagneticSpaceGroupType,
  MoyoMagneticHallSymbolEntry,
} from '@spglib/moyo-wasm'
import { getMoyo } from './wasm'
import { SPACE_GROUP_COUNT, LAYER_GROUP_COUNT, MAGNETIC_SG_COUNT } from './catalog'
import { constructTypeLabel } from './format'

export interface SpaceGroupRow extends MoyoSpaceGroupType {
  searchText: string
}

export interface LayerGroupRow extends MoyoLayerGroupType {
  searchText: string
}

export interface MagneticSpaceGroupRow extends MoyoMagneticSpaceGroupType {
  magnetic_hall_symbol: string
  parent_hm_short: string
  construct_label: string
  searchText: string
}

let cachedSg: Promise<SpaceGroupRow[]> | null = null
let cachedLg: Promise<LayerGroupRow[]> | null = null
let cachedMsg: Promise<MagneticSpaceGroupRow[]> | null = null

export function getAllSpaceGroups(): Promise<SpaceGroupRow[]> {
  if (!cachedSg) {
    cachedSg = getMoyo().then((m) => {
      const rows: SpaceGroupRow[] = []
      for (let n = 1; n <= SPACE_GROUP_COUNT; n++) {
        const t = m.space_group_type(n)
        rows.push({ ...t, searchText: spaceGroupSearchText(t) })
      }
      return rows
    })
  }
  return cachedSg
}

export function getAllLayerGroups(): Promise<LayerGroupRow[]> {
  if (!cachedLg) {
    cachedLg = getMoyo().then((m) => {
      const rows: LayerGroupRow[] = []
      for (let n = 1; n <= LAYER_GROUP_COUNT; n++) {
        const t = m.layer_group_type(n)
        rows.push({ ...t, searchText: layerGroupSearchText(t) })
      }
      return rows
    })
  }
  return cachedLg
}

export function getAllMagneticSpaceGroups(): Promise<MagneticSpaceGroupRow[]> {
  if (!cachedMsg) {
    cachedMsg = getMoyo().then((m) => {
      const parents = new Map<number, string>()
      const rows: MagneticSpaceGroupRow[] = []
      for (let n = 1; n <= MAGNETIC_SG_COUNT; n++) {
        const t = m.magnetic_space_group_type(n)
        const hall = m.magnetic_hall_symbol_entry(n)
        if (!parents.has(t.number)) {
          parents.set(t.number, m.space_group_type(t.number).hm_short)
        }
        const parent_hm_short = parents.get(t.number) ?? ''
        const construct_label = constructTypeLabel(t.construct_type)
        rows.push({
          ...t,
          magnetic_hall_symbol: hall.magnetic_hall_symbol,
          parent_hm_short,
          construct_label,
          searchText: magneticSpaceGroupSearchText(t, hall, parent_hm_short, construct_label),
        })
      }
      return rows
    })
  }
  return cachedMsg
}

/** Stripped form alongside the original, so a search for either
 *  `Fd-3m` or `F d -3 m` matches `hm_short: "F d -3 m"`. */
function withCompact(s: string): string {
  const compact = s.replace(/\s+/g, '')
  return compact === s ? s : `${s} ${compact}`
}

function spaceGroupSearchText(t: MoyoSpaceGroupType): string {
  return [
    t.number,
    withCompact(t.hm_short),
    withCompact(t.hm_full),
    t.arithmetic_number,
    t.arithmetic_symbol,
    t.geometric_crystal_class,
    t.crystal_system,
    t.crystal_family,
    t.lattice_system,
    t.bravais_class,
  ]
    .join(' ')
    .toLowerCase()
}

function layerGroupSearchText(t: MoyoLayerGroupType): string {
  return [
    t.number,
    withCompact(t.hm_short),
    withCompact(t.hm_full),
    t.arithmetic_number,
    t.arithmetic_symbol,
    t.geometric_crystal_class,
    t.lattice_system,
    t.bravais_class,
  ]
    .join(' ')
    .toLowerCase()
}

function magneticSpaceGroupSearchText(
  t: MoyoMagneticSpaceGroupType,
  hall: MoyoMagneticHallSymbolEntry,
  parent_hm_short: string,
  construct_label: string
): string {
  return [
    t.uni_number,
    t.litvin_number,
    t.bns_number,
    t.og_number,
    t.number,
    withCompact(parent_hm_short),
    withCompact(hall.magnetic_hall_symbol),
    construct_label,
  ]
    .join(' ')
    .toLowerCase()
}

export function filterRows<T extends { searchText: string }>(rows: T[], query: string): T[] {
  const q = query.trim().toLowerCase()
  if (!q) return rows
  const tokens = q.split(/\s+/)
  return rows.filter((r) => tokens.every((tok) => r.searchText.includes(tok)))
}
