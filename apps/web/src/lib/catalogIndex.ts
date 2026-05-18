import type { MoyoSpaceGroupType } from '@spglib/moyo-wasm'
import { getMoyo } from './wasm'
import { SPACE_GROUP_COUNT } from './catalog'

export interface SpaceGroupRow extends MoyoSpaceGroupType {
  searchText: string
}

let cached: Promise<SpaceGroupRow[]> | null = null

/** Eagerly load every 1..230 space-group-type entry, indexed for filtering. */
export function getAllSpaceGroups(): Promise<SpaceGroupRow[]> {
  if (!cached) {
    cached = getMoyo().then((m) => {
      const rows: SpaceGroupRow[] = []
      for (let n = 1; n <= SPACE_GROUP_COUNT; n++) {
        const t = m.space_group_type(n)
        rows.push({ ...t, searchText: searchText(t) })
      }
      return rows
    })
  }
  return cached
}

function searchText(t: MoyoSpaceGroupType): string {
  return [
    t.number,
    t.hm_short,
    t.hm_full,
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

export function filterSpaceGroups(rows: SpaceGroupRow[], query: string): SpaceGroupRow[] {
  const q = query.trim().toLowerCase()
  if (!q) return rows
  const tokens = q.split(/\s+/)
  return rows.filter((r) => tokens.every((tok) => r.searchText.includes(tok)))
}
