import type { Moyo } from './wasm'
import { getMoyo } from './wasm'

export type HallByNumber = Map<number, number[]>

function buildIndex(
  count: number,
  fetch: (m: Moyo, h: number) => { number: number; hall_number: number }
): Promise<HallByNumber> {
  return getMoyo().then((m) => {
    const map: HallByNumber = new Map()
    for (let h = 1; h <= count; h++) {
      const entry = fetch(m, h)
      const list = map.get(entry.number) ?? []
      list.push(entry.hall_number)
      map.set(entry.number, list)
    }
    return map
  })
}

let cachedSg: Promise<HallByNumber> | null = null
let cachedLg: Promise<HallByNumber> | null = null

/** Group every Hall number (1..530) by ITA space-group number. */
export function getHallNumbersBySpaceGroup(): Promise<HallByNumber> {
  if (!cachedSg) cachedSg = buildIndex(530, (m, h) => m.hall_symbol_entry(h))
  return cachedSg
}

/** Group every layer Hall number (1..116) by layer-group number. */
export function getLayerHallNumbersByGroup(): Promise<HallByNumber> {
  if (!cachedLg) cachedLg = buildIndex(116, (m, h) => m.layer_hall_symbol_entry(h))
  return cachedLg
}
