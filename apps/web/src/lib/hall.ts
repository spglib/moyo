import { getMoyo } from './wasm'

type HallByNumber = Map<number, number[]>

let cached: Promise<HallByNumber> | null = null

/** Group every Hall number (1..530) by ITA space-group number. */
export function getHallNumbersBySpaceGroup(): Promise<HallByNumber> {
  if (!cached) {
    cached = getMoyo().then((m) => {
      const map: HallByNumber = new Map()
      for (let h = 1; h <= 530; h++) {
        const entry = m.hall_symbol_entry(h)
        const list = map.get(entry.number) ?? []
        list.push(entry.hall_number)
        map.set(entry.number, list)
      }
      return map
    })
  }
  return cached
}

let cachedLayer: Promise<HallByNumber> | null = null

export function getLayerHallNumbersByGroup(): Promise<HallByNumber> {
  if (!cachedLayer) {
    cachedLayer = getMoyo().then((m) => {
      const map: HallByNumber = new Map()
      for (let h = 1; h <= 116; h++) {
        const entry = m.layer_hall_symbol_entry(h)
        const list = map.get(entry.number) ?? []
        list.push(entry.hall_number)
        map.set(entry.number, list)
      }
      return map
    })
  }
  return cachedLayer
}
