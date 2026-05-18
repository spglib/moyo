export function parseQuery(qs: string | undefined): URLSearchParams {
  return new URLSearchParams(qs ?? '')
}

export function buildQuery(params: Record<string, string | number | boolean | null | undefined>): string {
  const out = new URLSearchParams()
  for (const [k, v] of Object.entries(params)) {
    if (v === null || v === undefined || v === '' || v === false) continue
    out.set(k, v === true ? '1' : String(v))
  }
  return out.toString()
}
