/**
 * Parser for Hermann-Mauguin short symbols as written in the moyo database.
 *
 * Tokens are space-separated. Within a token:
 *   - `-X` denotes a rotoinversion: render X with an overline.
 *   - `_X` denotes a subscript: attach to the preceding glyph.
 *   - everything else is plain text (letters, digits, `/`).
 *
 * Examples:
 *   "F d -3 m"      -> F  d  3̄  m
 *   "P 6_3/m m c"   -> P  6₃/m  m  c
 *   "P 2_1"         -> P  2₁
 */

export type Segment =
  | { kind: 'plain'; text: string }
  | { kind: 'over'; text: string }
  | { kind: 'sub'; text: string }

export type Token = Segment[]

export function parseHmShort(symbol: string): Token[] {
  return symbol.trim().split(/\s+/).filter(Boolean).map(parseToken)
}

function parseToken(token: string): Segment[] {
  const out: Segment[] = []
  let i = 0
  while (i < token.length) {
    const c = token[i]
    if (c === '-' && i + 1 < token.length) {
      out.push({ kind: 'over', text: token[i + 1] })
      i += 2
    } else if (c === '_' && i + 1 < token.length) {
      let j = i + 1
      while (j < token.length && /[A-Za-z0-9]/.test(token[j])) j++
      out.push({ kind: 'sub', text: token.slice(i + 1, j) })
      i = j
    } else {
      let j = i + 1
      while (j < token.length && token[j] !== '-' && token[j] !== '_') j++
      out.push({ kind: 'plain', text: token.slice(i, j) })
      i = j
    }
  }
  return out
}
