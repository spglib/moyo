type Vec3 = readonly [number, number, number]
type Mat9 = readonly [number, number, number, number, number, number, number, number, number]

const COMMON_DENOMS = [1, 2, 3, 4, 6, 8, 12]
const EPS = 1e-6

/** Render a symmetry operation as a comma-separated coordinate triplet
 *  (Jones-faithful notation), e.g. `-x+1/2, y, -z+1/2`. When `timeReversal`
 *  is given, a trailing `, +1` or `, -1` is appended for magnetic ops.
 *
 *  `rotation` is the nalgebra-flattened 3x3 in column-major order, so the
 *  element at row `i`, column `j` is `rotation[j * 3 + i]`. */
export function formatOperationXyz(
  rotation: Mat9,
  translation: Vec3,
  timeReversal?: boolean
): string {
  const vars = ['x', 'y', 'z']
  const triplet = [0, 1, 2]
    .map((i) => {
      const parts: string[] = []
      for (let j = 0; j < 3; j++) {
        const c = rotation[j * 3 + i]
        if (Math.abs(c) < EPS) continue
        const sign = c < 0 ? '-' : parts.length > 0 ? '+' : ''
        const mag = Math.abs(c)
        const mantissa = Math.abs(mag - 1) < EPS ? '' : formatFraction(mag)
        parts.push(`${sign}${mantissa}${vars[j]}`)
      }
      const t = translation[i]
      if (Math.abs(t) > EPS) {
        const sign = t < 0 ? '-' : parts.length > 0 ? '+' : ''
        parts.push(`${sign}${formatFraction(Math.abs(t))}`)
      }
      return parts.length > 0 ? parts.join('') : '0'
    })
    .join(', ')
  if (timeReversal === undefined) return triplet
  return `${triplet}, ${timeReversal ? '-1' : '+1'}`
}

export function formatFraction(x: number): string {
  if (!isFinite(x)) return String(x)
  const sign = x < 0 ? '-' : ''
  const a = Math.abs(x)
  const whole = Math.floor(a + EPS)
  const frac = a - whole
  if (frac < EPS) return `${sign}${whole}`
  for (const d of COMMON_DENOMS) {
    const n = Math.round(frac * d)
    if (Math.abs(frac - n / d) < EPS && n > 0 && n < d) {
      return whole === 0 ? `${sign}${n}/${d}` : `${sign}${whole} ${n}/${d}`
    }
  }
  return `${sign}${a.toFixed(4)}`
}

function fmtInt(x: number): string {
  if (Math.abs(x - Math.round(x)) < EPS) return String(Math.round(x))
  return x.toFixed(3)
}

export function constructTypeLabel(t: number): string {
  switch (t) {
    case 1:
      return 'Type I (colorless)'
    case 2:
      return 'Type II (grey)'
    case 3:
      return 'Type III (BW, equi-translation)'
    case 4:
      return 'Type IV (BW, equi-class)'
    default:
      return `Type ${t}`
  }
}
