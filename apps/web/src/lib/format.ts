type Vec3 = readonly [number, number, number]
type Mat9 = readonly [number, number, number, number, number, number, number, number, number]

export function formatRotationRow(m: Mat9, row: 0 | 1 | 2): string {
  const i = row * 3
  return `[${fmtInt(m[i])}, ${fmtInt(m[i + 1])}, ${fmtInt(m[i + 2])}]`
}

export function formatRotation(m: Mat9): string {
  return [formatRotationRow(m, 0), formatRotationRow(m, 1), formatRotationRow(m, 2)].join(' / ')
}

export function formatTranslation(t: Vec3): string {
  return `(${formatFraction(t[0])}, ${formatFraction(t[1])}, ${formatFraction(t[2])})`
}

const COMMON_DENOMS = [1, 2, 3, 4, 6, 8, 12]
const EPS = 1e-6

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
      return 'Type I (Fedorov)'
    case 2:
      return 'Type II (grey)'
    case 3:
      return 'Type III (BW, equi-class)'
    case 4:
      return 'Type IV (BW, equi-translation)'
    default:
      return `Type ${t}`
  }
}
