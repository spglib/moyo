import { describe, it, expect } from 'vitest'
import { formatFraction, formatOperationXyz, constructTypeLabel } from '../src/lib/format'

/** Build a column-major 9-element slice from row-major arguments,
 *  matching what nalgebra::Matrix3 emits via .as_slice(). */
function colMajor(
  r00: number, r01: number, r02: number,
  r10: number, r11: number, r12: number,
  r20: number, r21: number, r22: number
): readonly [number, number, number, number, number, number, number, number, number] {
  return [r00, r10, r20, r01, r11, r21, r02, r12, r22] as const
}

describe('formatFraction', () => {
  it('renders integer translations', () => {
    expect(formatFraction(0)).toBe('0')
    expect(formatFraction(1)).toBe('1')
    expect(formatFraction(-2)).toBe('-2')
  })

  it('renders common fractions', () => {
    expect(formatFraction(0.5)).toBe('1/2')
    expect(formatFraction(1 / 3)).toBe('1/3')
    expect(formatFraction(2 / 3)).toBe('2/3')
    expect(formatFraction(0.25)).toBe('1/4')
  })

  it('renders signed fractions', () => {
    expect(formatFraction(-0.5)).toBe('-1/2')
  })

  it('falls back to decimal for unusual values', () => {
    expect(formatFraction(0.13)).toMatch(/0\.13/)
  })
})

describe('formatOperationXyz', () => {
  const ident = colMajor(1, 0, 0, 0, 1, 0, 0, 0, 1)

  it('renders the identity', () => {
    expect(formatOperationXyz(ident, [0, 0, 0])).toBe('x, y, z')
  })

  it('renders inversion', () => {
    const inv = colMajor(-1, 0, 0, 0, -1, 0, 0, 0, -1)
    expect(formatOperationXyz(inv, [0, 0, 0])).toBe('-x, -y, -z')
  })

  it('renders a screw axis with translation', () => {
    const r2 = colMajor(-1, 0, 0, 0, -1, 0, 0, 0, 1)
    expect(formatOperationXyz(r2, [0, 0, 0.5])).toBe('-x, -y, z+1/2')
  })

  it('renders the P3 three-fold rotation operations', () => {
    // R3 = [[0,-1,0],[1,-1,0],[0,0,1]]
    const r3 = colMajor(0, -1, 0, 1, -1, 0, 0, 0, 1)
    expect(formatOperationXyz(r3, [0, 0, 0])).toBe('-y, x-y, z')

    // R3^2 = [[-1,1,0],[-1,0,0],[0,0,1]]
    const r3sq = colMajor(-1, 1, 0, -1, 0, 0, 0, 0, 1)
    expect(formatOperationXyz(r3sq, [0, 0, 0])).toBe('-x+y, -x, z')
  })

  it('renders a negative translation', () => {
    expect(formatOperationXyz(ident, [-0.25, 0, 0])).toBe('x-1/4, y, z')
  })

  it('appends the time-reversal sign for magnetic operations', () => {
    expect(formatOperationXyz(ident, [0, 0, 0], false)).toBe('x, y, z, +1')
    expect(formatOperationXyz(ident, [0, 0, 0], true)).toBe('x, y, z, -1')
  })
})

describe('constructTypeLabel', () => {
  it('labels magnetic construct types', () => {
    expect(constructTypeLabel(1)).toMatch(/Fedorov/)
    expect(constructTypeLabel(2)).toMatch(/grey/)
    expect(constructTypeLabel(3)).toMatch(/BW/)
    expect(constructTypeLabel(4)).toMatch(/BW/)
  })
})
