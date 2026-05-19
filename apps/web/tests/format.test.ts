import { describe, it, expect } from 'vitest'
import {
  formatFraction,
  formatRotationRow,
  formatOperationXyz,
  constructTypeLabel,
} from '../src/lib/format'

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

describe('formatRotationRow', () => {
  it('formats a row of a 3x3 matrix', () => {
    const m = [1, 0, 0, 0, -1, 0, 0, 0, 1] as const
    expect(formatRotationRow(m, 0)).toBe('[1, 0, 0]')
    expect(formatRotationRow(m, 1)).toBe('[0, -1, 0]')
    expect(formatRotationRow(m, 2)).toBe('[0, 0, 1]')
  })
})

describe('formatOperationXyz', () => {
  const ident = [1, 0, 0, 0, 1, 0, 0, 0, 1] as const

  it('renders the identity', () => {
    expect(formatOperationXyz(ident, [0, 0, 0])).toBe('x, y, z')
  })

  it('renders inversion', () => {
    const inv = [-1, 0, 0, 0, -1, 0, 0, 0, -1] as const
    expect(formatOperationXyz(inv, [0, 0, 0])).toBe('-x, -y, -z')
  })

  it('renders a screw axis with translation', () => {
    const r2 = [-1, 0, 0, 0, -1, 0, 0, 0, 1] as const
    expect(formatOperationXyz(r2, [0, 0, 0.5])).toBe('-x, -y, z+1/2')
  })

  it('renders mixed-sign rotation entries', () => {
    const r = [0, -1, 0, 1, -1, 0, 0, 0, 1] as const
    expect(formatOperationXyz(r, [0, 0, 0])).toBe('-y, x-y, z')
  })

  it('renders a negative translation', () => {
    expect(formatOperationXyz(ident, [-0.25, 0, 0])).toBe('x-1/4, y, z')
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
