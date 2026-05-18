import { describe, it, expect } from 'vitest'
import { formatFraction, formatRotationRow, constructTypeLabel } from '../src/lib/format'

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

describe('constructTypeLabel', () => {
  it('labels magnetic construct types', () => {
    expect(constructTypeLabel(1)).toMatch(/Fedorov/)
    expect(constructTypeLabel(2)).toMatch(/grey/)
    expect(constructTypeLabel(3)).toMatch(/BW/)
    expect(constructTypeLabel(4)).toMatch(/BW/)
  })
})
