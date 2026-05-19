import { describe, it, expect } from 'vitest'
import { parseHmShort } from '../src/lib/hmShort'

describe('parseHmShort', () => {
  it('parses a token with an overline', () => {
    expect(parseHmShort('F d -3 m')).toEqual([
      [{ kind: 'plain', text: 'F' }],
      [{ kind: 'plain', text: 'd' }],
      [{ kind: 'over', text: '3' }],
      [{ kind: 'plain', text: 'm' }],
    ])
  })

  it('parses a screw-axis subscript', () => {
    expect(parseHmShort('P 2_1')).toEqual([
      [{ kind: 'plain', text: 'P' }],
      [
        { kind: 'plain', text: '2' },
        { kind: 'sub', text: '1' },
      ],
    ])
  })

  it('parses a screw-axis + glide combo', () => {
    expect(parseHmShort('P 6_3/m m c')).toEqual([
      [{ kind: 'plain', text: 'P' }],
      [
        { kind: 'plain', text: '6' },
        { kind: 'sub', text: '3' },
        { kind: 'plain', text: '/m' },
      ],
      [{ kind: 'plain', text: 'm' }],
      [{ kind: 'plain', text: 'c' }],
    ])
  })

  it('parses an overline followed by a glide', () => {
    expect(parseHmShort('-3 1 m')).toEqual([
      [{ kind: 'over', text: '3' }],
      [{ kind: 'plain', text: '1' }],
      [{ kind: 'plain', text: 'm' }],
    ])
  })

  it('handles a trivial symbol', () => {
    expect(parseHmShort('P 1')).toEqual([
      [{ kind: 'plain', text: 'P' }],
      [{ kind: 'plain', text: '1' }],
    ])
  })
})
