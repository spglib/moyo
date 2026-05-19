import init, * as moyo from '@spglib/moyo-wasm'
// Vite resolves `?url` imports to a hashed asset URL at build time.
import wasmUrl from '@spglib/moyo-wasm/moyo_wasm_bg.wasm?url'

export type Moyo = typeof moyo

let pending: Promise<Moyo> | null = null

export function getMoyo(): Promise<Moyo> {
  if (!pending) {
    pending = init({ module_or_path: wasmUrl }).then(() => moyo)
  }
  return pending
}

export function formatErr(err: unknown): string {
  if (typeof err === 'string') return err
  if (err instanceof Error) return err.message
  try {
    return JSON.stringify(err)
  } catch {
    return String(err)
  }
}
