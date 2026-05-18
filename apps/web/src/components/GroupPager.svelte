<script lang="ts">
  import { push } from 'svelte-spa-router'

  interface Props {
    value: number
    min: number
    max: number
    basePath: string
    label?: string
    query?: string
  }

  let { value, min, max, basePath, label = '#', query = '' }: Props = $props()

  let input = $state('')
  $effect(() => {
    input = String(value)
  })

  function go(n: number) {
    const clamped = Math.min(Math.max(n, min), max)
    const url = `${basePath}/${clamped}${query ? `?${query}` : ''}`
    push(url)
  }

  function submit(e: Event) {
    e.preventDefault()
    const n = Number(input)
    if (Number.isFinite(n)) go(Math.trunc(n))
  }
</script>

<form onsubmit={submit} class="flex items-center gap-2 text-sm">
  <button
    type="button"
    class="rounded border border-slate-300 dark:border-slate-700 px-2 py-1 hover:bg-slate-100 dark:hover:bg-slate-800 disabled:opacity-50"
    disabled={value <= min}
    onclick={() => go(value - 1)}
    aria-label="Previous"
  >
    &larr;
  </button>
  <label class="flex items-center gap-1">
    <span class="text-slate-500">{label}</span>
    <input
      class="w-20 rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1 font-mono"
      type="number"
      min={min}
      max={max}
      bind:value={input}
    />
    <span class="text-slate-500">/ {max}</span>
  </label>
  <button
    type="button"
    class="rounded border border-slate-300 dark:border-slate-700 px-2 py-1 hover:bg-slate-100 dark:hover:bg-slate-800 disabled:opacity-50"
    disabled={value >= max}
    onclick={() => go(value + 1)}
    aria-label="Next"
  >
    &rarr;
  </button>
</form>
