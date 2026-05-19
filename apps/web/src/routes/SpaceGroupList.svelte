<script lang="ts">
  import { push, link } from 'svelte-spa-router'
  import { getAllSpaceGroups, filterRows } from '../lib/catalogIndex'
  import { CRYSTAL_SYSTEMS, type CrystalSystem } from '../lib/catalog'
  import { formatErr } from '../lib/wasm'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  let query = $state('')
  let system = $state<CrystalSystem | 'All'>('All')
  let geomClass = $state<string>('All')
  let arithSymbol = $state<string>('All')

  const data = getAllSpaceGroups()

  function onRowKey(e: KeyboardEvent, n: number) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault()
      push(`/sg/${n}`)
    }
  }

  function uniqueSorted(values: Iterable<string>): string[] {
    return Array.from(new Set(values)).sort()
  }
</script>

<header
  class="flex flex-wrap items-end justify-between gap-4 py-2 border-b border-slate-200 dark:border-slate-800 mb-4"
>
  <div>
    <div class="text-xs uppercase tracking-wide text-slate-500">Space groups</div>
    <h1 class="text-2xl font-semibold">All 230 ITA space-group types</h1>
    <p class="text-sm text-slate-600 dark:text-slate-400">
      Search by number, HM symbol, arithmetic class, crystal system, Bravais class, ...
    </p>
  </div>
</header>

{#await data}
  <LoadingDots />
{:then rows}
  {@const geomOptions = uniqueSorted(rows.map((r) => r.geometric_crystal_class))}
  {@const arithOptions = uniqueSorted(rows.map((r) => r.arithmetic_symbol))}
  {@const filtered = filterRows(rows, query).filter(
    (r) =>
      (system === 'All' || r.crystal_system === system) &&
      (geomClass === 'All' || r.geometric_crystal_class === geomClass) &&
      (arithSymbol === 'All' || r.arithmetic_symbol === arithSymbol)
  )}

  <section class="space-y-4">
    <div class="flex flex-wrap items-center gap-3">
      <label class="flex-1 min-w-[16rem]">
        <span class="sr-only">Search</span>
        <input
          type="search"
          placeholder="Search... e.g. 'Fd-3m', '227', 'cubic m-3m'"
          class="w-full rounded border border-slate-300 dark:border-slate-700 bg-transparent px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-sky-500"
          bind:value={query}
        />
      </label>
      <span class="text-xs text-slate-500">{filtered.length} / {rows.length}</span>
    </div>

    <div class="flex flex-wrap items-center gap-x-4 gap-y-2 text-sm">
      <label class="flex items-center gap-2">
        <span class="text-slate-500">Arithmetic crystal class:</span>
        <select
          class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1 font-mono"
          bind:value={arithSymbol}
        >
          <option value="All">All</option>
          {#each arithOptions as a}
            <option value={a}>{a}</option>
          {/each}
        </select>
      </label>
      <label class="flex items-center gap-2">
        <span class="text-slate-500">Geometric crystal class:</span>
        <select
          class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1 font-mono"
          bind:value={geomClass}
        >
          <option value="All">All</option>
          {#each geomOptions as g}
            <option value={g}>{g}</option>
          {/each}
        </select>
      </label>
      <label class="flex items-center gap-2">
        <span class="text-slate-500">Crystal system:</span>
        <select
          class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1"
          bind:value={system}
        >
          <option value="All">All</option>
          {#each CRYSTAL_SYSTEMS as s}
            <option value={s}>{s}</option>
          {/each}
        </select>
      </label>
      {#if system !== 'All' || geomClass !== 'All' || arithSymbol !== 'All'}
        <button
          type="button"
          class="text-xs text-slate-500 hover:underline"
          onclick={() => {
            system = 'All'
            geomClass = 'All'
            arithSymbol = 'All'
          }}
        >
          Reset filters
        </button>
      {/if}
    </div>

    <div class="overflow-x-auto rounded border border-slate-200 dark:border-slate-800">
      <table class="min-w-full text-sm">
        <thead class="bg-slate-50 dark:bg-slate-900 text-xs uppercase tracking-wide">
          <tr>
            <th class="px-3 py-2 text-left">#</th>
            <th class="px-3 py-2 text-left">Short Hermann-Mauguin symbol</th>
            <th class="px-3 py-2 text-left">Full Hermann-Mauguin symbol</th>
            <th class="px-3 py-2 text-left">Arithmetic crystal class</th>
            <th class="px-3 py-2 text-left">Geometric crystal class</th>
            <th class="px-3 py-2 text-left">Crystal system</th>
          </tr>
        </thead>
        <tbody>
          {#each filtered as r (r.number)}
            <tr
              class="border-t border-slate-100 dark:border-slate-800 hover:bg-slate-50 dark:hover:bg-slate-900/60 cursor-pointer"
              tabindex="0"
              role="link"
              onclick={() => push(`/sg/${r.number}`)}
              onkeydown={(e) => onRowKey(e, r.number)}
            >
              <td class="px-3 py-1.5 font-mono">{r.number}</td>
              <td class="px-3 py-1.5 font-mono">
                <a use:link href={`/sg/${r.number}`} class="hover:underline">{r.hm_short}</a>
              </td>
              <td class="px-3 py-1.5 font-mono text-slate-600 dark:text-slate-400">{r.hm_full}</td>
              <td class="px-3 py-1.5 font-mono">{r.arithmetic_symbol}</td>
              <td class="px-3 py-1.5 font-mono">{r.geometric_crystal_class}</td>
              <td class="px-3 py-1.5">{r.crystal_system}</td>
            </tr>
          {/each}
          {#if filtered.length === 0}
            <tr>
              <td colspan="6" class="px-3 py-6 text-center text-sm text-slate-500">
                No space groups match the current filter.
              </td>
            </tr>
          {/if}
        </tbody>
      </table>
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load space-group catalog: ${formatErr(err)}`} />
{/await}
