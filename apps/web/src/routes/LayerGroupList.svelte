<script lang="ts">
  import { push, link } from 'svelte-spa-router'
  import { getAllLayerGroups, filterRows } from '../lib/catalogIndex'
  import { CRYSTAL_SYSTEMS, type CrystalSystem } from '../lib/catalog'
  import { formatErr } from '../lib/wasm'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  let query = $state('')
  let arithSymbol = $state<string>('All')
  let geomClass = $state<string>('All')
  let system = $state<CrystalSystem | 'All'>('All')
  let lattice = $state<string>('All')

  const data = getAllLayerGroups()

  function onRowKey(e: KeyboardEvent, n: number) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault()
      push(`/lg/${n}`)
    }
  }

  function uniqueSorted(values: Iterable<string>): string[] {
    return Array.from(new Set(values)).sort()
  }
</script>

<header class="page-header">
  <div>
    <div class="eyebrow">Layer groups</div>
    <h1 class="text-2xl font-semibold">All 80 layer-group types</h1>
    <p class="text-sm text-slate-600 dark:text-slate-400">
      Search by number, HM symbol, arithmetic class, lattice system, ...
    </p>
  </div>
</header>

{#await data}
  <LoadingDots />
{:then rows}
  {@const arithOptions = uniqueSorted(rows.map((r) => r.arithmetic_symbol))}
  {@const geomOptions = uniqueSorted(rows.map((r) => r.geometric_crystal_class))}
  {@const latticeOptions = uniqueSorted(rows.map((r) => r.lattice_system))}
  {@const systemsInUse = new Set(rows.map((r) => r.crystal_system))}
  {@const filtered = filterRows(rows, query).filter(
    (r) =>
      (arithSymbol === 'All' || r.arithmetic_symbol === arithSymbol) &&
      (geomClass === 'All' || r.geometric_crystal_class === geomClass) &&
      (system === 'All' || r.crystal_system === system) &&
      (lattice === 'All' || r.lattice_system === lattice)
  )}

  <section class="space-y-4">
    <div class="flex flex-wrap items-center gap-3">
      <label class="flex-1 min-w-[16rem]">
        <span class="sr-only">Search</span>
        <input
          type="search"
          placeholder="Search... e.g. 'p4mm', '15', 'oblique'"
          class="search-input"
          bind:value={query}
        />
      </label>
      <span class="text-xs text-slate-500">{filtered.length} / {rows.length}</span>
    </div>

    <div class="flex flex-wrap items-center gap-x-4 gap-y-2 text-sm">
      <label class="flex items-center gap-2">
        <span class="text-slate-500">Arithmetic crystal class:</span>
        <select
          class="filter-select font-mono"
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
          class="filter-select font-mono"
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
          class="filter-select"
          bind:value={system}
        >
          <option value="All">All</option>
          {#each CRYSTAL_SYSTEMS as s}
            {#if systemsInUse.has(s)}
              <option value={s}>{s}</option>
            {/if}
          {/each}
        </select>
      </label>
      <label class="flex items-center gap-2">
        <span class="text-slate-500">Lattice system:</span>
        <select
          class="filter-select"
          bind:value={lattice}
        >
          <option value="All">All</option>
          {#each latticeOptions as l}
            <option value={l}>{l}</option>
          {/each}
        </select>
      </label>
      {#if arithSymbol !== 'All' || geomClass !== 'All' || system !== 'All' || lattice !== 'All'}
        <button
          type="button"
          class="link-button"
          onclick={() => {
            arithSymbol = 'All'
            geomClass = 'All'
            system = 'All'
            lattice = 'All'
          }}
        >
          Reset filters
        </button>
      {/if}
    </div>

    <div class="table-shell">
      <table class="min-w-full text-sm">
        <thead class="table-head">
          <tr>
            <th class="px-3 py-2 text-left">#</th>
            <th class="px-3 py-2 text-left">Short Hermann-Mauguin symbol</th>
            <th class="px-3 py-2 text-left">Full Hermann-Mauguin symbol</th>
            <th class="px-3 py-2 text-left">Arithmetic crystal class</th>
            <th class="px-3 py-2 text-left">Geometric crystal class</th>
            <th class="px-3 py-2 text-left">Crystal system</th>
            <th class="px-3 py-2 text-left">Lattice system</th>
          </tr>
        </thead>
        <tbody>
          {#each filtered as r (r.number)}
            <tr
              class="table-row-link"
              tabindex="0"
              role="link"
              onclick={() => push(`/lg/${r.number}`)}
              onkeydown={(e) => onRowKey(e, r.number)}
            >
              <td class="px-3 py-1.5 font-mono">{r.number}</td>
              <td class="px-3 py-1.5 font-mono">
                <a use:link href={`/lg/${r.number}`} class="hover:underline">{r.hm_short}</a>
              </td>
              <td class="px-3 py-1.5 font-mono text-slate-600 dark:text-slate-400">{r.hm_full}</td>
              <td class="px-3 py-1.5 font-mono">{r.arithmetic_symbol}</td>
              <td class="px-3 py-1.5 font-mono">{r.geometric_crystal_class}</td>
              <td class="px-3 py-1.5">{r.crystal_system}</td>
              <td class="px-3 py-1.5">{r.lattice_system}</td>
            </tr>
          {/each}
          {#if filtered.length === 0}
            <tr>
              <td colspan="7" class="px-3 py-6 text-center text-sm text-slate-500">
                No layer groups match the current filter.
              </td>
            </tr>
          {/if}
        </tbody>
      </table>
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load layer-group catalog: ${formatErr(err)}`} />
{/await}
