<script lang="ts">
  import { push, link } from 'svelte-spa-router'
  import { getAllMagneticSpaceGroups, filterRows } from '../lib/catalogIndex'
  import { formatErr } from '../lib/wasm'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  let query = $state('')
  let constructType = $state<'All' | '1' | '2' | '3' | '4'>('All')

  const data = getAllMagneticSpaceGroups()

  function onRowKey(e: KeyboardEvent, n: number) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault()
      push(`/msg/${n}`)
    }
  }
</script>

<header
  class="flex flex-wrap items-end justify-between gap-4 py-2 border-b border-slate-200 dark:border-slate-800 mb-4"
>
  <div>
    <div class="text-xs uppercase tracking-wide text-slate-500">Magnetic space groups</div>
    <h1 class="text-2xl font-semibold">All 1651 magnetic space-group types</h1>
    <p class="text-sm text-slate-600 dark:text-slate-400">
      Search by UNI / Litvin / BNS / OG number, parent SG, magnetic Hall symbol, ...
    </p>
  </div>
</header>

{#await data}
  <LoadingDots />
{:then rows}
  {@const filtered = filterRows(rows, query).filter(
    (r) => constructType === 'All' || String(r.construct_type) === constructType
  )}

  <section class="space-y-4">
    <div class="flex flex-wrap items-center gap-3">
      <label class="flex-1 min-w-[16rem]">
        <span class="sr-only">Search</span>
        <input
          type="search"
          placeholder="Search... e.g. '227', 'Fd-3m', 'BNS 227.131', 'grey'"
          class="w-full rounded border border-slate-300 dark:border-slate-700 bg-transparent px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-sky-500"
          bind:value={query}
        />
      </label>
      <label class="flex items-center gap-2 text-sm">
        <span class="text-slate-500">Type:</span>
        <select
          class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1"
          bind:value={constructType}
        >
          <option value="All">All</option>
          <option value="1">I (Fedorov)</option>
          <option value="2">II (grey)</option>
          <option value="3">III (BW, equi-class)</option>
          <option value="4">IV (BW, equi-translation)</option>
        </select>
      </label>
      <span class="text-xs text-slate-500">{filtered.length} / {rows.length}</span>
    </div>

    <div class="overflow-x-auto rounded border border-slate-200 dark:border-slate-800">
      <table class="min-w-full text-sm">
        <thead class="bg-slate-50 dark:bg-slate-900 text-xs uppercase tracking-wide">
          <tr>
            <th class="px-3 py-2 text-left">UNI</th>
            <th class="px-3 py-2 text-left">Litvin</th>
            <th class="px-3 py-2 text-left">BNS</th>
            <th class="px-3 py-2 text-left">OG</th>
            <th class="px-3 py-2 text-left">Mag. Hall symbol</th>
            <th class="px-3 py-2 text-left">Parent SG</th>
            <th class="px-3 py-2 text-left">Type</th>
          </tr>
        </thead>
        <tbody>
          {#each filtered as r (r.uni_number)}
            <tr
              class="border-t border-slate-100 dark:border-slate-800 hover:bg-slate-50 dark:hover:bg-slate-900/60 cursor-pointer"
              tabindex="0"
              role="link"
              onclick={() => push(`/msg/${r.uni_number}`)}
              onkeydown={(e) => onRowKey(e, r.uni_number)}
            >
              <td class="px-3 py-1.5 font-mono">{r.uni_number}</td>
              <td class="px-3 py-1.5 font-mono">{r.litvin_number}</td>
              <td class="px-3 py-1.5 font-mono">{r.bns_number}</td>
              <td class="px-3 py-1.5 font-mono">{r.og_number}</td>
              <td class="px-3 py-1.5 font-mono">
                <a use:link href={`/msg/${r.uni_number}`} class="hover:underline"
                  >{r.magnetic_hall_symbol}</a
                >
              </td>
              <td class="px-3 py-1.5 font-mono">#{r.number} {r.parent_hm_short}</td>
              <td class="px-3 py-1.5">{r.construct_label}</td>
            </tr>
          {/each}
          {#if filtered.length === 0}
            <tr>
              <td colspan="7" class="px-3 py-6 text-center text-sm text-slate-500">
                No magnetic space groups match the current filter.
              </td>
            </tr>
          {/if}
        </tbody>
      </table>
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load magnetic-space-group catalog: ${formatErr(err)}`} />
{/await}
