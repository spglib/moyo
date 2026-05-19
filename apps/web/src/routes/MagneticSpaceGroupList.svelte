<script lang="ts">
  import { push, link } from 'svelte-spa-router'
  import { getAllMagneticSpaceGroups, filterRows } from '../lib/catalogIndex'
  import { formatErr } from '../lib/wasm'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  const PAGE_SIZE = 200

  let query = $state('')
  let reference = $state<string>('All')
  let limit = $state(PAGE_SIZE)

  const data = getAllMagneticSpaceGroups()

  $effect(() => {
    void query
    void reference
    limit = PAGE_SIZE
  })

  function onRowKey(e: KeyboardEvent, n: number) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault()
      push(`/msg/${n}`)
    }
  }
</script>

<header class="page-header">
  <div>
    <div class="eyebrow">Magnetic space groups</div>
    <h1 class="text-2xl font-semibold">All 1651 magnetic space-group types</h1>
    <p class="text-sm text-stone-600 dark:text-stone-400">
      Search by UNI / Litvin / BNS / OG number, magnetic Hall symbol, ...
    </p>
  </div>
</header>

{#await data}
  <LoadingDots />
{:then rows}
  {@const refOptions = Array.from(
    new Map(rows.map((r) => [r.number, r.parent_hm_short])).entries()
  ).sort((a, b) => a[0] - b[0])}
  {@const filtered = filterRows(rows, query).filter(
    (r) => reference === 'All' || String(r.number) === reference
  )}
  {@const visible = filtered.slice(0, limit)}

  <section class="space-y-4">
    <div class="flex flex-wrap items-center gap-3">
      <label class="flex-1 min-w-[16rem]">
        <span class="sr-only">Search</span>
        <input
          type="search"
          placeholder="Search... e.g. '227', 'Fd-3m', 'BNS 227.131', 'grey'"
          class="search-input"
          bind:value={query}
        />
      </label>
      <span class="text-xs text-stone-500"
        >{visible.length === filtered.length
          ? `${filtered.length} / ${rows.length}`
          : `${visible.length} of ${filtered.length} (of ${rows.length})`}</span
      >
    </div>

    <div class="flex flex-wrap items-center gap-x-4 gap-y-2 text-sm">
      <label class="flex items-center gap-2">
        <span class="text-stone-500">Reference space group in BNS setting:</span>
        <select class="filter-select font-mono" bind:value={reference}>
          <option value="All">All</option>
          {#each refOptions as [n, hm]}
            <option value={String(n)}>#{n} {hm}</option>
          {/each}
        </select>
      </label>
      {#if reference !== 'All'}
        <button
          type="button"
          class="link-button"
          onclick={() => (reference = 'All')}
        >
          Reset filters
        </button>
      {/if}
    </div>

    <div class="table-shell">
      <table class="min-w-full text-sm">
        <thead class="table-head">
          <tr>
            <th class="px-3 py-2 text-left">UNI number</th>
            <th class="px-3 py-2 text-left">Litvin number</th>
            <th class="px-3 py-2 text-left">BNS number</th>
            <th class="px-3 py-2 text-left">OG number</th>
            <th class="px-3 py-2 text-left">Magnetic Hall symbol</th>
            <th class="px-3 py-2 text-left">Reference space group in BNS setting</th>
            <th class="px-3 py-2 text-left">Construct type</th>
          </tr>
        </thead>
        <tbody>
          {#each visible as r (r.uni_number)}
            <tr
              class="table-row-link"
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
              <td colspan="7" class="px-3 py-6 text-center text-sm text-stone-500">
                No magnetic space groups match the current filter.
              </td>
            </tr>
          {/if}
        </tbody>
      </table>
    </div>

    {#if visible.length < filtered.length}
      <div class="text-center">
        <button
          type="button"
          class="rounded border border-stone-300 dark:border-stone-700 px-3 py-1 text-sm hover:bg-stone-100 dark:hover:bg-stone-800"
          onclick={() => (limit += PAGE_SIZE)}
        >
          Show {Math.min(PAGE_SIZE, filtered.length - visible.length)} more
        </button>
      </div>
    {/if}
  </section>
{:catch err}
  <ErrorCard message={`Failed to load magnetic-space-group catalog: ${formatErr(err)}`} />
{/await}
