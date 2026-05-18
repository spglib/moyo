<script lang="ts">
  import { replace } from 'svelte-spa-router'
  import type {
    MoyoMagneticHallSymbolEntry,
    MoyoMagneticOperation,
    MoyoMagneticSpaceGroupType,
    MoyoSpaceGroupType,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { MAGNETIC_SG_COUNT, clampInt } from '../lib/catalog'
  import { parseQuery, buildQuery } from '../lib/url'
  import { constructTypeLabel } from '../lib/format'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import GroupPager from '../components/GroupPager.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  interface Loaded {
    type: MoyoMagneticSpaceGroupType
    hall: MoyoMagneticHallSymbolEntry
    parent: MoyoSpaceGroupType
    operations: MoyoMagneticOperation[]
  }

  let { params, querystring }: { params: { uni_number: string }; querystring?: string } = $props()

  const uni = $derived(clampInt(Number(params.uni_number), 1, MAGNETIC_SG_COUNT))
  const query = $derived(parseQuery(querystring))
  const primitive = $derived(query.get('primitive') === '1')

  const data = $derived(load(uni, primitive))

  async function load(n: number, prim: boolean): Promise<Loaded> {
    const m = await getMoyo()
    const type = m.magnetic_space_group_type(n)
    const hall = m.magnetic_hall_symbol_entry(n)
    const parent = m.space_group_type(type.number)
    const operations = m.magnetic_operations_from_uni_number(n, prim)
    return { type, hall, parent, operations }
  }

  function pushQuery(next: { primitive?: boolean }) {
    const merged = { primitive: next.primitive ?? primitive }
    const q = buildQuery({ primitive: merged.primitive ? '1' : '' })
    replace(`/msg/${uni}${q ? `?${q}` : ''}`)
  }
</script>

{#await data}
  <LoadingDots />
{:then d}
  <header
    class="flex flex-wrap items-end justify-between gap-4 py-2 border-b border-slate-200 dark:border-slate-800 mb-4"
  >
    <div>
      <div class="text-xs uppercase tracking-wide text-slate-500">Magnetic space group</div>
      <h1 class="text-2xl font-semibold">
        <span class="font-mono">UNI #{d.type.uni_number}</span>
        <span class="ml-2 font-mono">{d.hall.magnetic_hall_symbol}</span>
      </h1>
      <p class="text-sm text-slate-600 dark:text-slate-400 font-mono">
        BNS {d.type.bns_number} &middot; OG {d.type.og_number} &middot; parent SG #{d.type.number}
        {d.parent.hm_short}
      </p>
    </div>
    <GroupPager value={uni} min={1} max={MAGNETIC_SG_COUNT} basePath="/msg" label="UNI" />
  </header>

  <section class="space-y-6">
    <InfoGrid
      rows={[
        { label: 'UNI number', value: d.type.uni_number, mono: true },
        { label: 'Litvin number', value: d.type.litvin_number, mono: true },
        { label: 'BNS number', value: d.type.bns_number, mono: true },
        { label: 'OG number', value: d.type.og_number, mono: true },
        { label: 'Magnetic Hall symbol', value: d.hall.magnetic_hall_symbol, mono: true },
        { label: 'Construct type', value: constructTypeLabel(d.type.construct_type) },
        {
          label: 'Parent space group',
          value: `#${d.parent.number} ${d.parent.hm_short}`,
          mono: true,
        },
        { label: 'Parent crystal system', value: d.parent.crystal_system },
      ]}
    />

    <div
      class="rounded border border-slate-200 dark:border-slate-800 p-4 bg-slate-50 dark:bg-slate-900/40 flex flex-wrap items-center gap-4 text-sm"
    >
      <label class="flex items-center gap-2">
        <input
          type="checkbox"
          checked={primitive}
          onchange={(e) =>
            pushQuery({ primitive: (e.currentTarget as HTMLInputElement).checked })}
        />
        <span>Primitive cell</span>
      </label>
    </div>

    <OperationsTable operations={d.operations} hasTimeReversal />
  </section>
{:catch err}
  <ErrorCard message={`Failed to load magnetic space group ${uni}: ${formatErr(err)}`} />
{/await}
