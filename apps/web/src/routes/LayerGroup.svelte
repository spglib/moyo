<script lang="ts">
  import type {
    MoyoOperation,
    MoyoLayerGroupType,
    MoyoLayerHallSymbolEntry,
    MoyoLayerArithmeticCrystalClass,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { getLayerHallNumbersByGroup } from '../lib/hall'
  import { LAYER_GROUP_COUNT, clampInt } from '../lib/catalog'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import GroupPager from '../components/GroupPager.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  interface Loaded {
    type: MoyoLayerGroupType
    hall: MoyoLayerHallSymbolEntry | null
    arith: MoyoLayerArithmeticCrystalClass
    operations: MoyoOperation[]
  }

  let { params }: { params: { number: string } } = $props()

  const number = $derived(clampInt(Number(params.number), 1, LAYER_GROUP_COUNT))
  const data = $derived(load(number))

  async function load(n: number): Promise<Loaded> {
    const m = await getMoyo()
    const hallList = (await getLayerHallNumbersByGroup()).get(n) ?? []
    const type = m.layer_group_type(n)
    const operations = m.operations_from_layer_number(n, { type: 'Standard' }, false)
    const hall = hallList.length > 0 ? m.layer_hall_symbol_entry(hallList[0]) : null
    const arith = m.layer_arithmetic_crystal_class(type.arithmetic_number)
    return { type, hall, arith, operations }
  }
</script>

{#await data}
  <LoadingDots />
{:then d}
  <header
    class="flex flex-wrap items-end justify-between gap-4 py-2 border-b border-slate-200 dark:border-slate-800 mb-4"
  >
    <div>
      <div class="text-xs uppercase tracking-wide text-slate-500">Layer group</div>
      <h1 class="text-2xl font-semibold">
        <span class="font-mono">#{d.type.number}</span>
        <span class="ml-2">{d.type.hm_short}</span>
      </h1>
      <p class="text-sm text-slate-600 dark:text-slate-400 font-mono">
        {d.type.hm_full}
        {#if d.hall}&middot; Hall: {d.hall.hall_symbol}{/if}
      </p>
    </div>
    <GroupPager value={number} min={1} max={LAYER_GROUP_COUNT} basePath="/lg" label="#" />
  </header>

  <section class="grid grid-cols-1 lg:grid-cols-3 gap-6">
    <div class="lg:col-span-2 space-y-6">
      <InfoGrid
        rows={[
          { label: 'Layer number', value: d.type.number, mono: true },
          { label: 'Hall number', value: d.hall?.hall_number ?? '-', mono: true },
          { label: 'HM short', value: d.type.hm_short, mono: true },
          { label: 'HM full', value: d.type.hm_full, mono: true },
          { label: 'Hall symbol', value: d.hall?.hall_symbol ?? '-', mono: true },
          { label: 'Setting (Hall row)', value: d.hall?.setting ?? '-', mono: true },
          { label: 'Lattice system', value: d.type.lattice_system },
          { label: 'Bravais class', value: d.type.bravais_class, mono: true },
          { label: 'Centering', value: d.hall?.centering ?? '-', mono: true },
          {
            label: 'Arithmetic crystal class',
            value: `${d.arith.arithmetic_number} (${d.arith.symbol})`,
            mono: true,
          },
          { label: 'Geometric crystal class', value: d.type.geometric_crystal_class, mono: true },
        ]}
      />

      <OperationsTable operations={d.operations} />
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load layer group ${number}: ${formatErr(err)}`} />
{/await}
