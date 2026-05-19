<script lang="ts">
  import type {
    MoyoOperation,
    MoyoSpaceGroupType,
    MoyoHallSymbolEntry,
    MoyoArithmeticCrystalClass,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { getHallNumbersBySpaceGroup } from '../lib/hall'
  import { SPACE_GROUP_COUNT, clampInt } from '../lib/catalog'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import GroupPager from '../components/GroupPager.svelte'
  import HmSymbol from '../components/HmSymbol.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  interface Loaded {
    type: MoyoSpaceGroupType
    hall: MoyoHallSymbolEntry | null
    arith: MoyoArithmeticCrystalClass
    operations: MoyoOperation[]
  }

  let { params }: { params: { number: string } } = $props()

  const number = $derived(clampInt(Number(params.number), 1, SPACE_GROUP_COUNT))
  const data = $derived(load(number))

  async function load(n: number): Promise<Loaded> {
    const m = await getMoyo()
    const hallList = (await getHallNumbersBySpaceGroup()).get(n) ?? []
    const type = m.space_group_type(n)
    const operations = m.operations_from_number(n, { type: 'Standard' }, false)
    const hall = hallList.length > 0 ? m.hall_symbol_entry(hallList[0]) : null
    const arith = m.arithmetic_crystal_class(type.arithmetic_number)
    return { type, hall, arith, operations }
  }

  const centeringLabels: Record<string, string> = {
    P: 'P (primitive)',
    A: 'A',
    B: 'B',
    C: 'C',
    I: 'I (body-centred)',
    R: 'R (rhombohedral)',
    F: 'F (face-centred)',
  }
</script>

{#await data}
  <LoadingDots />
{:then d}
  <header
    class="flex flex-wrap items-end justify-between gap-4 py-2 border-b border-slate-200 dark:border-slate-800 mb-4"
  >
    <div>
      <div class="text-xs uppercase tracking-wide text-slate-500">Space group</div>
      <h1 class="text-2xl font-semibold flex items-baseline gap-2">
        <span class="font-mono">#{d.type.number}</span>
        <HmSymbol symbol={d.type.hm_short} />
      </h1>
      <p class="text-sm text-slate-600 dark:text-slate-400 font-mono">
        {d.type.hm_full}
        {#if d.hall}&middot; Hall: {d.hall.hall_symbol}{/if}
      </p>
    </div>
    <GroupPager value={number} min={1} max={SPACE_GROUP_COUNT} basePath="/sg" label="#" />
  </header>

  <section class="space-y-6">
    <InfoGrid
      rows={[
        { label: 'ITA number', value: d.type.number, mono: true },
        { label: 'Hall number', value: d.hall?.hall_number ?? '-', mono: true },
        { label: 'HM short', value: d.type.hm_short, mono: true },
        { label: 'HM full', value: d.type.hm_full, mono: true },
        { label: 'Hall symbol', value: d.hall?.hall_symbol ?? '-', mono: true },
        { label: 'Setting (Hall row)', value: d.hall?.setting ?? '-', mono: true },
        { label: 'Crystal family', value: d.type.crystal_family },
        { label: 'Crystal system', value: d.type.crystal_system },
        { label: 'Lattice system', value: d.type.lattice_system },
        { label: 'Bravais class', value: d.type.bravais_class, mono: true },
        {
          label: 'Centering',
          value: d.hall ? (centeringLabels[d.hall.centering] ?? d.hall.centering) : '-',
        },
        {
          label: 'Arithmetic crystal class',
          value: `${d.arith.arithmetic_number} (${d.arith.symbol})`,
          mono: true,
        },
        { label: 'Geometric crystal class', value: d.type.geometric_crystal_class, mono: true },
      ]}
    />

    <OperationsTable operations={d.operations} />
  </section>
{:catch err}
  <ErrorCard message={`Failed to load space group ${number}: ${formatErr(err)}`} />
{/await}
