<script lang="ts">
  import { replace } from 'svelte-spa-router'
  import type {
    MoyoLayerSetting,
    MoyoOperation,
    MoyoLayerGroupType,
    MoyoLayerHallSymbolEntry,
    MoyoLayerArithmeticCrystalClass,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { getLayerHallNumbersByGroup } from '../lib/hall'
  import { LAYER_GROUP_COUNT, clampInt } from '../lib/catalog'
  import { parseQuery, buildQuery } from '../lib/url'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import GroupPager from '../components/GroupPager.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  type SettingKind = 'Standard' | 'Spglib' | 'HallNumber'
  interface Loaded {
    type: MoyoLayerGroupType
    hall: MoyoLayerHallSymbolEntry | null
    arith: MoyoLayerArithmeticCrystalClass
    operations: MoyoOperation[]
    hallList: number[]
  }

  let { params, querystring }: { params: { number: string }; querystring?: string } = $props()

  const number = $derived(clampInt(Number(params.number), 1, LAYER_GROUP_COUNT))
  const query = $derived(parseQuery(querystring))
  const settingKind = $derived<SettingKind>(((): SettingKind => {
    const s = query.get('setting')
    if (s === 'Spglib' || s === 'HallNumber') return s
    return 'Standard'
  })())
  const hallParam = $derived(Number(query.get('hall')))
  const primitive = $derived(query.get('primitive') === '1')

  const setting = $derived<MoyoLayerSetting>(
    settingKind === 'HallNumber' && Number.isFinite(hallParam)
      ? { type: 'HallNumber', value: hallParam }
      : settingKind === 'Spglib'
        ? { type: 'Spglib' }
        : { type: 'Standard' }
  )

  const data = $derived(load(number, setting, primitive))

  async function load(n: number, s: MoyoLayerSetting, prim: boolean): Promise<Loaded> {
    const m = await getMoyo()
    const hallList = (await getLayerHallNumbersByGroup()).get(n) ?? []
    const type = m.layer_group_type(n)
    const effective: MoyoLayerSetting =
      s.type === 'HallNumber' && !hallList.includes(s.value) ? { type: 'Standard' } : s
    const operations = m.operations_from_layer_number(n, effective, prim)
    const activeHall =
      effective.type === 'HallNumber'
        ? effective.value
        : (hallList[0] ?? 0)
    const hall = activeHall ? m.layer_hall_symbol_entry(activeHall) : null
    const arith = m.layer_arithmetic_crystal_class(type.arithmetic_number)
    return { type, hall, arith, operations, hallList }
  }

  function pushQuery(next: Partial<{ setting: SettingKind; hall: number | null; primitive: boolean }>) {
    const merged = {
      setting: next.setting ?? settingKind,
      hall:
        next.hall === null
          ? null
          : (next.hall ?? (settingKind === 'HallNumber' ? hallParam : null)),
      primitive: next.primitive ?? primitive,
    }
    const q = buildQuery({
      setting: merged.setting === 'Standard' ? '' : merged.setting,
      hall: merged.setting === 'HallNumber' ? (merged.hall ?? '') : '',
      primitive: merged.primitive ? '1' : '',
    })
    replace(`/lg/${number}${q ? `?${q}` : ''}`)
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

      <div
        class="rounded border border-slate-200 dark:border-slate-800 p-4 bg-slate-50 dark:bg-slate-900/40 flex flex-wrap items-center gap-4 text-sm"
      >
        <label class="flex items-center gap-2">
          <span class="text-slate-500">Setting:</span>
          <select
            class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1"
            value={settingKind}
            onchange={(e) =>
              pushQuery({
                setting: (e.currentTarget as HTMLSelectElement).value as SettingKind,
                hall: null,
              })}
          >
            <option value="Standard">Standard</option>
            <option value="Spglib">Spglib</option>
            <option value="HallNumber" disabled={d.hallList.length === 0}>HallNumber</option>
          </select>
        </label>

        {#if settingKind === 'HallNumber'}
          <label class="flex items-center gap-2">
            <span class="text-slate-500">Hall #:</span>
            <select
              class="rounded border border-slate-300 dark:border-slate-700 bg-transparent px-2 py-1 font-mono"
              value={d.hall?.hall_number ?? d.hallList[0]}
              onchange={(e) =>
                pushQuery({ hall: Number((e.currentTarget as HTMLSelectElement).value) })}
            >
              {#each d.hallList as h}
                <option value={h}>{h}</option>
              {/each}
            </select>
          </label>
        {/if}

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

      <OperationsTable operations={d.operations} />
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load layer group ${number}: ${formatErr(err)}`} />
{/await}
