<script lang="ts">
  import { push, replace } from 'svelte-spa-router'
  import type {
    MoyoSetting,
    MoyoOperation,
    MoyoSpaceGroupType,
    MoyoHallSymbolEntry,
    MoyoArithmeticCrystalClass,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { getHallNumbersBySpaceGroup } from '../lib/hall'
  import { SPACE_GROUP_COUNT, clampInt, systemOfSpaceGroup } from '../lib/catalog'
  import { parseQuery, buildQuery } from '../lib/url'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import GroupPager from '../components/GroupPager.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  type SettingKind = 'Standard' | 'Spglib' | 'HallNumber'
  interface Loaded {
    type: MoyoSpaceGroupType
    hall: MoyoHallSymbolEntry | null
    arith: MoyoArithmeticCrystalClass
    operations: MoyoOperation[]
    hallList: number[]
  }

  let { params, querystring }: { params: { number: string }; querystring?: string } = $props()

  const number = $derived(clampInt(Number(params.number), 1, SPACE_GROUP_COUNT))
  const query = $derived(parseQuery(querystring))
  const settingKind = $derived<SettingKind>(((): SettingKind => {
    const s = query.get('setting')
    if (s === 'Spglib' || s === 'HallNumber') return s
    return 'Standard'
  })())
  const hallParam = $derived(Number(query.get('hall')))
  const primitive = $derived(query.get('primitive') === '1')

  const setting = $derived<MoyoSetting>(
    settingKind === 'HallNumber' && Number.isFinite(hallParam)
      ? { type: 'HallNumber', value: hallParam }
      : settingKind === 'Spglib'
        ? { type: 'Spglib' }
        : { type: 'Standard' }
  )

  const data = $derived(load(number, setting, primitive))

  async function load(n: number, s: MoyoSetting, prim: boolean): Promise<Loaded> {
    const m = await getMoyo()
    const hallList = (await getHallNumbersBySpaceGroup()).get(n) ?? []
    const type = m.space_group_type(n)
    const effective: MoyoSetting =
      s.type === 'HallNumber' && !hallList.includes(s.value) ? { type: 'Standard' } : s
    const operations = m.operations_from_number(n, effective, prim)
    const activeHall =
      effective.type === 'HallNumber'
        ? effective.value
        : (hallList[0] ?? 0)
    const hall = activeHall ? m.hall_symbol_entry(activeHall) : null
    const arith = m.arithmetic_crystal_class(type.arithmetic_number)
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
    replace(`/sg/${number}${q ? `?${q}` : ''}`)
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
      <h1 class="text-2xl font-semibold">
        <span class="font-mono">#{d.type.number}</span>
        <span class="ml-2">{d.type.hm_short}</span>
      </h1>
      <p class="text-sm text-slate-600 dark:text-slate-400 font-mono">
        {d.type.hm_full}
        {#if d.hall}&middot; Hall: {d.hall.hall_symbol}{/if}
      </p>
    </div>
    <GroupPager value={number} min={1} max={SPACE_GROUP_COUNT} basePath="/sg" label="#" />
  </header>

  <section class="grid grid-cols-1 lg:grid-cols-3 gap-6">
    <div class="lg:col-span-2 space-y-6">
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
          { label: 'Centering', value: d.hall ? (centeringLabels[d.hall.centering] ?? d.hall.centering) : '-' },
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

        {#if systemOfSpaceGroup(number)}
          <span class="ml-auto text-xs text-slate-500">{systemOfSpaceGroup(number)} system</span>
        {/if}
      </div>

      <OperationsTable operations={d.operations} />
    </div>

    <aside class="space-y-3 lg:sticky lg:top-20 self-start">
      <h3 class="text-xs uppercase tracking-wide text-slate-500">Crystal-system boundaries</h3>
      <div class="flex flex-wrap gap-1">
        {#each [1, 2, 3, 15, 16, 74, 75, 142, 143, 167, 168, 194, 195, 230] as n}
          <button
            type="button"
            class="rounded border border-slate-300 dark:border-slate-700 px-2 py-1 font-mono text-xs hover:bg-slate-100 dark:hover:bg-slate-800"
            onclick={() => push(`/sg/${n}`)}
          >
            {n}
          </button>
        {/each}
      </div>
    </aside>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load space group ${number}: ${formatErr(err)}`} />
{/await}
