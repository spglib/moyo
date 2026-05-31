<script lang="ts">
  import type {
    MoyoOperation,
    MoyoSpaceGroupType,
    MoyoHallSymbolEntry,
    MoyoArithmeticCrystalClass,
    MoyoWyckoffPosition,
  } from '@spglib/moyo-wasm'
  import { getMoyo, formatErr } from '../lib/wasm'
  import { SPACE_GROUP_COUNT, clampInt } from '../lib/catalog'
  import { settingDescription } from '../lib/format'
  import InfoGrid from '../components/InfoGrid.svelte'
  import OperationsTable from '../components/OperationsTable.svelte'
  import WyckoffTable from '../components/WyckoffTable.svelte'
  import CollapsibleSection from '../components/CollapsibleSection.svelte'
  import HmSymbol from '../components/HmSymbol.svelte'
  import ErrorCard from '../components/ErrorCard.svelte'
  import LoadingDots from '../components/LoadingDots.svelte'

  interface Setting {
    hall: MoyoHallSymbolEntry
    operations: MoyoOperation[]
    wyckoffs: MoyoWyckoffPosition[]
  }

  interface Loaded {
    type: MoyoSpaceGroupType
    arith: MoyoArithmeticCrystalClass
    defaultHallNumber: number
    settings: Setting[]
  }

  let { params }: { params: { number: string } } = $props()

  const number = $derived(clampInt(Number(params.number), 1, SPACE_GROUP_COUNT))
  const data = $derived(load(number))

  async function load(n: number): Promise<Loaded> {
    const m = await getMoyo()
    const type = m.space_group_type(n)
    const arith = m.arithmetic_crystal_class(type.arithmetic_number)
    const settings: Setting[] = m.hall_symbol_entries_from_number(n).map((hall) => ({
      hall,
      operations: m.operations_from_number(
        n,
        { type: 'HallNumber', value: hall.hall_number },
        false
      ),
      wyckoffs: m.wyckoff_positions(hall.hall_number),
    }))
    return { type, arith, defaultHallNumber: type.hall_number, settings }
  }

  function settingTitle(setting: string): string {
    const desc = settingDescription(setting)
    return setting ? `${desc} (${setting})` : desc
  }
</script>

{#await data}
  <LoadingDots />
{:then d}
  <header class="py-2 border-b border-stone-300 dark:border-stone-700 mb-4">
    <div class="eyebrow">Space group</div>
    <h1 class="text-2xl font-semibold flex items-baseline gap-2">
      <span class="font-mono">#{d.type.number}</span>
      <HmSymbol symbol={d.type.hm_short} />
    </h1>
    <p class="text-sm text-stone-600 dark:text-stone-400 font-mono">{d.type.hm_full}</p>
  </header>

  <section class="space-y-6">
    <InfoGrid
      rows={[
        { label: 'ITA number', value: d.type.number, mono: true },
        { label: 'Short Hermann-Mauguin symbol', value: d.type.hm_short, mono: true },
        { label: 'Full Hermann-Mauguin symbol', value: d.type.hm_full, mono: true },
        { label: 'Crystal family', value: d.type.crystal_family },
        { label: 'Crystal system', value: d.type.crystal_system },
        { label: 'Lattice system', value: d.type.lattice_system },
        { label: 'Bravais class', value: d.type.bravais_class, mono: true },
        {
          label: 'Arithmetic crystal class',
          value: `${d.arith.arithmetic_number} (${d.arith.symbol})`,
          mono: true,
        },
        { label: 'Geometric crystal class', value: d.type.geometric_crystal_class, mono: true },
      ]}
    />

    <div>
      <h2 class="text-lg font-semibold mb-2">
        Settings <span class="text-sm font-normal text-stone-500">({d.settings.length})</span>
      </h2>
      <div class="space-y-3">
        {#each d.settings as s (s.hall.hall_number)}
          {@const isStandard = s.hall.hall_number === d.defaultHallNumber}
          <div
            class="rounded border p-3 {isStandard
              ? 'border-emerald-400 dark:border-emerald-600'
              : 'border-stone-200 dark:border-stone-800'}"
          >
            <CollapsibleSection
              title={settingTitle(s.hall.setting)}
              badge={isStandard ? 'ITA standard' : undefined}
              open={isStandard}
            >
              <div class="space-y-4">
                <InfoGrid
                  rows={[
                    { label: 'Hall number', value: s.hall.hall_number, mono: true },
                    { label: 'Hall symbol', value: s.hall.hall_symbol, mono: true },
                    { label: 'Centering', value: s.hall.centering, mono: true },
                  ]}
                />
                <CollapsibleSection title="Symmetry operations" count={s.operations.length}>
                  <OperationsTable operations={s.operations} />
                </CollapsibleSection>
                <CollapsibleSection title="Wyckoff positions" count={s.wyckoffs.length}>
                  <WyckoffTable positions={s.wyckoffs} />
                </CollapsibleSection>
              </div>
            </CollapsibleSection>
          </div>
        {/each}
      </div>
    </div>
  </section>
{:catch err}
  <ErrorCard message={`Failed to load space group ${number}: ${formatErr(err)}`} />
{/await}
