<script lang="ts">
  import { link } from 'svelte-spa-router'

  const cards = [
    {
      href: '/sg',
      title: 'Space groups',
      range: '1 - 230',
      blurb: 'The 230 crystallographic space-group types in 3-D.',
    },
    {
      href: '/lg',
      title: 'Layer groups',
      range: '1 - 80',
      blurb: 'The 80 layer-group types (subperiodic groups in 3-D with 2-D translations).',
    },
    {
      href: '/msg',
      title: 'Magnetic space groups',
      range: '1 - 1651',
      blurb: 'The 1651 Shubnikov / magnetic space-group types, indexed by UNI number.',
    },
  ]

  const sources: { range: string; title: string; refs: { text: string; href?: string }[] }[] = [
    {
      range: '1 - 230',
      title: 'Space groups',
      refs: [
        { text: 'Hall-symbol database curated by Dr. Yusuke Seto', href: 'https://yseto.net/en/sg/sg1' },
        {
          text: 'Hermann-Mauguin symbols and settings: International Tables for Crystallography, Vol. A (2016)',
        },
      ],
    },
    {
      range: '1 - 80',
      title: 'Layer groups',
      refs: [
        { text: 'Group tabulation: International Tables for Crystallography, Vol. E (2010)' },
        {
          text: 'Hall-symbol assignments: Fu et al., "Symmetry classification of 2D materials: layer groups versus space groups", 2D Mater. 11, 035009 (2024)',
          href: 'https://doi.org/10.1088/2053-1583/ad3e0c',
        },
      ],
    },
    {
      range: '1 - 1651',
      title: 'Magnetic space groups',
      refs: [
        {
          text: 'Magnetic Hall symbols: Gonzalez-Platas, Katcho & Rodriguez-Carvajal, "Extension of Hall symbols of crystallographic space groups to magnetic space groups", J. Appl. Cryst. (2021)',
          href: 'https://journals.iucr.org/paper?tu5004',
        },
        {
          text: 'UNI/BNS/OG numbering and Litvin numbers: Litvin, Magnetic Group Tables (IUCr, 2013)',
          href: 'https://www.iucr.org/publ/978-0-9553602-2-0',
        },
      ],
    },
  ]
</script>

<section class="py-4">
  <h1 class="text-2xl font-semibold tracking-tight">Browse crystallographic groups</h1>
  <p class="mt-1 text-stone-600 dark:text-stone-400 max-w-2xl">
    Static viewer powered by <code class="font-mono text-sm">moyo</code> compiled to WebAssembly.
    All look-ups run in your browser; no server, no analytics.
  </p>

  <div
    class="mt-6 border-y border-stone-300 dark:border-stone-700 divide-y divide-stone-200 dark:divide-stone-800"
  >
    {#each cards as c}
      <a
        use:link
        href={c.href}
        class="grid grid-cols-[7rem_1fr_auto] items-baseline gap-4 px-1 py-4 hover:bg-moyo-50 dark:hover:bg-moyo-950/30"
      >
        <span class="font-mono text-xs text-moyo-700 dark:text-moyo-400">{c.range}</span>
        <span>
          <span class="block font-medium">{c.title}</span>
          <span class="block text-sm text-stone-600 dark:text-stone-400">{c.blurb}</span>
        </span>
        <span class="font-mono text-xs text-stone-400">&rarr;</span>
      </a>
    {/each}
  </div>
</section>

<section class="py-4">
  <h2 class="text-lg font-semibold tracking-tight">Data sources</h2>
  <p class="mt-1 text-sm text-stone-600 dark:text-stone-400 max-w-2xl">
    Crystallographic data embedded in this viewer is sourced from the following references.
  </p>

  <div
    class="mt-4 border-y border-stone-300 dark:border-stone-700 divide-y divide-stone-200 dark:divide-stone-800"
  >
    {#each sources as s}
      <div class="grid grid-cols-[7rem_1fr] items-baseline gap-4 px-1 py-4">
        <span class="font-mono text-xs text-moyo-700 dark:text-moyo-400">{s.range}</span>
        <div>
          <span class="block font-medium">{s.title}</span>
          <ul class="mt-1 space-y-1 text-sm text-stone-600 dark:text-stone-400 list-disc pl-5">
            {#each s.refs as r}
              <li>
                {#if r.href}
                  <a
                    href={r.href}
                    target="_blank"
                    rel="noopener noreferrer"
                    class="underline hover:text-moyo-700 dark:hover:text-moyo-400"
                  >
                    {r.text}
                  </a>
                {:else}
                  {r.text}
                {/if}
              </li>
            {/each}
          </ul>
        </div>
      </div>
    {/each}
  </div>

  <p class="mt-4 text-xs text-stone-500 dark:text-stone-400">
    <span class="font-medium">Source code:</span>
    <a
      href="https://github.com/spglib/moyo"
      target="_blank"
      rel="noopener noreferrer"
      class="underline hover:text-moyo-700 dark:hover:text-moyo-400">github.com/spglib/moyo</a
    >
    &middot; <span class="font-medium">License:</span>
    <a
      href="https://github.com/spglib/moyo/blob/main/LICENSE-MIT"
      target="_blank"
      rel="noopener noreferrer"
      class="underline hover:text-moyo-700 dark:hover:text-moyo-400">MIT</a
    >
    or
    <a
      href="https://github.com/spglib/moyo/blob/main/LICENSE-APACHE"
      target="_blank"
      rel="noopener noreferrer"
      class="underline hover:text-moyo-700 dark:hover:text-moyo-400">Apache-2.0</a
    >
  </p>
</section>
