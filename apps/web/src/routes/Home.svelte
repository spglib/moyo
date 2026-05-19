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
