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
  <p class="mt-1 text-slate-600 dark:text-slate-400 max-w-2xl">
    Static viewer powered by <code class="font-mono text-sm">moyo</code> compiled to WebAssembly. All
    look-ups run in your browser; no server, no analytics.
  </p>

  <div class="mt-6 grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
    {#each cards as c}
      <a
        use:link
        href={c.href}
        class="block rounded-lg border border-slate-200 dark:border-slate-800 p-5 hover:border-slate-400 dark:hover:border-slate-600 transition-colors"
      >
        <div class="flex items-baseline justify-between">
          <h2 class="text-lg font-medium">{c.title}</h2>
          <span class="text-xs text-slate-500 dark:text-slate-400 font-mono">{c.range}</span>
        </div>
        <p class="mt-2 text-sm text-slate-600 dark:text-slate-400">{c.blurb}</p>
      </a>
    {/each}
  </div>
</section>
