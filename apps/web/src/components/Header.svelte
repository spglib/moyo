<script lang="ts">
  import { link, location } from 'svelte-spa-router'
  import ThemeToggle from './ThemeToggle.svelte'

  const tabs = [
    { href: '/', label: 'Home', match: /^\/?$/ },
    { href: '/sg', label: 'Space groups', match: /^\/sg(\/|$)/ },
    { href: '/lg', label: 'Layer groups', match: /^\/lg(\/|$)/ },
    { href: '/msg', label: 'Magnetic SGs', match: /^\/msg(\/|$)/ },
  ]
</script>

<header
  class="border-b border-stone-300 dark:border-stone-700 bg-white/80 dark:bg-stone-950/80 backdrop-blur sticky top-0 z-10"
>
  <div class="mx-auto max-w-6xl px-4 py-3 flex items-center gap-6">
    <a use:link href="/" class="font-semibold tracking-tight">
      <span class="text-moyo-700 dark:text-moyo-400">moyo</span> · Crystallographic Group Viewer
    </a>
    <nav class="ml-2 flex gap-4 text-sm">
      {#each tabs as tab}
        {@const active = tab.match.test($location)}
        <a
          use:link
          href={tab.href}
          class="border-b-2 -mb-px py-1 transition-colors {active
            ? 'border-moyo-600 text-moyo-700 dark:text-moyo-400'
            : 'border-transparent text-stone-600 hover:text-stone-950 dark:text-stone-300 dark:hover:text-white'}"
        >
          {tab.label}
        </a>
      {/each}
    </nav>
    <div class="ml-auto"><ThemeToggle /></div>
  </div>
</header>
