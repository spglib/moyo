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

<header class="border-b border-slate-200 dark:border-slate-800 bg-white/80 dark:bg-slate-950/80 backdrop-blur sticky top-0 z-10">
  <div class="mx-auto max-w-6xl px-4 py-3 flex items-center gap-4">
    <a use:link href="/" class="font-semibold tracking-tight">moyo · Crystallographic Group Viewer</a>
    <nav class="ml-2 flex gap-1 text-sm">
      {#each tabs as tab}
        {@const active = tab.match.test($location)}
        <a
          use:link
          href={tab.href}
          class="px-2 py-1 rounded transition-colors {active
            ? 'bg-slate-900 text-white dark:bg-slate-100 dark:text-slate-900'
            : 'text-slate-600 hover:bg-slate-100 dark:text-slate-300 dark:hover:bg-slate-800'}"
        >
          {tab.label}
        </a>
      {/each}
    </nav>
    <div class="ml-auto"><ThemeToggle /></div>
  </div>
</header>
