<script lang="ts">
  import Router from 'svelte-spa-router'
  import { routes } from './router'
  import Header from './components/Header.svelte'

  let { } = $props()
</script>

<div class="min-h-screen flex flex-col">
  <Header />
  <main class="flex-1 mx-auto w-full max-w-6xl px-4 py-6">
    <Router {routes} />
  </main>
  <footer class="mx-auto w-full max-w-6xl px-4 py-6 text-xs text-slate-500 dark:text-slate-400">
    Powered by <a class="underline" href="https://github.com/spglib/moyo">moyo</a> via WebAssembly.
  </footer>
</div>
