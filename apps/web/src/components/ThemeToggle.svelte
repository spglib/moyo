<script lang="ts">
  let dark = $state(initial())

  function initial(): boolean {
    if (typeof window === 'undefined') return false
    const stored = localStorage.getItem('moyo-theme')
    if (stored === 'dark') return true
    if (stored === 'light') return false
    return window.matchMedia('(prefers-color-scheme: dark)').matches
  }

  $effect(() => {
    const root = document.documentElement
    root.classList.toggle('dark', dark)
    localStorage.setItem('moyo-theme', dark ? 'dark' : 'light')
  })
</script>

<button
  type="button"
  class="rounded border border-stone-300 dark:border-stone-700 px-2 py-1 text-xs hover:bg-stone-100 dark:hover:bg-stone-800"
  onclick={() => (dark = !dark)}
  aria-label="Toggle dark mode"
>
  {dark ? 'Light' : 'Dark'}
</button>
