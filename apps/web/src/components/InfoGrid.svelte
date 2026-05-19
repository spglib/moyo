<script lang="ts">
  interface Row {
    label: string
    value: string | number
    mono?: boolean
  }
  let { rows }: { rows: Row[] } = $props()
</script>

<dl class="grid grid-cols-1 sm:grid-cols-2 gap-x-6 gap-y-2 text-sm">
  {#each rows as r}
    <div class="flex flex-col border-b border-stone-200 dark:border-stone-800 py-1">
      <dt class="text-xs uppercase tracking-wide text-stone-500 dark:text-stone-400">{r.label}</dt>
      <dd class={r.mono ? 'font-mono text-base' : 'text-base'}>{r.value}</dd>
    </div>
  {/each}
</dl>
