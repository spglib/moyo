<script lang="ts">
  import type { Snippet } from 'svelte'

  let {
    title,
    count,
    open = false,
    children,
  }: { title: string; count?: number; open?: boolean; children: Snippet } = $props()
</script>

<details class="group" {open}>
  <summary
    class="flex items-center gap-2 mb-2 cursor-pointer list-none text-lg font-semibold [&::-webkit-details-marker]:hidden"
  >
    <svg
      class="w-3 h-3 shrink-0 transition-transform group-open:rotate-90"
      viewBox="0 0 12 12"
      fill="currentColor"
      aria-hidden="true"
    >
      <path d="M3 1l6 5-6 5z" />
    </svg>
    {title}
    {#if count !== undefined}
      <span class="text-sm font-normal text-stone-500">({count})</span>
    {/if}
  </summary>
  {@render children()}
</details>
