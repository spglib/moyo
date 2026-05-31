<script lang="ts">
  import type { Snippet } from 'svelte'

  let {
    title,
    count,
    open = false,
    children,
  }: { title: string; count?: number; open?: boolean; children: Snippet } = $props()
</script>

<details {open}>
  <summary class="flex items-center gap-2 mb-2 cursor-pointer text-lg font-semibold">
    <svg
      class="chevron w-3 h-3 shrink-0 transition-transform"
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

<style>
  /* Hide the native disclosure marker; we render our own chevron. */
  summary {
    list-style: none;
  }
  summary::-webkit-details-marker {
    display: none;
  }
  /* Child combinator keeps nested sections from rotating each other's chevron. */
  details[open] > summary .chevron {
    transform: rotate(90deg);
  }
</style>
