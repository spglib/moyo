<script lang="ts">
  let { message }: { message: string } = $props()
</script>

<div class="rounded border border-red-300 bg-red-50 text-red-900 dark:bg-red-950/40 dark:text-red-200 dark:border-red-900 p-3 text-sm">
  {message}
</div>
