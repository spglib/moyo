<script lang="ts">
  import { link } from 'svelte-spa-router'
</script>

<section class="py-10 text-center">
  <h1 class="text-2xl font-semibold">Not found</h1>
  <p class="mt-2 text-stone-600 dark:text-stone-400">
    That route does not exist. <a use:link href="/" class="underline">Back to home</a>.
  </p>
</section>
