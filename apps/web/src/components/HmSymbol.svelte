<script lang="ts">
  import katex from 'katex'
  import 'katex/dist/katex.min.css'
  import { hmShortToLatex } from '../lib/hmShort'

  let { symbol }: { symbol: string } = $props()
  const html = $derived(
    katex.renderToString(hmShortToLatex(symbol), {
      throwOnError: false,
      displayMode: false,
      output: 'html',
    })
  )
</script>

<span class="hm-symbol">{@html html}</span>

<style>
  .hm-symbol {
    white-space: nowrap;
  }
  .hm-symbol :global(.katex) {
    font-size: 1em;
  }
</style>
