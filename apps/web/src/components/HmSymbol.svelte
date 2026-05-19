<script lang="ts">
  import { parseHmShort } from '../lib/hmShort'

  let { symbol }: { symbol: string } = $props()
  const tokens = $derived(parseHmShort(symbol))
</script>

<span class="hm-symbol font-serif italic">
  {#each tokens as token, ti}
    {#if ti > 0}<span class="hm-gap"></span>{/if}<span class="hm-token"
      >{#each token as seg}{#if seg.kind === 'over'}<span class="hm-bar">{seg.text}</span
          >{:else if seg.kind === 'sub'}<sub class="hm-sub">{seg.text}</sub
          >{:else}{seg.text}{/if}{/each}</span
    >
  {/each}
</span>

<style>
  .hm-symbol {
    white-space: nowrap;
    font-feature-settings: 'tnum';
  }
  .hm-gap {
    display: inline-block;
    width: 0.25em;
  }
  .hm-bar {
    text-decoration: overline;
    text-decoration-thickness: 0.06em;
    text-underline-offset: 0.1em;
  }
  .hm-sub {
    font-size: 0.7em;
    line-height: 0;
    position: relative;
    vertical-align: baseline;
    bottom: -0.25em;
    margin-left: 0.02em;
  }
</style>
