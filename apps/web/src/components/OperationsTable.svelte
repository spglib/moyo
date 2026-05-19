<script lang="ts">
  import type { MoyoOperation, MoyoMagneticOperation } from '@spglib/moyo-wasm'
  import { formatOperationXyz } from '../lib/format'

  type Op = MoyoOperation | MoyoMagneticOperation
  let { operations, hasTimeReversal = false }: { operations: Op[]; hasTimeReversal?: boolean } =
    $props()
</script>

<div class="overflow-x-auto rounded border border-slate-200 dark:border-slate-800">
  <table class="min-w-full text-sm font-mono">
    <thead class="bg-slate-50 dark:bg-slate-900 text-xs uppercase tracking-wide">
      <tr>
        <th class="px-2 py-1 text-left">#</th>
        <th class="px-2 py-1 text-left">Coordinate triplet</th>
        {#if hasTimeReversal}
          <th class="px-2 py-1 text-left">1'</th>
        {/if}
      </tr>
    </thead>
    <tbody>
      {#each operations as op, i}
        <tr class="border-t border-slate-100 dark:border-slate-800">
          <td class="px-2 py-1 align-top text-slate-500">{i + 1}</td>
          <td class="px-2 py-1 align-top">{formatOperationXyz(op.rotation, op.translation)}</td>
          {#if hasTimeReversal}
            <td class="px-2 py-1 align-top"
              >{(op as MoyoMagneticOperation).time_reversal ? '-1' : '+1'}</td
            >
          {/if}
        </tr>
      {/each}
    </tbody>
  </table>
</div>
<p class="mt-2 text-xs text-slate-500 dark:text-slate-400">
  {operations.length} operation{operations.length === 1 ? '' : 's'}. Note: the ordering of
  operations is not guaranteed to match the conventions used by ITA or the Bilbao Crystallographic
  Server.
</p>
