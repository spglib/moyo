<script lang="ts">
  import { formatRotationRow, formatFraction } from '../lib/format'

  type Mat9 = [number, number, number, number, number, number, number, number, number]
  type Vec3 = [number, number, number]
  interface Op {
    rotation: Mat9
    translation: Vec3
    time_reversal?: boolean
  }

  let { operations, hasTimeReversal = false }: { operations: Op[]; hasTimeReversal?: boolean } =
    $props()
</script>

<div class="overflow-x-auto rounded border border-slate-200 dark:border-slate-800">
  <table class="min-w-full text-sm font-mono">
    <thead class="bg-slate-50 dark:bg-slate-900 text-xs uppercase tracking-wide">
      <tr>
        <th class="px-2 py-1 text-left">#</th>
        <th class="px-2 py-1 text-left">Rotation</th>
        <th class="px-2 py-1 text-left">Translation</th>
        {#if hasTimeReversal}
          <th class="px-2 py-1 text-left">1'</th>
        {/if}
      </tr>
    </thead>
    <tbody>
      {#each operations as op, i}
        <tr class="border-t border-slate-100 dark:border-slate-800">
          <td class="px-2 py-1 align-top text-slate-500">{i + 1}</td>
          <td class="px-2 py-1 align-top whitespace-pre">
            {formatRotationRow(op.rotation, 0)}
            {'\n'}{formatRotationRow(op.rotation, 1)}
            {'\n'}{formatRotationRow(op.rotation, 2)}
          </td>
          <td class="px-2 py-1 align-top">
            ({formatFraction(op.translation[0])}, {formatFraction(op.translation[1])}, {formatFraction(
              op.translation[2]
            )})
          </td>
          {#if hasTimeReversal}
            <td class="px-2 py-1 align-top">{op.time_reversal ? '-1' : '+1'}</td>
          {/if}
        </tr>
      {/each}
    </tbody>
  </table>
</div>
<p class="mt-2 text-xs text-slate-500 dark:text-slate-400">
  {operations.length} operation{operations.length === 1 ? '' : 's'}.
</p>
