import { defineConfig } from 'vite'
import { svelte } from '@sveltejs/vite-plugin-svelte'
import tailwindcss from '@tailwindcss/vite'

// On GitHub Pages the site is served from /<repo>/.
// Override locally / per environment via VITE_BASE.
const base = process.env.VITE_BASE ?? '/'

export default defineConfig({
  base,
  plugins: [svelte(), tailwindcss()],
  test: {
    environment: 'jsdom',
    globals: true,
  },
})
