import { mount } from 'svelte'
import App from './App.svelte'
import './app.css'
import {
  getAllSpaceGroups,
  getAllLayerGroups,
  getAllMagneticSpaceGroups,
} from './lib/catalogIndex'

const target = document.getElementById('app')
if (!target) throw new Error('#app element not found')

const app = mount(App, { target })

queueMicrotask(() => {
  void getAllSpaceGroups()
  void getAllLayerGroups()
  void getAllMagneticSpaceGroups()
})

export default app
