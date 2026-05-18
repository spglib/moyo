import Home from './routes/Home.svelte'
import SpaceGroup from './routes/SpaceGroup.svelte'
import LayerGroup from './routes/LayerGroup.svelte'
import MagneticSpaceGroup from './routes/MagneticSpaceGroup.svelte'
import NotFound from './routes/NotFound.svelte'

export const routes = {
  '/': Home,
  '/sg/:number': SpaceGroup,
  '/lg/:number': LayerGroup,
  '/msg/:uni_number': MagneticSpaceGroup,
  '*': NotFound,
}
