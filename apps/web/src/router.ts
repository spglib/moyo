import Home from './routes/Home.svelte'
import SpaceGroupList from './routes/SpaceGroupList.svelte'
import SpaceGroup from './routes/SpaceGroup.svelte'
import LayerGroupList from './routes/LayerGroupList.svelte'
import LayerGroup from './routes/LayerGroup.svelte'
import MagneticSpaceGroupList from './routes/MagneticSpaceGroupList.svelte'
import MagneticSpaceGroup from './routes/MagneticSpaceGroup.svelte'
import NotFound from './routes/NotFound.svelte'

export const routes = {
  '/': Home,
  '/sg': SpaceGroupList,
  '/sg/:number': SpaceGroup,
  '/lg': LayerGroupList,
  '/lg/:number': LayerGroup,
  '/msg': MagneticSpaceGroupList,
  '/msg/:uni_number': MagneticSpaceGroup,
  '*': NotFound,
}
