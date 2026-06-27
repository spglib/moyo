---
project: moyof
summary: Fortran interface to moyoc, the C bindings for moyo.
author: Kohei Shinohara
src_dir: ./fortran
exclude_dir: ./fortran/tests
output_dir: ./docs/fortran-api
project_github: https://github.com/spglib/moyo
display: public
source: false
graph: false
predocmark: >
md_extensions: markdown.extensions.toc
---

Fortran interface (`module moyo`) to **moyoc**, the C bindings for the
[moyo](https://github.com/spglib/moyo) crystal symmetry finder.

See the [C API reference](../api/index.html) for the underlying C functions, and
the [Migration from spglib](../migration/) guide for usage.
