set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

prek:
    prek run --all-files

prek-all:
    prek run --all-files --hook-stage manual

clean:
    rm moyopy/python/moyopy/*.so
    rm -rf moyopy/site

################################################################################
# Rust
################################################################################

[group('rust')]
test:
    cargo test

[group('rust')]
doc:
    cargo doc --open

[group('rust')]
upgrade:
    cargo upgrade -i

[group('rust')]
[working-directory: 'moyo']
profile:
    CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --test test_moyo_dataset --root

################################################################################
# Python
################################################################################

[group('python')]
py-build:
    uv run maturin develop --release --manifest-path moyopy/Cargo.toml

[group('python')]
py-install:
    uv sync --all-extras --directory moyopy

[group('python')]
py-test:
    python -m pytest -v moyopy/python/tests

[group('python')]
[working-directory: 'moyopy']
py-docs:
    uv run zensical serve

################################################################################
# C
################################################################################

[group('c')]
[working-directory: 'moyoc']
c-build:
    cmake -S . -B build
    cmake --build build

[group('c')]
[working-directory: 'moyoc']
c-test: c-build
    ctest --test-dir build --output-on-failure

[group('c')]
[working-directory: 'moyoc']
c-docs:
    cbindgen --config cbindgen.toml --output build/moyoc.h
    doxygen Doxyfile
    uvx ford ford.md
    uvx zensical serve

[group('c')]
[working-directory: 'moyoc']
c-clean:
    rm -rf build

################################################################################
# JS (WASM)
################################################################################

[group('js')]
[working-directory: 'moyo-wasm']
js-install:
    npm install

[group('js')]
[working-directory: 'moyo-wasm']
js-build:
    npm run build

[group('js')]
[working-directory: 'moyo-wasm']
js-test:
    npm test

################################################################################
# Web app (apps/web)
################################################################################

[group('web')]
[working-directory: 'apps/web']
web-install: js-build
    npm install

[group('web')]
[working-directory: 'apps/web']
web-dev: js-build
    npm run dev

[group('web')]
[working-directory: 'apps/web']
web-build:
    npm run build

[group('web')]
[working-directory: 'apps/web']
web-preview: web-build
    npm run preview

[group('web')]
[working-directory: 'apps/web']
web-test:
    npm test

[group('web')]
[working-directory: 'apps/web']
web-check:
    npm run check

################################################################################
# Landing hub (apps/landing)
################################################################################

[group('landing')]
[working-directory: 'apps/landing']
landing-build:
    uv run zensical build --clean --strict

################################################################################
# All documents (apps/landing + moyopy + moyoc)
#
# `docs` builds every doc site, assembles them under apps/landing/site in the
# GitHub Pages layout (root = landing, /python = moyopy, /c = moyoc), and serves
# the result at http://localhost:8000. Cross-interface switcher links resolve to
# the production URLs; the /viewer/ path is not built here.
################################################################################

# Build all doc sites and assemble them under apps/landing/site.
[group('docs')]
docs-build:
    cd moyopy && uv run zensical build --clean --strict
    cd apps/landing && uv run zensical build --clean --strict
    cd moyoc && cbindgen --config cbindgen.toml --output build/moyoc.h && doxygen Doxyfile && uvx ford ford.md && uvx zensical build --clean --strict
    rm -rf apps/landing/site/python apps/landing/site/c
    cp -r moyopy/site apps/landing/site/python
    cp -r moyoc/site apps/landing/site/c

# Build all docs, then serve the assembled hub at http://localhost:8000.
[group('docs')]
docs: docs-build
    python3 -m http.server -d apps/landing/site 8000
