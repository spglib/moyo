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
    rm -r moyopy/docs/_build

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
py-docs:
    sphinx-autobuild moyopy/docs moyopy/docs/_build

################################################################################
# C
################################################################################

[group('c')]
[working-directory: 'moyoc']
c-build:
    make build

[group('c')]
[working-directory: 'moyoc']
c-test:
    make tests

[group('c')]
[working-directory: 'moyoc']
c-clean:
    make clean

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
web-install:
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
