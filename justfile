set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

pre-commit:
    pre-commit run --all-files

pre-commit-all:
    pre-commit run --all-files --hook-stage manual

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
[working-directory: 'moyo']
profile:
    CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --test test_moyo_dataset --root

################################################################################
# Python
################################################################################

[group('python')]
py-build:
    maturin develop --release --manifest-path moyopy/Cargo.toml

[group('python')]
py-install: py-build
    python -m pip install -e "moyopy[dev]"

[group('python')]
py-test:
    python -m pytest -v moyopy/python/tests

[group('python')]
py-docs:
    sphinx-autobuild moyopy/docs moyopy/docs/_build
