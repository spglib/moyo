set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

build-python:
    maturin develop --release --manifest-path moyopy/Cargo.toml

install-python:
    python -m pip install uv maturin
    maturin develop --release --manifest-path moyopy/Cargo.toml
    python -m uv pip install -e "moyopy[dev]"
    pre-commit install

test-python:
    pytest -v moyopy/python/tests

build-python-docs:
    sphinx-autobuild moyopy/docs moyopy/docs/_build

# at moyo/
profile:
    CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --test test_moyo_dataset --root

doc:
    cargo doc --open

pre-commit:
    pre-commit run --all-files

pre-commit-all:
    pre-commit run --all-files --hook-stage manual

clean:
    rm moyopy/python/moyopy/*.so
