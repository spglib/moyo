set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

install-python:
    rm moyo_py/python/moyo/*.so
    maturin develop --release --manifest-path moyo_py/Cargo.toml
    pip install -e moyo_py

test-python:
    pytest -v moyo_py/python/tests

pre-commit:
    pre-commit run --all-files
