set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

build-python:
    maturin develop --release --manifest-path moyo_py/Cargo.toml

install-python:
    maturin develop --release --manifest-path moyo_py/Cargo.toml
    pip install -e "moyo_py[dev]"

test-python:
    pytest -v moyo_py/python/tests

pre-commit:
    pre-commit run --all-files

clean:
    rm moyo_py/python/moyo/*.so
