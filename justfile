set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

build-python:
    maturin develop -m moyo_py/Cargo.toml

pre-commit:
    pre-commit run --all-files
