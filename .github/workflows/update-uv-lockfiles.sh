#!/bin/bash

for version in 3.{8..13}; do
    uv pip compile --python-version="$version" pyproject.toml requirements-sympy.txt tests/requirements.in > ".github/requirements-$version.txt"
done
