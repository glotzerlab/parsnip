for version in 3.{9..13}; do
    uv pip compile --python-version="$version" pyproject.toml tests/requirements.in > ".github/requirements-$version.txt"
done
