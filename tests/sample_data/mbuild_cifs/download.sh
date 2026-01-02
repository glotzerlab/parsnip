#!/usr/bin/env sh

SITE="https://api.github.com/repos/mosdef-hub/mbuild/contents/mbuild/utils/reference"

curl -s $SITE \
| jq -r '.[] 
         | select(.name | endswith(".cif")) 
         | .download_url' \
| xargs -n 1 curl -O
