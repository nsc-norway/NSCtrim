#!/bin/bash

VERSION=$( git describe --tags --dirty )
if [ -z "$VERSION" ]
then
    VERSION=latest
fi
docker build -t nsctrim:$VERSION --build-arg VERSION="$VERSION" .
echo "--- NSCtrim Docker image name: nsctrim:$VERSION ---"
