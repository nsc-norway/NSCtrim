#!/bin/bash

VERSION=$( git describe --tags --dirty )
docker build -t nsctrim:$VERSION --build-arg VERSION="$VERSION" .
echo "--- NSCtrim Docker image name: nsctrim:$VERSION ---"