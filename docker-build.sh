#!/bin/bash

VERSION=$(git describe --tags --dirty)
docker build -t trimmer:$VERSION --build-arg VERSION="$VERSION" .
