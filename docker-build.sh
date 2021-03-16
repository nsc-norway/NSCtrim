#!/bin/bash

VERSION=$(git describe --tags --dirty)
docker build -t trimmer --build-arg VERSION="$VERSION" .
