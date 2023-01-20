#!/bin/bash

ROOT=$(pwd)
BUILDDIR=release-build
PROJECT=$1

cd $PROJECT
rm -rf $BUILDDIR
cd ..

