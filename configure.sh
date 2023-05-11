#!/bin/bash
cmake -B build -S .
cmake --build build
cmake --graphviz==./doc/test.dot
