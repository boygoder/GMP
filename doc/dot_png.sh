#!/bin/bash

for file in ./*
do
   if [[ -f "$file" ]] && [[ "$file" == "./test.dot"* ]] && [[ "$file" != *.png ]]; then
        filename=$(basename -- "$file")
        extension="${filename##*.}"
        filename="${filename%.*}"
        echo $file
        dot -Tpng "$file" -o "./${filename}.png"
    fi
done

