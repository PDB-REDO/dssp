#!/bin/bash

set -e

file="$1"

echo -n "[["
while IFS= read -r line; do
	echo $line | sed -e 's/\[/@<:@/g' -e 's/\]/@:>@/g' -e 's/#/@%:@/g' -e 's/\$/@S|@/g'
	echo
done < "$file"

echo -n "]]"
