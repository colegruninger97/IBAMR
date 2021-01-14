#! /bin/bash
FILES="*.cpp *.h"
for file in $FILES; do
  echo $file
  clang-format -i -style=file $file
done
