#!/bin/bash

# BASE_CMD="python main.py -d TOX -m DrugNCF -f 0"
BASE_CMD="python main.py -d TOX -m DrugNCF -f 2"

for k in 16 32 64 128 256 512 1024 2048
do
    echo "Running with -k: $k"
    $BASE_CMD -k $k
    if [ $? -ne 0 ]; then
        echo "Error with -k $k"
    fi
done

echo "Done!"
