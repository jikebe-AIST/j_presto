#!/bin/bash

dirs=("src_GEprep" "src_Ens_Ana" "src_PCAaxis" "src_PCAproj" "src_distrib" "src_md_run")

for dir in "${dirs[@]}"; do
    if [ -d "$dir" ]; then
        echo "Cleaning $dir"
        (cd "$dir" && make clean)
    else
        echo "Warning: Directory $dir not found"
    fi
done

if [ -f "sp/calc_Efield.so" ]; then
    echo "Removing sp/calc_Efield.so"
    rm "sp/calc_Efield.so"
else
    echo "Warning: sp/calc_Efield.so not found"
fi

if [ -n "$J_PRESTO_PATH" ]; then
    if [ -d "$J_PRESTO_PATH" ]; then
        echo "Removing directory $J_PRESTO_PATH"
        rm -rf "$J_PRESTO_PATH"
    else
        echo "Warning: Directory $J_PRESTO_PATH not found"
    fi
else
    echo "Error: J_PRESTO_PATH is not set"
fi

echo "Uninstallation complete."

