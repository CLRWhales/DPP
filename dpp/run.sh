#!/bin/bash
# ================================================================
# Script to CREATE .ini files for multiple dates, 
# based on example.ini
# ================================================================

TEMPLATE="example.ini"
CONFIG_DIR="./configs"
START_DATE=20250812
END_DATE=20250821

mkdir -p "$CONFIG_DIR"

current_date="$START_DATE"
while [ "$current_date" -le "$END_DATE" ]; do
    echo "=== Creation of the .ini file for $current_date ==="

    config_file="$CONFIG_DIR/config_${current_date}.ini"

    cp "$TEMPLATE" "$config_file"

    sed -i "s|/raw_data/[0-9]\{8\}/|/raw_data/${current_date}/|g" "$config_file"
    sed -i "s|/cleaning/[0-9]\{8\}/|/cleaning/${current_date}/|g" "$config_file" 2>/dev/null || true
    sed -i "s/^run_name *= *.*/run_name = run_${current_date}/" "$config_file"

    echo "Created : $config_file"
    echo ""

    current_date=$(date -d "$current_date +1 day" +"%Y%m%d")
done

echo "=== All .ini files have been created in $CONFIG_DIR ==="
