#!/bin/bash

# Compile Macro_all.m with all dependencies (ALL files except Git artifacts)
# Usage: ./compile_macro_all.sh

# Configuration
MATLAB_COMPILER="/opt/matlab/R2024b/bin/mcc"
SOURCE_FILE="Macro_all.m"
OUTPUT_DIR="compiled_output"
EXCLUDE_PATTERNS=(".git/*" "*.gitignore" "*.gitmodules")

echo "=== Starting compilation (including all non-Git files) ==="

# Clean previous compilation if output directory exists
if [ -d "$OUTPUT_DIR" ]; then
    echo "Removing existing output directory..."
    rm -rf "$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_DIR"

# Build find command for exclusion
find_exclude=()
for pattern in "${EXCLUDE_PATTERNS[@]}"; do
    find_exclude+=(-not -path "./$pattern")
done

# Collect all files except Git-related
echo "Scanning for files to include..."
ALL_FILES=()
while IFS= read -r -d $'\0' file; do
    ALL_FILES+=("-a" "$file")
done < <(find . -type f "${find_exclude[@]}" -print0)

# Get the user paths added during runtime
MCR_ARGS=()
while IFS= read -r path; do
    # Clean carriage returns (if any) and add to list
    clean_path=$(echo "$path" | tr -d '\r')
    MCR_ARGS+=("-a" "$clean_path")
done < Macro_all_paths.txt
echo "User paths added: ${#MCR_ARGS[@]}"

# Compile command
echo "Compiling with mcc (${#ALL_FILES[@]} files included)..."
$MATLAB_COMPILER -m "$SOURCE_FILE" \
    "${ALL_FILES[@]}" \
    "${MCR_ARGS[@]}" \
    -d "$OUTPUT_DIR" \
    -v \
    -R '-nodisplay' \
    -R '-nosplash'

# Check success
if [ $? -eq 0 ]; then
    echo -e "\n=== Success ==="
    echo "Executable: $OUTPUT_DIR/Macro_all"
    echo "Included files:"
    find "$OUTPUT_DIR" -type f | sed 's/^/  /'
else
    echo -e "\n!!! Compilation failed !!!"
    exit 1
fi

# Set executable permissions (Linux/Mac)
if [ "$OSTYPE" != "msys" ] && [ "$OSTYPE" != "win32" ]; then
    chmod +x "$OUTPUT_DIR/Macro_all"
fi

exit 0
