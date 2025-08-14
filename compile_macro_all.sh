#!/bin/bash

# Compile Macro_all.m with all dependencies (ALL files except Git artifacts)
# Usage: ./compile_macro_all.sh

# Configuration
MATLAB_COMPILER="/opt/matlab/R2024b/bin/mcc"
SOURCE_FILE="Macro_all.m"
OUTPUT_DIR="compiled_matlab"
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

# Collect all folders/files except Git-related
echo "Scanning for files to include..."
ALL_DIRS=()
ALL_FILES=()
while IFS= read -r -d '' file; do
    ext="${file##*.}"
    if [[ "$ext" == "m" ]]; then        # For m-files -> add their parent folder via -I
        ALL_DIRS+=("$(dirname "$file")")
    fi
    ALL_FILES+=("-a" "$file")
done < <(find "$(dirname "$0")" \
              /home/common/matlab/sepia/external/MRI_susceptibility_calculation/MRI_susceptibility_calculation_20190912 \
              /home/common/matlab/sepia/external/SEGUE/SEGUE_28012021 \
              -type f "${find_exclude[@]}" -print0)

# Remove duplicate directories for -I
mapfile -t UNIQUE_DIRS < <(printf '%s\n' "${ALL_DIRS[@]}" | sort -u)

# Build -I args
INCLUDE_DIRS=()
for dir in "${UNIQUE_DIRS[@]}"; do
    INCLUDE_DIRS+=("-I" "$dir")
done

# Compile command
echo "Compiling with mcc (${#ALL_DIRS[@]} folders and ${#ALL_FILES[@]} files included)..."
$MATLAB_COMPILER -m "$SOURCE_FILE" \
    "${INCLUDE_DIRS[@]}" \
    "${ALL_FILES[@]}" \
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
