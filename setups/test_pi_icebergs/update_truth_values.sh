#!/bin/bash
# Script to update truth values in setup.yml based on CI test output
# Usage: 
#   ./update_truth_values.sh ci_output.txt [setup.yml]
#   ./update_truth_values.sh setup.yml ci_output.txt
#   ./update_truth_values.sh < ci_output.txt
#   ./update_truth_values.sh setup.yml "paste CI output here"

# Smart argument detection
if [ -z "$1" ]; then
    # No arguments, check stdin
    if [ ! -t 0 ]; then
        SETUP_FILE="setup.yml"
        INPUT=$(cat)
    else
        echo "Error: No input provided"
        echo "Usage: $0 ci_output.txt [setup.yml]"
        echo "   or: $0 setup.yml ci_output.txt"
        echo "   or: $0 < ci_output.txt"
        echo "   or: $0 setup.yml \"CI output text\""
        exit 1
    fi
elif [ -f "$1" ] && grep -q "Variable:.*current_value:" "$1" 2>/dev/null; then
    # First arg is a CI output file
    INPUT=$(cat "$1")
    SETUP_FILE="${2:-setup.yml}"
elif [ -f "$1" ]; then
    # First arg is setup file
    SETUP_FILE="$1"
    if [ -n "$2" ]; then
        if [ -f "$2" ]; then
            INPUT=$(cat "$2")
        else
            INPUT="$2"
        fi
    elif [ ! -t 0 ]; then
        INPUT=$(cat)
    else
        echo "Error: No CI output provided"
        echo "Usage: $0 setup.yml ci_output.txt"
        echo "   or: $0 setup.yml \"CI output text\""
        echo "   or: cat ci_output.txt | $0 setup.yml"
        exit 1
    fi
else
    echo "Error: File $1 not found"
    exit 1
fi

# Check if setup.yml exists
if [ ! -f "$SETUP_FILE" ]; then
    echo "Error: $SETUP_FILE not found"
    echo "Available setup directories:"
    find . -maxdepth 2 -name "setup.yml" -type f 2>/dev/null | sed 's|^\./||'
    exit 1
fi

# Parse and update values
declare -A updates
count=0

while IFS= read -r line; do
    # Match lines like: "Variable: a_ice, current_value: 0.19109711397780874, master_value: ..."
    if [[ "$line" =~ Variable:\ ([^,]+),\ current_value:\ ([^,]+), ]]; then
        var_name="${BASH_REMATCH[1]}"
        var_name=$(echo "$var_name" | xargs)  # trim whitespace
        current_val="${BASH_REMATCH[2]}"
        current_val=$(echo "$current_val" | xargs)  # trim whitespace
        
        updates["$var_name"]="$current_val"
        ((count++))
    fi
done <<< "$INPUT"

if [ $count -eq 0 ]; then
    echo "Error: No variable updates found in input"
    echo "Expected format: Variable: <name>, current_value: <value>, master_value: ..."
    rm "$BACKUP_FILE"
    exit 1
fi

echo "Found $count variable(s) to update:"
for var in "${!updates[@]}"; do
    echo "  $var: ${updates[$var]}"
done
echo ""

# Apply updates to setup.yml
for var in "${!updates[@]}"; do
    value="${updates[$var]}"
    # Escape special characters in value for sed
    escaped_value=$(echo "$value" | sed 's/[\/&]/\\&/g')
    sed -i "s/^\\([ ]*${var}:\\).*/\\1 ${escaped_value}/" "$SETUP_FILE"
done

echo "✓ Truth values updated successfully in $SETUP_FILE"
