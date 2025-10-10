#!/bin/bash
#===============================================================================
# test/ecbundle/test_ecbundle_setup.sh - Test ecbundle setup
#===============================================================================

set -e

echo "Testing ecbundle setup..."

# Check if ecbundle-create is available
if command -v ecbundle-create >/dev/null 2>&1; then
    echo "✓ ecbundle-create found: $(which ecbundle-create)"
else
    echo "✗ ecbundle-create not found in PATH"
    echo "  Please install ecbundle from https://github.com/ecmwf/ecbundle"
    exit 1
fi

# Check if ecbundle-build is available
if command -v ecbundle-build >/dev/null 2>&1; then
    echo "✓ ecbundle-build found: $(which ecbundle-build)"
else
    echo "✗ ecbundle-build not found in PATH"
    echo "  Please install ecbundle from https://github.com/ecmwf/ecbundle"
    exit 1
fi

# Check if bundle.yml exists
if [ -f "bundle.yml" ]; then
    echo "✓ bundle.yml found"
else
    echo "✗ bundle.yml not found in current directory"
    echo "  Please run this script from the FESOM2 root directory"
    exit 1
fi

# Test ecbundle-create help
echo "Testing ecbundle-create help..."
if ecbundle-create --help >/dev/null 2>&1; then
    echo "✓ ecbundle-create --help works"
else
    echo "✗ ecbundle-create --help failed"
fi

# Test ecbundle-build help
echo "Testing ecbundle-build help..."
if ecbundle-build --help >/dev/null 2>&1; then
    echo "✓ ecbundle-build --help works"
else
    echo "✗ ecbundle-build --help failed"
fi

echo ""
echo "ecbundle setup test completed successfully!"
echo "You can now run the ecbundle tests with:"
echo "  ctest -L ecbundle"
echo "  or"
echo "  make run_ecbundle_tests" 