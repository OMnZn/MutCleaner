#!/bin/bash
set -e

VERSION=$1
if [ -z "$VERSION" ]; then
    echo "Usage: ./tools/release.sh 0.5.0"
    exit 1
fi

echo "Starting release process for v$VERSION"
echo "========================================"

if [ -z "$GITHUB_TOKEN" ]; then
    echo "Error: Please set GITHUB_TOKEN"
    exit 1
fi

echo "[1/8] Running tests..."
pytest tests/ -v


echo "[2/8] Updating version to $VERSION..."
OLD_VERSION=$(grep 'version = ' pyproject.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')
sed -i "s/version = \"$OLD_VERSION\"/version = \"$VERSION\"/" pyproject.toml
sed -i "s/__version__ = \"$OLD_VERSION\"/__version__ = \"$VERSION\"/" mutcleaner/__init__.py
sed -i "s/release = \"$OLD_VERSION\"/release = \"$VERSION\"/" doc/source/conf.py

echo "[3/8] Generating CHANGELOG..."
LAST_TAG=$(git describe --tags --abbrev=0)

python tools/changelog.py \
    "$GITHUB_TOKEN" \
    "xulab-research/MutCleaner" \
    "$LAST_TAG..HEAD" \
    --template keepachangelog \
    --output "./doc/changelog/CHANGELOG_$1.md"

echo "[4/8] Building distribution..."
rm -rf dist/ build/
python -m build

echo "[5/8] Checking distribution..."
twine check dist/*

echo "[6/8] Committing changes..."
git add pyproject.toml mutcleaner/__init__.py doc/source/conf.py "doc/changelog/CHANGELOG_$VERSION.md"
git commit -m "chore: bump version to $VERSION"
git push origin main

echo "[7/8] Creating and pushing tag v$VERSION..."
git tag -a "v$VERSION" -m "Release v$VERSION"
git push origin "v$VERSION"


echo "[8/8] Uploading to PyPI..."
twine upload dist/*

echo "========================================"
echo "Release v$VERSION completed!"
echo ""
echo "Next steps:"
echo "1. Create GitHub Release at:"
echo "   https://github.com/xulab-research/MutCleaner/releases/new?tag=v$VERSION"
echo "2. Verify PyPI: https://pypi.org/project/mutcleaner/$VERSION/"
echo "3. Verify docs: https://xulab-research.github.io/MutCleaner/"