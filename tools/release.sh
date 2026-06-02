#!/bin/bash
set -euo pipefail

VERSION=${1:-}

if [ -z "$VERSION" ]; then
    echo "Usage: ./tools/release.sh 0.3.0"
    exit 1
fi

if [ -z "${GITHUB_TOKEN:-}" ]; then
    echo "Error: Please set GITHUB_TOKEN"
    exit 1
fi

if [ "$(git branch --show-current)" != "main" ]; then
    echo "Error: release must be run from main branch"
    exit 1
fi

if [ -n "$(git status --porcelain)" ]; then
    echo "Error: working tree is not clean"
    git status --short
    exit 1
fi

echo "Starting release v$VERSION"
echo "========================================"

echo "[1/9] Pulling latest main..."
git pull --ff-only origin main

echo "[2/9] Running tests..."
pytest tests/ -v

echo "[3/9] Updating version..."
OLD_VERSION=$(grep 'version = ' pyproject.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')

sed -i "s/version = \"$OLD_VERSION\"/version = \"$VERSION\"/" pyproject.toml
sed -i "s/__version__ = \"$OLD_VERSION\"/__version__ = \"$VERSION\"/" mutcleaner/__init__.py
sed -i "s/release = \"$OLD_VERSION\"/release = \"$VERSION\"/" doc/source/conf.py

echo "[4/9] Generating changelog..."
LAST_TAG=$(git describe --tags --abbrev=0)
CHANGELOG="doc/changelog/CHANGELOG_$VERSION.md"

python tools/changelog.py \
    "$GITHUB_TOKEN" \
    "xulab-research/MutCleaner" \
    "$LAST_TAG..HEAD" \
    --template keepachangelog \
    --output "$CHANGELOG"

sed -i "s/## \[HEAD\]/## [$VERSION]/" "$CHANGELOG"
sed -i "s/$LAST_TAG..HEAD/$LAST_TAG..v$VERSION/g" "$CHANGELOG"

echo "[5/9] Committing release changes..."
git add pyproject.toml mutcleaner/__init__.py doc/source/conf.py "$CHANGELOG"
git commit -m "chore: release v$VERSION"

echo "[6/9] Creating local tag..."
git tag -a "v$VERSION" -m "Release v$VERSION"

echo "[7/9] Building distribution..."
rm -rf dist/ build/ *.egg-info
python -m build

python - <<EOF
from pathlib import Path

version = "$VERSION"
files = [p.name for p in Path("dist").iterdir()]
if not files or any(version not in f for f in files):
    raise SystemExit(f"Unexpected dist files: {files}")
print("\\n".join(files))
EOF

echo "[8/9] Checking distribution..."
twine check dist/*

echo "[9/9] Pushing and uploading..."
git push origin main
git push origin "v$VERSION"
twine upload dist/*

echo "========================================"
echo "Release v$VERSION completed!"
echo ""
echo "Next steps:"
echo "1. Create GitHub Release:"
echo "   https://github.com/xulab-research/MutCleaner/releases/new?tag=v$VERSION"
echo "2. Verify PyPI:"
echo "   https://pypi.org/project/mutcleaner/$VERSION/"
echo "3. Verify docs:"
echo "   https://xulab-research.github.io/MutCleaner/"