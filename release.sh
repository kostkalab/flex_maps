#!/bin/bash
# Create a GitHub release for packaged flex_maps artifacts.

set -euo pipefail

RELEASE_ONLY=false

if [ "${1:-}" = "--release-only" ]; then
    RELEASE_ONLY=true
fi

VERSION=$(tr -d ' \n' < VERSION)
if [ -z "$VERSION" ]; then
    echo "ERROR: Could not read version from VERSION"
    exit 1
fi

if [ "$RELEASE_ONLY" = false ]; then
    python scripts/package_release.py
fi

if [ ! -f release/LATEST_RELEASE ]; then
    echo "ERROR: release/LATEST_RELEASE not found. Run scripts/package_release.py first."
    exit 1
fi

TIMESTAMP=$(tr -d ' \n' < release/LATEST_RELEASE)
TAG="v${VERSION}"

expected=(
    "release/flex_maps_maps.${TIMESTAMP}.tar.gz"
    "release/flex_maps_reaction_tables.${TIMESTAMP}.tar.gz"
    "release/flex_maps_reaction_embeddings.${TIMESTAMP}.tar.gz"
    "release/flex_maps_validation.${TIMESTAMP}.tar.gz"
)

for file in "${expected[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Missing release archive: $file"
        exit 1
    fi
done

if [ -n "$(git status --porcelain)" ]; then
    echo "ERROR: Git working tree is not clean. Commit or stash changes before release."
    git status --short
    exit 1
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    tag_commit=$(git rev-list -n 1 "$TAG")
    head_commit=$(git rev-parse HEAD)
    if [ "$tag_commit" != "$head_commit" ]; then
        echo "ERROR: Tag $TAG already exists and does not point to HEAD."
        exit 1
    fi
fi

if ! command -v gh >/dev/null 2>&1; then
    echo "ERROR: GitHub CLI (gh) is not installed."
    exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
    echo "ERROR: Not authenticated with GitHub. Run: gh auth login"
    exit 1
fi

echo "Release: $TAG"
echo "Timestamp: $TIMESTAMP"
echo ""
echo "Files to release:"
for file in "${expected[@]}"; do
    sha=$(python -c "import hashlib, pathlib; p=pathlib.Path('$file'); print(hashlib.sha256(p.read_bytes()).hexdigest())")
    size=$(python -c "import pathlib; print(pathlib.Path('$file').stat().st_size)")
    echo "$file  $size bytes  sha256=$sha"
done
echo ""

read -p "Continue with release $TAG? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Release cancelled."
    exit 0
fi

git tag "$TAG" 2>/dev/null || echo "Tag $TAG already exists, continuing."
git push origin "$TAG" 2>/dev/null || echo "Tag $TAG already pushed, continuing."

gh release create "$TAG" \
    --title "Flex Maps ${TAG}" \
    --notes "Cross-species SMILES-filtered metabolic maps and reaction embeddings.

## Contents
- SMILES-filtered GraphML maps and PDF reports
- Reaction tables with oriented I_to_O and O_to_I reaction strings
- Reaction embeddings: Morgan side means, DRFP, ReactionT5 4vec, RXNGraphormer 2vec
- Validation archive with manifest and cross-species sanity checks

Pipeline timestamp: ${TIMESTAMP}" \
    "${expected[@]}"

echo "Release $TAG created successfully."
