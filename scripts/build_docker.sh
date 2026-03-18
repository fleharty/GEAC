#!/usr/bin/env bash
# Build and optionally push the geac container image to Google Container Registry.
# Uses Podman (docker-compatible, no daemon required).
#
# Usage:
#   bash scripts/build_docker.sh <GCP_PROJECT_ID> [--push]
#
# Examples:
#   bash scripts/build_docker.sh my-gcp-project
#   bash scripts/build_docker.sh my-gcp-project --push
#
# Prerequisites for pushing to GCR:
#   gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(tr -d '[:space:]' < "${REPO_ROOT}/VERSION")"

GCP_PROJECT="${1:?Usage: $0 <GCP_PROJECT_ID> [--push]}"
PUSH="${2:-}"

IMAGE="gcr.io/${GCP_PROJECT}/geac"

echo "==> Building geac ${VERSION} from ${REPO_ROOT}"
podman build \
    -f "${REPO_ROOT}/docker/Dockerfile" \
    -t "${IMAGE}:${VERSION}" \
    -t "${IMAGE}:latest" \
    "${REPO_ROOT}"

echo ""
echo "Built:"
echo "  ${IMAGE}:${VERSION}"
echo "  ${IMAGE}:latest"

if [[ "${PUSH}" == "--push" ]]; then
    echo ""
    echo "==> Pushing to GCR..."
    podman push "${IMAGE}:${VERSION}"
    podman push "${IMAGE}:latest"
    echo ""
    echo "Pushed:"
    echo "  ${IMAGE}:${VERSION}"
    echo "  ${IMAGE}:latest"
else
    echo ""
    echo "To push to GCR, first authenticate:"
    echo "  gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io"
    echo ""
    echo "Then push:"
    echo "  podman push ${IMAGE}:${VERSION}"
    echo "  podman push ${IMAGE}:latest"
    echo ""
    echo "Or re-run with --push (after logging in):"
    echo "  bash scripts/build_docker.sh ${GCP_PROJECT} --push"
fi
