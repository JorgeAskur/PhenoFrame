#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_SO="${ROOT}/libmaize_c_api.so"
CXX="${CXX:-g++}"

echo "Building ${OUT_SO}"
"${CXX}" -shared -fPIC -O2 -std=c++17 -DMAIZE_C_API_BUILD \
    -I"${ROOT}" \
    "${ROOT}/maize_c_api.cpp" \
    "${ROOT}/Maize.cpp" \
    "${ROOT}/Descriptor.cpp" \
    "${ROOT}/vect3d.cpp" \
    "${ROOT}/tinyxml2.cpp" \
    -o "${OUT_SO}"

if [[ ! -f "${OUT_SO}" ]]; then
    echo "Build failed: expected output not created: ${OUT_SO}" >&2
    exit 1
fi

echo "Build complete: ${OUT_SO}"
