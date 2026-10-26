#!/usr/bin/env bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

echo ""
echo "=========================================================================="
echo -e "\033[1;33m ⭐ Support Open Science: Please Star the Repository on GitHub!\033[0m"
echo -e "\033[4;34m 👉 https://github.com/Lukecele/spatial-massive-validation-study3\033[0m"
echo "=========================================================================="
echo ""
echo "🔬 Launching Massive Clinical Survival Validation Battery (N = 854)..."
echo ""

python3 "$DIR/src/run_survival_validation.py" "$@"
