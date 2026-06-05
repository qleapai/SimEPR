#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3.11 -m streamlit run app.py --server.port 8502
