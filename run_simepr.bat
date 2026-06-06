@echo off
cd /d "%~dp0"
echo Starting SimEPR (E-drive venv)...
".venv\Scripts\streamlit.exe" run app.py --server.port 8502