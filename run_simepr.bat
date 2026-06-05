@echo off
cd /d "%~dp0"
py -3.11 -m streamlit run app.py --server.port 8502
