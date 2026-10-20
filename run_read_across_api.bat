@echo off
cd /d "%~dp0"

echo.
echo ========================================
echo   PRIME-DS Read-Across Flask API
echo   Browser: http://localhost:5000
echo ========================================
echo.

pip show flask >nul 2>&1
if errorlevel 1 (
    echo Installing flask...
    pip install "flask>=3.0.0"
)

set PRIME_FLASK_PORT=5000
python read_across_flask_p00_app.py

if errorlevel 1 (
    echo.
    echo [오류] Flask 시작 실패.
    pause
)
