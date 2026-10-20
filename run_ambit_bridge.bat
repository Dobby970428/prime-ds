@echo off

cd /d "%~dp0"

echo.

echo ========================================

echo   PRIME-DS AMBIT2 Bridge Server

echo   URL: http://localhost:8080/ambit2

echo ========================================

echo.

echo PRIME-DS Streamlit(run_app.bat)과 별도 창에서 실행하세요.

echo 이 창을 닫거나 Ctrl+C 하면 Connected 상태가 proxy로 돌아갑니다.

echo.

echo Stopping old bridge on port 8080...

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do taskkill /F /PID %%a >nul 2>&1

echo Starting bridge... (Ctrl+C to stop)

echo.

python p55_ambit_bridge_server.py

