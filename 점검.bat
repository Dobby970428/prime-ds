@echo off
chcp 65001 >nul
title PRIME-DS 자체 점검

echo.
echo  ==============================================================
echo   PRIME-DS  ·  2026-09-03 수정본 자체 점검
echo  ==============================================================
echo.
echo   아무것도 고치지 않습니다. 돌려 보고 결과만 적습니다.
echo   1~3분쯤 걸립니다.
echo.

cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
    py -3.11 _selftest_20260903.py 2>nul || py _selftest_20260903.py
) else (
    python _selftest_20260903.py
)

if errorlevel 1 (
    echo.
    echo  ------------------------------------------------------------
    echo   점검 스크립트 자체가 실행되지 않았습니다.
    echo   위에 찍힌 오류를 그대로 Claude 에게 보여 주세요.
    echo  ------------------------------------------------------------
    echo.
) else (
    echo.
    echo  ------------------------------------------------------------
    echo   결과 파일 두 개가 이 폴더에 생겼습니다.
    echo     _selftest_결과.txt    ^(사람이 읽는 보고서^)
    echo     _selftest_결과.json   ^(Claude 가 읽는 원자료^)
    echo.
    echo   Claude 에게 "점검 끝났어" 라고만 말씀하시면
    echo   제가 이 파일들을 읽고 판단하겠습니다.
    echo  ------------------------------------------------------------
    echo.
)

pause
