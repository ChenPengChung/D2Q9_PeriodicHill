@echo off
REM Git 自動同步控制腳本 (Windows版本)
REM 用法:
REM   scripts\git-sync.bat start [間隔秒數]  - 啟動自動同步 (預設300秒)
REM   scripts\git-sync.bat stop              - 停止自動同步
REM   scripts\git-sync.bat status            - 查看狀態

setlocal enabledelayedexpansion

set "PROJECT_DIR=%~dp0.."
set "PID_FILE=%PROJECT_DIR%\scripts\.sync.pid"
set "LOG_FILE=%PROJECT_DIR%\auto-sync.log"

cd /d "%PROJECT_DIR%"

if "%1"=="" goto status
if "%1"=="start" goto start
if "%1"=="stop" goto stop
if "%1"=="status" goto status
goto usage

:start
set INTERVAL=%2
if "%INTERVAL%"=="" set INTERVAL=300

REM 檢查是否已在運行
if exist "%PID_FILE%" (
    set /p OLD_PID=<"%PID_FILE%"
    tasklist /FI "PID eq !OLD_PID!" 2>nul | find "!OLD_PID!" >nul
    if !errorlevel! equ 0 (
        echo ⚠️  自動同步已在運行中 (PID: !OLD_PID!)
        exit /b 1
    )
)

REM 使用 PowerShell 啟動背景任務
powershell -Command "Start-Process -WindowStyle Hidden -FilePath 'powershell.exe' -ArgumentList '-ExecutionPolicy Bypass -File \"%PROJECT_DIR%\scripts\auto-sync-worker.ps1\" -Interval %INTERVAL%'"

timeout /t 2 /nobreak >nul

REM 取得最新的 PowerShell 進程 PID
for /f "tokens=2" %%a in ('tasklist /fi "imagename eq powershell.exe" /fo list ^| findstr /i "PID"') do (
    set "NEW_PID=%%a"
)

echo !NEW_PID!>"%PID_FILE%"
echo ✅ 自動同步已啟動 (PID: !NEW_PID!, 間隔: %INTERVAL%秒)
echo    日誌: %LOG_FILE%
exit /b 0

:stop
if not exist "%PID_FILE%" (
    echo ℹ️  自動同步未在運行
    exit /b 0
)

set /p PID=<"%PID_FILE%"
tasklist /FI "PID eq %PID%" 2>nul | find "%PID%" >nul
if %errorlevel% equ 0 (
    taskkill /PID %PID% /F >nul 2>&1
    del /q "%PID_FILE%" 2>nul
    echo ✅ 自動同步已停止 (PID: %PID%)
) else (
    del /q "%PID_FILE%" 2>nul
    echo ⚠️  進程已不存在，已清理 PID 檔案
)
exit /b 0

:status
if exist "%PID_FILE%" (
    set /p PID=<"%PID_FILE%"
    tasklist /FI "PID eq !PID!" 2>nul | find "!PID!" >nul
    if !errorlevel! equ 0 (
        echo ✅ 自動同步運行中 (PID: !PID!)
        echo.
        echo 最近日誌:
        powershell -Command "Get-Content '%LOG_FILE%' -Tail 5 -ErrorAction SilentlyContinue"
    ) else (
        echo ❌ 自動同步未運行
    )
) else (
    echo ❌ 自動同步未運行
)
exit /b 0

:usage
echo 用法: %~nx0 {start^|stop^|status} [間隔秒數]
echo.
echo   start [秒數]  啟動自動同步 (預設 300 秒)
echo   stop          停止自動同步
echo   status        查看運行狀態
exit /b 1
