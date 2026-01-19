@echo off
REM 自動同步腳本 (Windows版本) - 每隔指定時間自動 commit 並 push
REM 用法: scripts\auto-sync.bat [間隔秒數]
REM 預設間隔: 300 秒 (5分鐘)
REM 提示: 請使用 git-sync.bat 來控制啟動/停止

setlocal enabledelayedexpansion

set INTERVAL=%1
if "%INTERVAL%"=="" set INTERVAL=300

set "PROJECT_DIR=%~dp0.."
set "LOG_FILE=%PROJECT_DIR%\auto-sync.log"

cd /d "%PROJECT_DIR%"

echo ==========================================
echo   Git 自動同步已啟動
echo   同步間隔: %INTERVAL% 秒
echo   按 Ctrl+C 停止
echo ==========================================

:loop
REM 檢查是否有修改
git status --porcelain >nul 2>&1
if errorlevel 1 (
    echo [%date% %time:~0,8%] Git 命令錯誤
    timeout /t %INTERVAL% /nobreak >nul
    goto loop
)

for /f %%i in ('git status --porcelain') do (
    set HAS_CHANGES=1
    goto sync
)
set HAS_CHANGES=0

:sync
if "%HAS_CHANGES%"=="1" (
    echo.
    echo [%date% %time:~0,8%] 偵測到修改，開始同步...

    REM 暫存所有修改
    git add . >nul 2>&1

    REM 提交
    git commit -m "auto-sync: %date% %time:~0,8%" 2>&1 | findstr /v "^$"

    REM 推送
    git push 2>&1 | findstr /v "^$"

    echo [%date% %time:~0,8%] 同步完成!
) else (
    echo [%time:~0,8%] 無修改，等待中...
)

timeout /t %INTERVAL% /nobreak >nul
goto loop
