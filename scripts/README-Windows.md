# Windows 自動同步設定說明

## 快速開始

### 方法一：背景自動同步（推薦）

1. **啟動自動同步**（預設每 5 分鐘同步一次）
   ```cmd
   scripts\git-sync.bat start
   ```

2. **自訂同步間隔**（例如每 60 秒）
   ```cmd
   scripts\git-sync.bat start 60
   ```

3. **查看運行狀態**
   ```cmd
   scripts\git-sync.bat status
   ```

4. **停止自動同步**
   ```cmd
   scripts\git-sync.bat stop
   ```

### 方法二：前台執行（會佔用終端視窗）

```cmd
scripts\auto-sync.bat 300
```

按 `Ctrl+C` 停止。

## 檔案說明

| 檔案 | 用途 | 平台 |
|------|------|------|
| `git-sync.bat` | 控制腳本（啟動/停止/狀態） | Windows |
| `auto-sync-worker.ps1` | PowerShell 背景工作腳本 | Windows |
| `auto-sync.bat` | 前台執行版本 | Windows |
| `git-sync.sh` | 控制腳本（啟動/停止/狀態） | Mac/Linux |
| `auto-sync.sh` | 前台執行版本 | Mac/Linux |

## 日誌檔案

所有同步記錄會寫入：
```
E:\D2Q9_PeriodicHill\auto-sync.log
```

查看最近記錄：
```cmd
powershell -Command "Get-Content auto-sync.log -Tail 20"
```

## 系統需求

- Git 已安裝並可在命令列使用
- PowerShell（Windows 內建）
- 已設定 GitHub 遠端倉庫認證

## 注意事項

1. **首次使用前**請確認 Git 認證已設定完成
2. **背景執行**會在關機前持續運行
3. 可在**開機時自動啟動**：
   - 建立捷徑到 `git-sync.bat start`
   - 放入 `%APPDATA%\Microsoft\Windows\Start Menu\Programs\Startup`

## 進階設定

### 設定開機自動啟動

1. 建立啟動腳本：
   ```cmd
   echo @echo off > "%APPDATA%\Microsoft\Windows\Start Menu\Programs\Startup\start-git-sync.bat"
   echo cd /d "E:\D2Q9_PeriodicHill" >> "%APPDATA%\Microsoft\Windows\Start Menu\Programs\Startup\start-git-sync.bat"
   echo scripts\git-sync.bat start 300 >> "%APPDATA%\Microsoft\Windows\Start Menu\Programs\Startup\start-git-sync.bat"
   ```

2. 下次開機會自動啟動

### 查看背景進程

```cmd
tasklist | findstr powershell
```

## 故障排除

### 問題：無法推送到 GitHub

**解決方法**：
```cmd
git config --global credential.helper wincred
git push
```
輸入 GitHub 帳號密碼一次，之後會自動記憶。

### 問題：PowerShell 執行政策錯誤

**解決方法**：
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

### 問題：腳本無法啟動

**檢查步驟**：
1. 確認當前目錄：`cd E:\D2Q9_PeriodicHill`
2. 測試 Git：`git status`
3. 手動執行：`scripts\auto-sync.bat 60`
