# Git 自動同步背景工作腳本 (Windows PowerShell)
param(
    [int]$Interval = 300
)

$ProjectDir = Split-Path -Parent $PSScriptRoot
$LogFile = Join-Path $ProjectDir "auto-sync.log"

Set-Location $ProjectDir

function Write-Log {
    param([string]$Message)
    $Timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $LogMessage = "[$Timestamp] $Message"
    Add-Content -Path $LogFile -Value $LogMessage
    Write-Host $LogMessage
}

Write-Log "=========================================="
Write-Log "  Git 自動同步已啟動 (PID: $PID)"
Write-Log "  同步間隔: $Interval 秒"
Write-Log "=========================================="

while ($true) {
    try {
        # 檢查是否有修改
        $status = git status --porcelain 2>&1

        if ($status) {
            $Timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
            Write-Log ""
            Write-Log "偵測到修改，開始同步..."

            # 暫存所有修改
            git add . 2>&1 | Out-Null

            # 提交
            $commitOutput = git commit -m "auto-sync: $Timestamp" 2>&1
            Write-Log $commitOutput

            # 推送
            $pushOutput = git push 2>&1
            Write-Log $pushOutput

            Write-Log "同步完成!"
        }
    }
    catch {
        Write-Log "錯誤: $_"
    }

    Start-Sleep -Seconds $Interval
}
