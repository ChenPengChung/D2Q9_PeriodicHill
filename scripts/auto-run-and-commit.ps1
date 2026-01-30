param([int]$IntervalMinutes = 30)
$WorkspaceDir = 'C:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill'
Set-Location $WorkspaceDir

function Write-Log {
    param($Message)
    $timestamp = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    $logMessage = "[$timestamp] $Message"
    Write-Host $logMessage
    Add-Content -Path "$WorkspaceDir\auto_run.log" -Value $logMessage
}

function Test-Divergence {
    param($OutputFile)
    if (Test-Path $OutputFile) {
        $tail = Get-Content $OutputFile -Tail 50 -ErrorAction SilentlyContinue
        if ($tail -match 'nan|NaN|inf|Inf') { return $true }
    }
    return $false
}

function Commit-And-Push {
    param($Message)
    try {
        $added = git add output/*.vtk *.log 2>&1
        $status = git status --short
        if ($status) {
            git commit -m "$Message" | Out-Null
            git push | Out-Null
            Write-Log "已提交並推送: $Message"
        } else {
            Write-Log "沒有新變更"
        }
    } catch {
        Write-Log "Git 操作失敗: $_"
    }
}

Write-Log '=== 開始自動執行模擬 ==='
Write-Log "間隔時間: $IntervalMinutes 分鐘"

# 啟動模擬為背景作業
$job = Start-Job -ScriptBlock {
    Set-Location 'C:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill'
    & .\main.exe > simulation_output.txt 2>&1
}

Write-Log "模擬作業已啟動 (Job ID: $($job.Id))"
$iteration = 0
$startTime = Get-Date

try {
    while ($true) {
        Start-Sleep -Seconds ($IntervalMinutes * 60)
        $iteration++
        
        $jobState = Get-Job -Id $job.Id
        if ($jobState.State -ne 'Running') {
            Write-Log "模擬作業已結束 (狀態: $($jobState.State))"
            break
        }
        
        $runningHours = [math]::Round((Get-Date - $startTime).TotalHours, 2)
        Write-Log "--- 檢查點 $iteration (運行時間: $runningHours 小時) ---"
        
        if (Test-Divergence "$WorkspaceDir\simulation_output.txt") {
            Write-Log '!!! 檢測到發散 (NaN/Inf) !!!'
            Stop-Job -Id $job.Id
            Remove-Job -Id $job.Id
            Write-Log '已停止模擬作業'
            Commit-And-Push "發散於檢查點 $iteration - $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
            break
        }
        
        # 查看最新輸出
        if (Test-Path "$WorkspaceDir\simulation_output.txt") {
            $tail = Get-Content "$WorkspaceDir\simulation_output.txt" -Tail 5 -ErrorAction SilentlyContinue
            Write-Log "最新輸出: $($tail[-1])"
        }
        
        $commitMessage = "自動檢查點 $iteration - 運行 $runningHours h - $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
        Commit-And-Push $commitMessage
        
        Write-Log '繼續監控...'
    }
} catch {
    Write-Log "錯誤: $_"
} finally {
    $jobState = Get-Job -Id $job.Id -ErrorAction SilentlyContinue
    if ($jobState -and $jobState.State -eq 'Running') {
        Stop-Job -Id $job.Id
        Write-Log '已停止模擬作業'
    }
    Remove-Job -Id $job.Id -ErrorAction SilentlyContinue
    Write-Log '=== 自動執行結束 ==='
}
