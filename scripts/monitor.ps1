# 監控腳本 - 每30分鐘檢查並提交
param([int]$IntervalMinutes = 30, [string]$OutputFile = 'simulation_output.txt')

$WorkspaceDir = 'C:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill'
Set-Location $WorkspaceDir
$LogFile = "$WorkspaceDir\monitor.log"

function Write-Log {
    param($Message)
    $timestamp = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    $logMessage = "[$timestamp] $Message"
    Write-Host $logMessage
    Add-Content -Path $LogFile -Value $logMessage
}

function Test-Divergence {
    if (Test-Path $OutputFile) {
        $tail = Get-Content $OutputFile -Tail 100 -ErrorAction SilentlyContinue | Out-String
        return ($tail -match 'nan|NaN|inf|Inf')
    }
    return $false
}

function Commit-And-Push {
    param($Message)
    try {
        git add output/*.vtk *.log 2>&1 | Out-Null
        $status = git status --short
        if ($status) {
            git commit -m "$Message" 2>&1 | Out-Null
            git push 2>&1 | Out-Null
            Write-Log "已提交: $Message"
        } else {
            Write-Log "無變更"
        }
    } catch {
        Write-Log "Git錯誤: $_"
    }
}

Write-Log '=== 開始監控 ==='
Write-Log "檢查間隔: $IntervalMinutes 分鐘"

$iteration = 0
$startTime = Get-Date

while ($true) {
    Start-Sleep -Seconds ($IntervalMinutes * 60)
    $iteration++
    
    # 檢查 main.exe 是否還在運行
    $process = Get-Process -Name 'main' -ErrorAction SilentlyContinue
    if (-not $process) {
        Write-Log "main.exe 已停止"
        Commit-And-Push "程序結束於檢查點 $iteration"
        break
    }
    
    $hours = [math]::Round((Get-Date - $startTime).TotalHours, 2)
    Write-Log "--- 檢查點 $iteration ($hours 小時) ---"
    
    # 檢查發散
    if (Test-Divergence) {
        Write-Log "!!! 檢測到發散 !!!"
        Stop-Process -Name 'main' -Force
        Commit-And-Push "發散於檢查點 $iteration - $(Get-Date -Format 'MM-dd HH:mm')"
        break
    }
    
    # 顯示進度
    if (Test-Path $OutputFile) {
        $lastLine = Get-Content $OutputFile -Tail 1
        Write-Log "進度: $lastLine"
    }
    
    # 定期提交
    Commit-And-Push "檢查點$iteration - $hours h - $(Get-Date -Format 'MM-dd HH:mm')"
}

Write-Log '=== 監控結束 ==='
