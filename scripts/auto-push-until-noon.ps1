# Auto Analysis and Push Script
# Run analysis and push every 30 minutes until noon

$repoPath = "c:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill"
$pythonExe = "C:/Users/88697.CHENPENGCHUNG12/D2Q9_PeriodicHill/.venv/Scripts/python.exe"
$analyzeScript = "c:/Users/88697.CHENPENGCHUNG12/D2Q9_PeriodicHill/analyze_velocity.py"
$intervalMinutes = 30
$stopHour = 12

Write-Host "========================================"
Write-Host "Auto Analysis and Push Script"
Write-Host "========================================"
Write-Host "Repo: $repoPath"
Write-Host "Interval: Every $intervalMinutes minutes"
Write-Host "Stop at: 12:00 noon"
$startTime = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
Write-Host "Start time: $startTime"
Write-Host "========================================"
Write-Host ""

Set-Location $repoPath

$cycleCount = 0

while ($true) {
    $now = Get-Date
    $currentHour = $now.Hour
    
    if ($currentHour -ge $stopHour) {
        Write-Host ""
        Write-Host "========================================"
        Write-Host "Reached 12:00, stopping automatic tasks"
        Write-Host "Total cycles: $cycleCount"
        $endTime = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        Write-Host "End time: $endTime"
        Write-Host "========================================"
        break
    }
    
    $cycleCount++
    $timeStr = $now.ToString("HH:mm:ss")
    Write-Host ""
    Write-Host "================================================"
    Write-Host "Cycle #$cycleCount - $timeStr"
    Write-Host "================================================"
    
    Write-Host ""
    Write-Host "[1/3] Running Python analysis script..."
    Write-Host "----------------------------------------------"
    
    try {
        & $pythonExe $analyzeScript
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  OK Python analysis completed"
        } else {
            Write-Host "  WARNING Python returned code: $LASTEXITCODE"
        }
    } catch {
        Write-Host "  ERROR Python execution failed"
    }
    
    Start-Sleep -Seconds 2
    
    Write-Host ""
    Write-Host "[2/3] Checking Git status..."
    Write-Host "----------------------------------------------"
    
    $status = git status --porcelain
    
    if ($status) {
        $changedFiles = ($status | Measure-Object).Count
        Write-Host "  Found $changedFiles changes"
        
        $status | ForEach-Object {
            Write-Host "    $_"
        }
        
        Write-Host ""
        Write-Host "[3/3] Committing and pushing..."
        Write-Host "----------------------------------------------"
        
        git add .
        
        $timestamp = $now.ToString("yyyy-MM-dd HH:mm:ss")
        
        $vtkPattern = "\.vtk"
        $figPattern = "figures.*\.png"
        
        $vtkMatches = $status | Select-String -Pattern $vtkPattern
        $figMatches = $status | Select-String -Pattern $figPattern
        
        $vtkCount = 0
        $figCount = 0
        
        if ($vtkMatches) {
            $vtkCount = ($vtkMatches | Measure-Object).Count
        }
        
        if ($figMatches) {
            $figCount = ($figMatches | Measure-Object).Count
        }
        
        if ($figCount -gt 0) {
            $commitMsg = "Update analysis figures at $timestamp"
        } elseif ($vtkCount -gt 0) {
            $commitMsg = "Update simulation data at $timestamp"
        } else {
            $commitMsg = "Auto-update at $timestamp"
        }
        
        Write-Host "  Commit: $commitMsg"
        git commit -m "$commitMsg"
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  OK Commit succeeded"
            
            Write-Host "  Pushing to remote..."
            git push
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host "  OK Push succeeded"
            } else {
                Write-Host "  ERROR Push failed"
            }
        } else {
            Write-Host "  ERROR Commit failed"
        }
    } else {
        Write-Host "  Working tree clean, nothing to push"
    }
    
    $nextRun = $now.AddMinutes($intervalMinutes)
    $nextStr = $nextRun.ToString("HH:mm:ss")
    Write-Host ""
    Write-Host "Next run: $nextStr"
    Write-Host "================================================"
    
    Write-Host ""
    Write-Host "Waiting $intervalMinutes minutes..."
    Start-Sleep -Seconds ($intervalMinutes * 60)
}
