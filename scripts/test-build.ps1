# PowerShell 编译测试脚本
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  D2Q9 PeriodicHill 编译测试" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# 检查编译器
Write-Host "[1/4] 检查 g++ 编译器..." -ForegroundColor Yellow
$gccVersion = g++ --version 2>&1 | Select-Object -First 1
if ($LASTEXITCODE -eq 0) {
    Write-Host "  ✓ g++ 找到: $gccVersion" -ForegroundColor Green
} else {
    Write-Host "  ✗ 找不到 g++ 编译器！" -ForegroundColor Red
    Write-Host "  请确认 C:\msys64\mingw64\bin 在 PATH 中" -ForegroundColor Yellow
    exit 1
}

# 清理旧文件
Write-Host ""
Write-Host "[2/4] 清理旧的编译文件..." -ForegroundColor Yellow
Remove-Item -Path "hill.exe" -ErrorAction SilentlyContinue
if (Test-Path "hill.exe") {
    Write-Host "  ✗ 无法删除 hill.exe（文件可能正在使用）" -ForegroundColor Red
} else {
    Write-Host "  ✓ 清理完成" -ForegroundColor Green
}

# Debug 编译
Write-Host ""
Write-Host "[3/4] Debug 模式编译..." -ForegroundColor Yellow
g++ -std=c++17 -g -Wall -Wextra -o hill.exe main.cpp -lm
if ($LASTEXITCODE -eq 0 -and (Test-Path "hill.exe")) {
    $fileSize = (Get-Item "hill.exe").Length / 1KB
    Write-Host "  ✓ Debug 编译成功！文件大小: $([math]::Round($fileSize, 2)) KB" -ForegroundColor Green
} else {
    Write-Host "  ✗ Debug 编译失败！" -ForegroundColor Red
    exit 1
}

# Release 编译
Write-Host ""
Write-Host "[4/4] Release 模式编译..." -ForegroundColor Yellow
Remove-Item -Path "hill.exe" -ErrorAction SilentlyContinue
g++ -std=c++17 -O3 -Wall -o hill.exe main.cpp -lm
if ($LASTEXITCODE -eq 0 -and (Test-Path "hill.exe")) {
    $fileSize = (Get-Item "hill.exe").Length / 1KB
    Write-Host "  ✓ Release 编译成功！文件大小: $([math]::Round($fileSize, 2)) KB" -ForegroundColor Green
} else {
    Write-Host "  ✗ Release 编译失败！" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  所有测试通过！" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "使用方法:" -ForegroundColor Yellow
Write-Host "  1. 在 VSCode 中按 Ctrl+Shift+B 编译" -ForegroundColor White
Write-Host "  2. 在 VSCode 中按 F5 开始调试" -ForegroundColor White
Write-Host "  3. 或在终端执行: .\hill.exe" -ForegroundColor White
