#!/bin/bash
#=============================================================================
# run_analysis.sh - 跨平台執行分析腳本
# 自動偵測並使用虛擬環境
#=============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="$SCRIPT_DIR/analyze_velocity.py"

# 清理 macOS 資源分支檔案（如果存在）
if [[ "$OSTYPE" == "darwin"* ]]; then
    find "$PARENT_DIR/.venv" -name "._*" -delete 2>/dev/null
    find "$PARENT_DIR/.venv" -name ".__*" -delete 2>/dev/null
fi

# 嘗試找到 Python
if [ -f "$PARENT_DIR/.venv/bin/python" ]; then
    PYTHON_BIN="$PARENT_DIR/.venv/bin/python"
elif [ -f "$PARENT_DIR/.venv/Scripts/python.exe" ]; then
    PYTHON_BIN="$PARENT_DIR/.venv/Scripts/python.exe"
elif command -v python3 &> /dev/null; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "Using Python: $PYTHON_BIN"
echo "Running: $PYTHON_SCRIPT"
echo "=============================================="

"$PYTHON_BIN" "$PYTHON_SCRIPT" "$@"
