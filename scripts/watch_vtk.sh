#!/bin/bash
#=============================================================================
# watch_vtk.sh - 監控 VTK 檔案並自動執行視覺化
# 使用方式：./scripts/watch_vtk.sh &
#=============================================================================

# 設定路徑
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"
PYTHON_SCRIPT="$SCRIPT_DIR/analyze_velocity.py"
LOG_FILE="$SCRIPT_DIR/output/watch_vtk.log"
PYTHON_BIN="/Volumes/Seagate/D2Q9_PeriodicHill/.venv/bin/python"

# 記錄已處理的檔案
PROCESSED_FILE="$SCRIPT_DIR/output/.processed_vtk"

# 初始化
echo "===== VTK Watcher Started at $(date) =====" >> "$LOG_FILE"
echo "Monitoring: $OUTPUT_DIR" >> "$LOG_FILE"
echo "Python script: $PYTHON_SCRIPT" >> "$LOG_FILE"

# 建立已處理檔案清單（如果不存在）
touch "$PROCESSED_FILE"

# 主監控迴圈
while true; do
    # 找出所有 VTK 檔案
    for vtk_file in "$OUTPUT_DIR"/flow_*.vtk; do
        # 跳過如果沒有符合的檔案
        [ -e "$vtk_file" ] || continue
        
        # 取得檔案基本名稱
        basename_vtk=$(basename "$vtk_file")
        
        # 檢查是否已處理過
        if ! grep -q "^$basename_vtk$" "$PROCESSED_FILE" 2>/dev/null; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] New VTK found: $basename_vtk" >> "$LOG_FILE"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] New VTK found: $basename_vtk"
            
            # 執行 Python 視覺化腳本
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running visualization..." >> "$LOG_FILE"
            "$PYTHON_BIN" "$PYTHON_SCRIPT" >> "$LOG_FILE" 2>&1
            
            if [ $? -eq 0 ]; then
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Visualization completed successfully" >> "$LOG_FILE"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Visualization completed successfully"
            else
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Visualization failed!" >> "$LOG_FILE"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Visualization failed!"
            fi
            
            # 標記為已處理
            echo "$basename_vtk" >> "$PROCESSED_FILE"
        fi
    done
    
    # 等待 10 秒後再檢查
    sleep 10
done
