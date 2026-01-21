#!/bin/bash
# Git 自動同步控制腳本
# 用法:
#   ./scripts/git-sync.sh start [間隔秒數]  - 啟動自動同步 (預設300秒)
#   ./scripts/git-sync.sh stop              - 停止自動同步
#   ./scripts/git-sync.sh status            - 查看狀態

PROJECT_DIR="/Volumes/Seagate/D2Q9_PeriodicHill/D2Q9_PeriodicHill"
PID_FILE="${PROJECT_DIR}/scripts/.sync.pid"
LOG_FILE="${PROJECT_DIR}/auto-sync.log"

cd "$PROJECT_DIR" || exit 1

start_sync() {
    INTERVAL=${1:-300}

    # 檢查是否已在運行
    if [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        echo "⚠️  自動同步已在運行中 (PID: $(cat "$PID_FILE"))"
        return 1
    fi

    # 背景啟動同步
    (
        echo "=========================================="
        echo "  Git 自動同步已啟動"
        echo "  同步間隔: ${INTERVAL} 秒"
        echo "=========================================="

        while true; do
            if [[ -n $(git status --porcelain) ]]; then
                TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
                echo ""
                echo "[${TIMESTAMP}] 偵測到修改，開始同步..."
                git add .
                git commit -m "auto-sync: ${TIMESTAMP}"
                git push 2>&1
                echo "[${TIMESTAMP}] 同步完成!"
            fi
            sleep "$INTERVAL"
        done
    ) >> "$LOG_FILE" 2>&1 &

    echo $! > "$PID_FILE"
    echo "✅ 自動同步已啟動 (PID: $!, 間隔: ${INTERVAL}秒)"
    echo "   日誌: $LOG_FILE"
}

stop_sync() {
    if [[ -f "$PID_FILE" ]]; then
        PID=$(cat "$PID_FILE")
        if kill -0 "$PID" 2>/dev/null; then
            kill "$PID"
            rm -f "$PID_FILE"
            echo "✅ 自動同步已停止 (PID: $PID)"
        else
            rm -f "$PID_FILE"
            echo "⚠️  進程已不存在，已清理 PID 檔案"
        fi
    else
        echo "ℹ️  自動同步未在運行"
    fi
}

status_sync() {
    if [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        echo "✅ 自動同步運行中 (PID: $(cat "$PID_FILE"))"
        echo ""
        echo "最近日誌:"
        tail -5 "$LOG_FILE" 2>/dev/null
    else
        echo "❌ 自動同步未運行"
    fi
}

case "${1:-status}" in
    start)
        start_sync "$2"
        ;;
    stop)
        stop_sync
        ;;
    status)
        status_sync
        ;;
    *)
        echo "用法: $0 {start|stop|status} [間隔秒數]"
        echo ""
        echo "  start [秒數]  啟動自動同步 (預設 300 秒)"
        echo "  stop          停止自動同步"
        echo "  status        查看運行狀態"
        ;;
esac
