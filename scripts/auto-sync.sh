#!/bin/bash
# 自動同步腳本 - 每隔指定時間自動 commit 並 push
# 用法: ./scripts/auto-sync.sh [間隔秒數]
# 預設間隔: 300 秒 (5分鐘)

INTERVAL=${1:-300}

echo "=========================================="
echo "  Git 自動同步已啟動"
echo "  同步間隔: ${INTERVAL} 秒"
echo "  按 Ctrl+C 停止"
echo "=========================================="

while true; do
    # 檢查是否有修改
    if [[ -n $(git status --porcelain) ]]; then
        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

        echo ""
        echo "[${TIMESTAMP}] 偵測到修改，開始同步..."

        # 暫存所有修改
        git add .

        # 提交
        git commit -m "auto-sync: ${TIMESTAMP}"

        # 推送
        git push

        echo "[${TIMESTAMP}] 同步完成!"
    else
        echo "[$(date '+%H:%M:%S')] 無修改，等待中..."
    fi

    sleep $INTERVAL
done
