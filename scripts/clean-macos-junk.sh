#!/bin/zsh
# =============================================================================
# clean-macos-junk.sh
# 清除 macOS 產生的 ._* 和 .DS_Store 垃圾檔案
# 用法：
#   ./scripts/clean-macos-junk.sh           # 單次清理
#   ./scripts/clean-macos-junk.sh --watch   # 背景監控模式（持續監聽並刪除）
# =============================================================================

set -euo pipefail

# 取得腳本所在目錄的父目錄（專案根目錄）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# 顏色輸出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# 單次清理函數
clean_once() {
    log_info "開始清理 macOS 垃圾檔案..."
    log_info "專案目錄: $PROJECT_ROOT"
    
    local count=0
    
    # 清理工作目錄中的 ._* 檔案
    while IFS= read -r -d '' file; do
        rm -f "$file"
        log_warn "已刪除: $file"
        ((count++))
    done < <(find "$PROJECT_ROOT" -name '._*' -type f -print0 2>/dev/null)
    
    # 清理 .DS_Store
    while IFS= read -r -d '' file; do
        rm -f "$file"
        log_warn "已刪除: $file"
        ((count++))
    done < <(find "$PROJECT_ROOT" -name '.DS_Store' -type f -print0 2>/dev/null)
    
    # 清理 .git/objects/pack 中的 ._* 檔案（這些會導致 git 錯誤）
    local git_pack="$PROJECT_ROOT/.git/objects/pack"
    if [[ -d "$git_pack" ]]; then
        while IFS= read -r -d '' file; do
            rm -f "$file"
            log_warn "已刪除 (git pack): $file"
            ((count++))
        done < <(find "$git_pack" -name '._*' -type f -print0 2>/dev/null)
    fi
    
    if [[ $count -eq 0 ]]; then
        log_info "沒有找到需要清理的檔案"
    else
        log_info "已清理 $count 個檔案"
    fi
}

# 背景監控模式（使用 fswatch 或 polling）
watch_mode() {
    log_info "進入背景監控模式..."
    log_info "專案目錄: $PROJECT_ROOT"
    log_info "按 Ctrl+C 停止監控"
    
    # 檢查是否有 fswatch
    if command -v fswatch &>/dev/null; then
        log_info "使用 fswatch 進行即時監控..."
        fswatch -0 --event Created --event Renamed "$PROJECT_ROOT" | while IFS= read -r -d '' event; do
            # 檢查是否為 ._* 或 .DS_Store
            if [[ "$event" == *"/._"* ]] || [[ "$event" == *"/.DS_Store" ]]; then
                if [[ -f "$event" ]]; then
                    rm -f "$event"
                    log_warn "即時刪除: $event"
                fi
            fi
        done
    else
        log_warn "未安裝 fswatch，使用輪詢模式（每 5 秒檢查一次）"
        log_info "建議安裝 fswatch 以獲得更好的效能: brew install fswatch"
        
        while true; do
            # 清理 ._* 檔案
            find "$PROJECT_ROOT" -name '._*' -type f -print0 2>/dev/null | while IFS= read -r -d '' file; do
                rm -f "$file"
                log_warn "已刪除: $file"
            done
            
            # 清理 .DS_Store
            find "$PROJECT_ROOT" -name '.DS_Store' -type f -print0 2>/dev/null | while IFS= read -r -d '' file; do
                rm -f "$file"
                log_warn "已刪除: $file"
            done
            
            # 清理 .git/objects/pack 中的垃圾
            local git_pack="$PROJECT_ROOT/.git/objects/pack"
            if [[ -d "$git_pack" ]]; then
                find "$git_pack" -name '._*' -type f -print0 2>/dev/null | while IFS= read -r -d '' file; do
                    rm -f "$file"
                    log_warn "已刪除 (git pack): $file"
                done
            fi
            
            sleep 5
        done
    fi
}

# 主程式
main() {
    case "${1:-}" in
        --watch|-w)
            watch_mode
            ;;
        --help|-h)
            echo "用法: $0 [選項]"
            echo ""
            echo "選項:"
            echo "  (無)        單次清理所有 ._* 和 .DS_Store 檔案"
            echo "  --watch, -w 背景監控模式，持續監聽並自動刪除"
            echo "  --help, -h  顯示此說明"
            ;;
        *)
            clean_once
            ;;
    esac
}

main "$@"
