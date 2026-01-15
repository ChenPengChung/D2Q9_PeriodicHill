# Git Codespace 工作流程指南

---

## 本地 vs Codespace 環境說明

```
┌─────────────────────────┐
│  你的 Mac（本地電腦）     │
│  ✅ 可啟用自動同步腳本    │  ← 主要開發環境
│  編輯 → 自動 push        │
└─────────────────────────┘
            ↕ 同步
┌─────────────────────────┐
│  GitHub 雲端              │  ← 備份中心
└─────────────────────────┘
            ↕ 獨立
┌─────────────────────────┐
│  Codespace（雲端編輯器）  │
│  ❌ 需手動 commit + push │  ← 外出/臨時編輯
└─────────────────────────┘
```

### 重要觀念
- **本地電腦** 和 **Codespace** 是兩個獨立環境
- 自動同步腳本只在啟動的環境中運作
- 切換環境時，記得先同步！

---

## 不同地點的工作流程

| 地點 | 做法 |
|------|------|
| **在家（本地）** | 啟動自動同步腳本，只管編輯 |
| **在外（Codespace）** | 編輯完手動 `git add . && git commit -m "..." && git push` |
| **回家後** | 先執行 `git pull` 把 Codespace 的修改拉下來 |

---

## 本地自動同步設定

### 啟動自動同步（每 5 分鐘）
```bash
cd /Volumes/Seagate/D2Q9_PeriodicHill
nohup ./scripts/auto-sync.sh 300 > auto-sync.log 2>&1 &
```

### 查看同步日誌
```bash
cat auto-sync.log
```

### 停止自動同步
```bash
pkill -f auto-sync.sh
```

### 確認腳本是否運行中
```bash
ps aux | grep auto-sync | grep -v grep
```

---

## 日常工作流程（每次編輯完程式碼後）

### 1. 查看修改狀態
```bash
git status
```
顯示哪些檔案被修改、新增或刪除。

### 2. 查看具體修改內容
```bash
git diff
```
查看所有未暫存的修改細節。

### 3. 暫存修改的檔案

暫存所有修改：
```bash
git add .
```

或暫存特定檔案：
```bash
git add 檔案名稱.h
```

### 4. 確認暫存內容
```bash
git status
```
確認要提交的檔案（顯示為綠色）。

### 5. 提交修改
```bash
git commit -m "簡短描述這次修改的內容"
```

範例：
```bash
git commit -m "修復 MRT 碰撞運算邏輯"
git commit -m "新增壁面邊界條件處理"
git commit -m "優化 streaming 步驟效能"
```

### 6. 推送到 GitHub
```bash
git push
```

---

## 快速一行指令（適合小修改）

```bash
git add . && git commit -m "修改描述" && git push
```

---

## 完整工作流程範例

```bash
# 1. 編輯完程式碼後，先查看狀態
git status

# 2. 查看修改了什麼
git diff

# 3. 暫存所有修改
git add .

# 4. 提交並寫說明
git commit -m "實作 main.cu 初始化流場功能"

# 5. 推送到 GitHub
git push
```

---

## 從 GitHub 同步最新程式碼（多裝置協作時）

```bash
git pull
```

---

## 其他常用指令

### 查看提交歷史
```bash
git log --oneline -10
```

### 取消未暫存的修改（還原檔案）
```bash
git checkout -- 檔案名稱.h
```

### 取消已暫存的檔案（從 staging 移除）
```bash
git reset HEAD 檔案名稱.h
```

### 修改上一次的 commit 訊息
```bash
git commit --amend -m "新的訊息"
```

---

## Codespace 專用

### 在 Codespace 中開啟終端機
- 快捷鍵：`Ctrl + `` ` 或 `Cmd + `` `
- 或點擊選單：Terminal → New Terminal

### Codespace 會自動保存，但仍需手動 commit & push！

---

## 建議的 Commit 訊息格式

```
類型: 簡短描述

類型包括：
- feat: 新功能
- fix: 修復 bug
- refactor: 重構程式碼
- docs: 文件更新
- test: 測試相關
```

範例：
```bash
git commit -m "feat: 新增 MRT 碰撞運算完整邏輯"
git commit -m "fix: 修復邊界條件判斷錯誤"
git commit -m "refactor: 重構插值函數提升可讀性"
```
