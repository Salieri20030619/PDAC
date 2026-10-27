#!/bin/bash
# 检查 Task A Figure L2 相关图是否已在服务器生成
# 用法: 在 sgs10:/mnt/BioDisk/wangmingyuan/ 下直接执行此脚本

cd /mnt/BioDisk/wangmingyuan/

echo "========================================"
echo "1. 列出工作目录所有候选 figure 目录"
echo "========================================"
find . -maxdepth 3 -type d \( -iname "*figure*" -o -iname "*fig*" -o -iname "*plot*" -o -iname "*L2*" -o -iname "*TaskA*" -o -iname "*taska*" \) 2>/dev/null

echo ""
echo "========================================"
echo "2. 关键词搜索 L2 相关图文件 (png/pdf/jpg/svg/tiff)"
echo "========================================"
echo ""
echo "--- (a) 文件名含 L2 ---"
find . -type f \( -iname "*L2*.png" -o -iname "*L2*.pdf" -o -iname "*L2*.jpg" -o -iname "*L2*.svg" -o -iname "*L2*.tiff" \) 2>/dev/null

echo ""
echo "--- (b) 文件名含 TGFb / tgfb ---"
find . -type f \( -iname "*tgfb*.png" -o -iname "*tgfb*.pdf" -o -iname "*tgfb*.jpg" -o -iname "*tgfb*.svg" \) 2>/dev/null

echo ""
echo "--- (c) 文件名含 Smad3 / smad3 ---"
find . -type f \( -iname "*smad3*.png" -o -iname "*smad3*.pdf" -o -iname "*smad3*.jpg" -o -iname "*smad3*.svg" \) 2>/dev/null

echo ""
echo "--- (d) 文件名含 host_zone 或 host-zone ---"
find . -type f \( -iname "*host_zone*" -o -iname "*host-zone*" -o -iname "*hostzone*" \) 2>/dev/null | grep -iE "\.(png|pdf|jpg|svg|tiff)$"

echo ""
echo "--- (e) 文件名含 pct heatmap / pathway heatmap ---"
find . -type f \( -iname "*pct*heatmap*" -o -iname "*pathway*heatmap*" -o -iname "*heatmap*tgfb*" \) 2>/dev/null | grep -iE "\.(png|pdf|jpg|svg|tiff)$"

echo ""
echo "--- (f) 文件名含 paracrine / autocrine / pearson ---"
find . -type f \( -iname "*paracrine*" -o -iname "*autocrine*" -o -iname "*pearson*" \) 2>/dev/null | grep -iE "\.(png|pdf|jpg|svg|tiff)$"

echo ""
echo "========================================"
echo "3. 列出 Task A / TaskA 相关目录最近修改的 figure"
echo "========================================"
for d in $(find . -maxdepth 3 -type d \( -iname "*TaskA*" -o -iname "*taska*" -o -iname "*task_a*" \) 2>/dev/null); do
  echo "--- 目录: $d ---"
  find "$d" -maxdepth 4 -type f \( -iname "*.png" -o -iname "*.pdf" -o -iname "*.jpg" -o -iname "*.svg" \) 2>/dev/null | head -40
done

echo ""
echo "========================================"
echo "4. 兜底: 最近 30 天修改且文件名含 figure/fig 的图"
echo "========================================"
find . -type f \( -iname "*figure*.png" -o -iname "*figure*.pdf" -o -iname "*fig*.png" -o -iname "*fig*.pdf" \) -mtime -30 2>/dev/null | head -50

echo ""
echo "========================================"
echo "完成. 如果以上几项都没结果, 说明 Figure L2 还没生成."
echo "========================================"
