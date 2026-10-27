#!/usr/bin/env Rscript
# 测试装 ferroviz (含 FerrDb v2 数据)

# 先看是否有 devtools / remotes
have_remotes <- requireNamespace("remotes", quietly = TRUE)
have_devtools <- requireNamespace("devtools", quietly = TRUE)
cat(sprintf("remotes: %s\ndevtools: %s\n", have_remotes, have_devtools))

# 若没有, 先装一个轻量级的 remotes
if (!have_remotes && !have_devtools) {
    cat("装 remotes (轻量级, 比 devtools 简单)...\n")
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

# 装 ferroviz
cat("\n开始装 ferroviz ...\n")
if (have_devtools) {
    devtools::install_github("cdesterke/ferroviz", upgrade = "never")
} else {
    remotes::install_github("cdesterke/ferroviz", upgrade = "never")
}

# 测试
cat("\n测试加载...\n")
library(ferroviz)
data(ferrdbhs)
cat(sprintf("ferrdbhs 维度: %d × %d\n", nrow(ferrdbhs), ncol(ferrdbhs)))
cat("ferrdbhs 列名:", paste(colnames(ferrdbhs), collapse=", "), "\n")
cat("前 3 行:\n")
print(head(ferrdbhs, 3))
