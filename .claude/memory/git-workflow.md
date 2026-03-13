# Git Workflow

## 推送规则

- 默认直接推送
- 如果推送失败，自动尝试使用代理推送
- 代理地址：`http://localhost:7897`

## 代理配置命令

```bash
git config http.proxy http://localhost:7897
git config https.proxy http://localhost:7897
```
