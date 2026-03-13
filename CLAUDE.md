# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SGP4py - SGP4/HPOP 轨道外推算法的 Python 实现，用于 TLE 数据轨道计算。

## Usage

**SGP4 模式（快速）:**
```bash
python sgp4.py              # 默认 ISS.tle
python sgp4.py <tle_file>   # 指定 TLE 文件
```

**HPOP 模式（高精度）:**
```bash
python sgp4.py --hpop       # 使用 HPOP 高精度外推
python hpop.py              # 直接运行 HPOP 模块
```

## Code Structure

- `sgp4.py` - SGP4 核心实现
  - `TLE` - TLE 数据解析和存储
  - `SGP4Propagator` - SGP4 轨道外推器
  - `PositionVelocity` - 位置和速度矢量 (ECI 坐标系)
  - `propagate_to_datetime` - 外推到指定时间

- `hpop.py` - HPOP 高精度外推模块
  - `HPOPPropagator` - HPOP 外推器（使用 RK45 数值积分）
  - `SpacecraftState` - 航天器物理参数（质量、面积）
  - `hprop_from_tle` - 从 TLE 开始 HPOP 外推
  - `compare_sgp4_hpop` - 对比 SGP4 和 HPOP 结果

- `*.tle` - TLE 数据文件 (ISS.tle, NOAA19.tle, shenzhou15.tle)

## SGP4 vs HPOP

| 模式 | 精度 | 速度 | 摄动模型 |
|------|------|------|----------|
| SGP4 | ~1km | 快 | J2-J4 + 简化阻力 |
| HPOP | ~10m | 慢 | J2-J10 + 日月引力 + 大气阻力 + 太阳辐射压 + 相对论 |

## Key Constants (WGS72)

- `R_EARTH = 6378.135 km`
- `MU_EARTH = 398600.8 km^3/s^2`
- `J2 = 1.082616e-3`
