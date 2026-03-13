# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SGP4py - SGP4 轨道外推算法的 Python 实现，用于 TLE 数据轨道计算。基于 Spacetrack Report #3 和 CelesTrak 参考实现。

## Usage

运行主程序（默认使用 ISS.tle）:
```
python sgp4.py
```

指定 TLE 文件:
```
python sgp4.py <tle_file>
```

## Code Structure

- `sgp4.py` - 核心实现，包含:
  - `TLE` - TLE 数据解析和存储
  - `SGP4Propagator` - SGP4 轨道外推器
  - `PositionVelocity` - 位置和速度矢量 (ECI 坐标系)
  - 辅助函数：`parse_tle`, `load_tle_file`, `propagate_to_datetime`, `get_lat_lon_alt`

- `*.tle` - TLE 数据文件 (ISS.tle, NOAA19.tle, shenzhou15.tle)

## Key Constants (WGS72)

- `R_EARTH = 6378.135 km`
- `MU_EARTH = 398600.8 km^3/s^2`
- `J2 = 1.082616e-3`
