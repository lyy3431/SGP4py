# SGP4py

SGP4 (Simplified General Perturbations 4) 轨道外推算法的 Python 实现，用于根据 TLE 数据计算卫星轨道位置和速度。

## 功能特性

- 解析 TLE (Two-Line Element) 轨道根数文件
- SGP4 轨道外推算法
- ECI (Earth-Centered Inertial) 坐标系位置和速度计算
- 星下点 (经纬度) 和轨道高度计算
- 未来轨道预报

## 使用方法

运行主程序（默认使用 ISS.tle）:

```bash
python sgp4.py
```

指定 TLE 文件:

```bash
python sgp4.py <tle_file>
```

## 输出示例

![程序输出示例](screenshot.png)

输出信息包括:
- TLE 历元时刻的轨道状态
- 当前时刻的卫星位置
- 未来 24 小时轨道预报（每 15 分钟）

## TLE 数据格式

TLE (Two-Line Element) 是 NASA/NORAD 提供的卫星轨道根数数据格式，包含:
- 卫星编号和名称
- 历元时间
- 轨道倾角、升交点赤经、偏心率
- 近地点幅角、平近点角
- 平均运动、阻力系数等

项目内置了以下 TLE 文件:
- `ISS.tle` - 国际空间站
- `NOAA19.tle` - NOAA-19 气象卫星
- `shenzhou15.tle` - 神舟十五号

## 技术实现

基于 Spacetrack Report #3 和 CelesTrak 参考实现，使用 WGS72 坐标系常数:

- 地球半径：6378.135 km
- 地球重力常数：398600.8 km³/s²
- J2 带谐系数：1.082616e-3

## 依赖

仅使用 Python 标准库，无需额外安装包。
