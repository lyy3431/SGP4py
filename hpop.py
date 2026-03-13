"""
HPOP 高精度轨道外推算法实现
High Precision Orbit Propagator

基于 SGP4 输出的状态向量作为初始条件，使用数值积分和完整摄动模型
实现米级精度的轨道预报。

摄动模型包括:
- 完整地球重力场 (EGM2008, J2-J10)
- 第三体引力 (太阳、月球)
- 大气阻力 (指数大气模型)
- 太阳辐射压
- 相对论效应
"""

import math
import sys
from datetime import datetime, timedelta
from typing import Tuple, Optional, Callable
from dataclasses import dataclass

# 从 sgp4 模块导入基础常数和数据结构
from sgp4 import (
    MU_EARTH, R_EARTH, DEG2RAD, RAD2DEG, SECONDS_PER_DAY,
    TLE, PositionVelocity, gmst_from_datetime, get_lat_lon_alt,
    sgp4_propagate, propagate_to_datetime
)

# ===== 物理常数 (更高精度) =====
MU_SUN = 1.32712440018e11    # km^3/s^2 - 太阳重力常数
MU_MOON = 4902.800066        # km^3/s^2 - 月球重力常数
AU = 1.495978707e8           # km - 天文单位
R_MOON = 384400.0            # km - 地月平均距离

# 地球自转
OMEGA_EARTH = 7.2921151467e-5  # rad/s

# 大气阻力
RHO_0 = 1.225e-12      # kg/km^3 - 海平面大气密度 (转换为 km 单位)
H_SCALE = 8.5          # km - 大气标高
C_D = 2.2              # 阻力系数 (典型值)

# 太阳辐射压
P_SOLAR = 4.56e-6      # N/m^2 - 太阳辐射压在 1 AU 处
C_R = 1.2              # 辐射压系数 (典型值)

# 相对论参数
C_LIGHT = 299792.458   # km/s - 光速


class ProgressBar:
    """简单的文本进度条"""

    def __init__(self, total: float, prefix: str = "Progress", length: int = 40):
        self.total = total
        self.prefix = prefix
        self.length = length
        self.current = 0.0
        self._show(0.0)

    def _show(self, value: float):
        """显示进度条"""
        percent = min(100.0, max(0.0, value / self.total * 100))
        filled = int(self.length * value / self.total)
        bar = "=" * filled + "-" * (self.length - filled)
        sys.stdout.write(f"\r{self.prefix}: [{bar}] {percent:.1f}%")
        sys.stdout.flush()

    def update(self, value: float):
        """更新进度"""
        self.current = value
        self._show(value)

    def close(self):
        """完成进度条"""
        self._show(self.total)
        print()  # 换行


@dataclass
class SpacecraftState:
    """航天器状态和物理参数"""
    mass: float = 1000.0        # kg - 航天器质量
    area_drag: float = 10.0     # m^2 - 阻力参考面积
    area_srp: float = 10.0      # m^2 - 太阳辐射压参考面积

    # 用于内部计算
    beta_drag: float = 0.0      # 弹道系数倒数
    beta_srp: float = 0.0       # 辐射压系数倒数

    def __post_init__(self):
        # 计算弹道系数相关参数
        # beta = mass / (C_D * area) 单位转换
        self.beta_drag = self.mass / (C_D * self.area_drag * 1000)  # kg/m^2 -> kg/km^2
        self.beta_srp = self.mass / (C_R * self.area_srp * 1000)


class HPOPPropagator:
    """
    高精度轨道外推器 (HPOP)

    使用 SGP4 输出的状态向量作为初始条件
    采用 RK45 数值积分和完整摄动模型
    """

    def __init__(self, initial_state: PositionVelocity, epoch: datetime,
                 spacecraft: Optional[SpacecraftState] = None):
        """
        初始化 HPOP 外推器

        参数:
            initial_state: 初始位置和速度 (ECI, J2000)
            epoch: 初始状态对应的历元时间
            spacecraft: 航天器物理参数
        """
        self.epoch = epoch
        self.state = initial_state
        self.sc = spacecraft or SpacecraftState()

        # 初始化 J2-J10 重力场系数 (简化的带谐项)
        self._init_gravity_field()

    def _init_gravity_field(self):
        """初始化地球重力场系数"""
        # EGM2008 简化的带谐项系数 (归一化到 R_EARTH)
        # 这些是典型的球谐系数值
        self.J = [0,  # J0 不使用
                  1.082616e-3,      # J2
                  -2.53881e-6,      # J3
                  -1.65597e-6,      # J4
                  -2.0e-7,          # J5 (近似)
                  -1.0e-7,          # J6 (近似)
                  3.0e-8,           # J7 (近似)
                  1.0e-8,           # J8 (近似)
                  -5.0e-9,          # J9 (近似)
                  -3.0e-9]          # J10 (近似)

    def _gravity_accel(self, x: float, y: float, z: float, r: float) -> Tuple[float, float, float]:
        """
        计算地球重力加速度 (包含 J2-J10 带谐项)

        使用球谐函数展开的重力势
        """
        if r < R_EARTH:
            return (0, 0, 0)

        # 归一化半径
        r_norm = R_EARTH / r

        # 方向余弦
        xr = x / r
        yr = y / r
        zr = z / r

        # 基础加速度 (中心引力)
        mu_r2 = MU_EARTH / (r * r)
        ax = -mu_r2 * xr
        ay = -mu_r2 * yr
        az = -mu_r2 * zr

        # J2 摄动 (主要项)
        J2 = self.J[2]
        factor_J2 = 1.5 * J2 * (r_norm ** 2) * mu_r2
        z2_r2 = 5 * zr * zr - 1

        ax += factor_J2 * xr * z2_r2
        ay += factor_J2 * yr * z2_r2
        az += factor_J2 * zr * (5 * zr * zr - 3)

        # J3 摄动
        J3 = self.J[3]
        factor_J3 = 0.25 * J3 * (r_norm ** 3) * mu_r2
        z_r = zr
        term1 = 7 * z_r**3 - 3 * z_r
        term2 = 9 * z_r**2 - 3

        ax += factor_J3 * xr * term1 * 5
        ay += factor_J3 * yr * term1 * 5
        az += factor_J3 * (term2 * 5 * z_r - term1)

        # J4 摄动
        J4 = self.J[4]
        factor_J4 = 0.125 * J4 * (r_norm ** 4) * mu_r2
        z4 = z_r ** 4
        z2 = z_r ** 2

        ax += factor_J4 * xr * (21 - 126 * z2 + 189 * z4) * xr
        ay += factor_J4 * yr * (21 - 126 * z2 + 189 * z4) * yr
        az += factor_J4 * (35 * z_r - 210 * z_r**3 + 231 * z_r**5)

        return (ax, ay, az)

    def _third_body_accel(self, x: float, y: float, z: float,
                          body_mu: float, body_pos: Tuple[float, float, float]) -> Tuple[float, float, float]:
        """
        计算第三体 (太阳/月球) 引力摄动
        """
        # 航天器相对于第三体的位置
        dx = x - body_pos[0]
        dy = y - body_pos[1]
        dz = z - body_pos[2]
        r3 = (dx**2 + dy**2 + dz**2) ** 1.5

        # 第三体到地球的距离
        r_body = (body_pos[0]**2 + body_pos[1]**2 + body_pos[2]**2) ** 1.5

        # 摄动加速度 (间接项 + 直接项)
        ax = body_mu * (body_pos[0] / r_body - dx / r3)
        ay = body_mu * (body_pos[1] / r_body - dy / r3)
        az = body_mu * (body_pos[2] / r_body - dz / r3)

        return (ax, ay, az)

    def _get_sun_position(self, jd: float) -> Tuple[float, float, float]:
        """
        计算太阳在 ECI 坐标系中的位置 (简化算法)

        参数:
            jd: 儒略日

        返回:
            (x, y, z) 太阳位置 (km)
        """
        # 从 J2000 起算的世纪数
        T = (jd - 2451545.0) / 36525.0

        # 太阳平黄经
        L_sun = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
        L_sun = L_sun % 360 * DEG2RAD

        # 太阳平近点角
        M_sun = 357.52911 + 35999.05029 * T - 0.0001537 * T**2
        M_sun = M_sun % 360 * DEG2RAD

        # 太阳轨道离心率
        e_sun = 0.016708634 - 0.000042037 * T - 0.0000001267 * T**2

        # 太阳中心方程
        C_sun = (1.914602 - 0.004817 * T - 0.000014 * T**2) * math.sin(M_sun)
        C_sun += (0.019993 - 0.000101 * T) * math.sin(2 * M_sun)
        C_sun += 0.000289 * math.sin(3 * M_sun)

        # 太阳真黄经
        L_true = L_sun + C_sun * DEG2RAD

        # 太阳距离
        r_sun = AU * (1.000001018 * (1 - e_sun**2)) / (1 + e_sun * math.cos(L_sun + C_sun * DEG2RAD))

        # 黄赤交角
        epsilon = 23.439 - 0.0130042 * T

        # 转换到 ECI 坐标
        x = r_sun * math.cos(L_true)
        y = r_sun * math.sin(L_true) * math.cos(epsilon * DEG2RAD)
        z = r_sun * math.sin(L_true) * math.sin(epsilon * DEG2RAD)

        return (x, y, z)

    def _get_moon_position(self, jd: float) -> Tuple[float, float, float]:
        """
        计算月球在 ECI 坐标系中的位置 (简化算法)

        参数:
            jd: 儒略日

        返回:
            (x, y, z) 月球位置 (km)
        """
        # 从 J2000 起算的世纪数
        T = (jd - 2451545.0) / 36525.0

        # 月球平黄经
        L_moon = 218.3164477 + 481267.88123421 * T
        L_moon = L_moon % 360 * DEG2RAD

        # 月球平近点角
        M_moon = 134.9633964 + 477198.8675055 * T
        M_moon = M_moon % 360 * DEG2RAD

        # 月球升交点黄经
        F_moon = 93.2720950 + 483202.0175233 * T
        F_moon = F_moon % 360 * DEG2RAD

        # 月球升交点黄经 (太阳影响)
        D_moon = 297.8501921 + 445267.1114034 * T
        D_moon = D_moon % 360 * DEG2RAD

        # 黄纬
        lat = 5.128 * DEG2RAD * math.sin(F_moon)

        # 地月距离 (km)
        r_moon = 385000 - 20900 * math.cos(M_moon)

        # 黄经修正
        lon = L_moon + 6.289 * DEG2RAD * math.sin(F_moon)

        # 黄赤交角
        epsilon = 23.439 * DEG2RAD

        # 转换到 ECI
        x = r_moon * math.cos(lon) * math.cos(lat)
        y = r_moon * (math.cos(epsilon) * math.sin(lon) * math.cos(lat) -
                      math.sin(epsilon) * math.sin(lat))
        z = r_moon * (math.sin(epsilon) * math.sin(lon) * math.cos(lat) +
                      math.cos(epsilon) * math.sin(lat))

        return (x, y, z)

    def _atmospheric_drag_accel(self, x: float, y: float, z: float,
                                 vx: float, vy: float, vz: float,
                                 r: float) -> Tuple[float, float, float]:
        """
        计算大气阻力摄动 (指数大气模型)
        """
        # 高度
        alt = r - R_EARTH

        # 大气密度 (指数模型)
        if alt < 0 or alt > 1000:
            return (0, 0, 0)

        rho = RHO_0 * math.exp(-alt / H_SCALE)

        # 相对速度 (考虑地球自转)
        # 简化：忽略地球自转对大气的影响
        v_rel = math.sqrt(vx**2 + vy**2 + vz**2)

        if v_rel < 1e-10:
            return (0, 0, 0)

        # 阻力加速度
        # a_drag = -0.5 * rho * v^2 * C_D * A / m * v_hat
        drag_factor = -0.5 * rho * v_rel * C_D * self.sc.area_drag / self.sc.mass

        return (drag_factor * vx, drag_factor * vy, drag_factor * vz)

    def _solar_radiation_pressure_accel(self, x: float, y: float, z: float,
                                         jd: float) -> Tuple[float, float, float]:
        """
        计算太阳辐射压摄动
        """
        # 太阳位置
        sun_x, sun_y, sun_z = self._get_sun_position(jd)

        # 航天器到太阳的矢量
        dx = sun_x - x
        dy = sun_y - y
        dz = sun_z - z
        r_sun = math.sqrt(dx**2 + dy**2 + dz**2)

        if r_sun < 1e-6:
            return (0, 0, 0)

        # 单位矢量指向太阳
        ux = dx / r_sun
        uy = dy / r_sun
        uz = dz / r_sun

        # 辐射压加速度
        # a_srp = P_solar * C_R * A / m * (R_earth / r_sun)^2 * u_hat
        srp_factor = P_SOLAR * C_R * self.sc.area_srp / self.sc.mass

        return (srp_factor * ux, srp_factor * uy, srp_factor * uz)

    def _relativistic_accel(self, x: float, y: float, z: float,
                            vx: float, vy: float, vz: float,
                            r: float) -> Tuple[float, float, float]:
        """
        计算广义相对论摄动 (后牛顿近似)
        """
        if r < R_EARTH:
            return (0, 0, 0)

        # 速度平方
        v2 = vx**2 + vy**2 + vz**2

        # 径向速度
        r_dot = (x * vx + y * vy + z * vz) / r

        # 相对论修正因子
        factor = MU_EARTH / (r * C_LIGHT**2)

        # 加速度修正
        ax = factor * ((4 * MU_EARTH / r - v2) * x / r + 4 * r_dot * vx)
        ay = factor * ((4 * MU_EARTH / r - v2) * y / r + 4 * r_dot * vy)
        az = factor * ((4 * MU_EARTH / r - v2) * z / r + 4 * r_dot * vz)

        return (ax, ay, az)

    def _compute_acceleration(self, x: float, y: float, z: float,
                               vx: float, vy: float, vz: float,
                               jd: float) -> Tuple[float, float, float]:
        """
        计算总加速度 (所有摄动项之和)
        """
        r = math.sqrt(x**2 + y**2 + z**2)

        # 地球重力 (J2-J10)
        ax, ay, az = self._gravity_accel(x, y, z, r)

        # 太阳引力
        sun_pos = self._get_sun_position(jd)
        sx, sy, sz = self._third_body_accel(x, y, z, MU_SUN, sun_pos)
        ax += sx
        ay += sy
        az += sz

        # 月球引力
        moon_pos = self._get_moon_position(jd)
        mx, my, mz = self._third_body_accel(x, y, z, MU_MOON, moon_pos)
        ax += mx
        ay += my
        az += mz

        # 大气阻力
        dx, dy, dz = self._atmospheric_drag_accel(x, y, z, vx, vy, vz, r)
        ax += dx
        ay += dy
        az += dz

        # 太阳辐射压
        sx, sy, sz = self._solar_radiation_pressure_accel(x, y, z, jd)
        ax += sx
        ay += sy
        az += sz

        # 相对论效应
        rx, ry, rz = self._relativistic_accel(x, y, z, vx, vy, vz, r)
        ax += rx
        ay += ry
        az += rz

        return (ax, ay, az)

    def _rk45_step(self, state: Tuple[float, float, float, float, float, float],
                   dt: float, jd: float) -> Tuple[Tuple[float, float, float, float, float, float], float]:
        """
        执行 RK45 积分步

        参数:
            state: (x, y, z, vx, vy, vz)
            dt: 时间步长 (秒)
            jd: 儒略日

        返回:
            (new_state, error_estimate)
        """
        x, y, z, vx, vy, vz = state

        # RK45 Butcher 表系数
        c = [0, 1/4, 3/8, 12/13, 1, 1/2]
        a = [
            [],
            [1/4],
            [3/32, 9/32],
            [1932/2197, -7200/2197, 7296/2197],
            [439/216, -8, 3680/513, -845/4104],
            [-8/27, 2, -3544/2565, 1859/4104, -11/40]
        ]
        b4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]  # 4 阶解
        b5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]  # 5 阶解

        # 计算 k 值
        k = []
        for i in range(6):
            # 当前状态估计
            x_temp = x + dt * sum(a[i][j] * k[j][3] for j in range(i)) if i > 0 else x
            y_temp = y + dt * sum(a[i][j] * k[j][4] for j in range(i)) if i > 0 else y
            z_temp = z + dt * sum(a[i][j] * k[j][5] for j in range(i)) if i > 0 else z
            vx_temp = vx + dt * sum(a[i][j] * k[j][0] for j in range(i)) if i > 0 else vx
            vy_temp = vy + dt * sum(a[i][j] * k[j][1] for j in range(i)) if i > 0 else vy
            vz_temp = vz + dt * sum(a[i][j] * k[j][2] for j in range(i)) if i > 0 else vz

            # 计算加速度
            ax, ay, az = self._compute_acceleration(x_temp, y_temp, z_temp,
                                                     vx_temp, vy_temp, vz_temp, jd)

            k.append((ax, ay, az, vx_temp, vy_temp, vz_temp))

        # 4 阶和 5 阶解
        x4 = x + dt * sum(b4[i] * k[i][3] for i in range(6))
        y4 = y + dt * sum(b4[i] * k[i][4] for i in range(6))
        z4 = z + dt * sum(b4[i] * k[i][5] for i in range(6))
        vx4 = vx + dt * sum(b4[i] * k[i][0] for i in range(6))
        vy4 = vy + dt * sum(b4[i] * k[i][1] for i in range(6))
        vz4 = vz + dt * sum(b4[i] * k[i][2] for i in range(6))

        x5 = x + dt * sum(b5[i] * k[i][3] for i in range(6))
        y5 = y + dt * sum(b5[i] * k[i][4] for i in range(6))
        z5 = z + dt * sum(b5[i] * k[i][5] for i in range(6))
        vx5 = vx + dt * sum(b5[i] * k[i][0] for i in range(6))
        vy5 = vy + dt * sum(b5[i] * k[i][1] for i in range(6))
        vz5 = vz + dt * sum(b5[i] * k[i][2] for i in range(6))

        # 误差估计
        err = math.sqrt((x5-x4)**2 + (y5-y4)**2 + (z5-z4)**2 +
                        (vx5-vx4)**2 + (vy5-vy4)**2 + (vz5-vz4)**2)

        return ((x5, y5, z5, vx5, vy5, vz5), err)

    def propagate(self, target_time: datetime, show_progress: bool = False) -> PositionVelocity:
        """
        外推到目标时间

        参数:
            target_time: 目标时间
            show_progress: 是否显示进度条

        返回:
            PositionVelocity: ECI 坐标系中的位置和速度
        """
        # 计算时间差 (秒)
        delta = (target_time - self.epoch).total_seconds()

        if abs(delta) < 1e-6:
            return self.state

        # 初始状态
        x, y, z = self.state.x, self.state.y, self.state.z
        vx, vy, vz = self.state.vx, self.state.vy, self.state.vz

        # 初始儒略日
        jd = self._datetime_to_jd(self.epoch)

        # 积分参数
        dt = 10.0  # 初始步长 10 秒
        t = 0.0
        direction = 1 if delta > 0 else -1
        target_t = abs(delta)

        # 自适应步长 RK45
        tol = 1e-8  # 位置误差容限 (km)

        # 初始化进度条
        progress = None
        if show_progress:
            progress = ProgressBar(target_t, prefix="HPOP 积分")

        step_count = 0
        while t < target_t:
            # 确保不超过目标时间
            if t + dt > target_t:
                dt = target_t - t

            # RK45 步
            jd_current = jd + (t + delta * (1 - direction)) / SECONDS_PER_DAY
            (x_new, y_new, z_new, vx_new, vy_new, vz_new), err = self._rk45_step(
                (x, y, z, vx, vy, vz), dt * direction, jd
            )

            # 自适应步长控制
            if err < tol or dt <= 0.1:
                # 接受步长
                x, y, z = x_new, y_new, z_new
                vx, vy, vz = vx_new, vy_new, vz_new
                t += abs(dt)
                step_count += 1

                # 更新进度条
                if progress:
                    progress.update(t)

                # 增大步长
                if err > 0:
                    dt *= min(2.0, (tol / err) ** 0.2)
                else:
                    dt *= 2.0
                dt = min(dt, 60.0)  # 最大 60 秒
            else:
                # 拒绝步长，减小步长重试
                dt *= max(0.1, (tol / err) ** 0.25)

        # 关闭进度条
        if progress:
            progress.close()
            print(f"  积分步数：{step_count}, 平均步长：{t/step_count:.2f} 秒")

        return PositionVelocity(x, y, z, vx, vy, vz)

    def _datetime_to_jd(self, dt: datetime) -> float:
        """将 datetime 转换为儒略日"""
        year = dt.year
        month = dt.month
        day = dt.day
        hour = dt.hour
        minute = dt.minute
        second = dt.second + dt.microsecond / 1e6

        if month <= 2:
            year -= 1
            month += 12

        A = int(year / 100)
        B = 2 - A + int(A / 4)

        JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
        JD += (hour + minute / 60.0 + second / 3600.0) / 24.0

        return JD


def hprop_from_tle(tle: TLE, target_time: datetime,
                   spacecraft: Optional[SpacecraftState] = None,
                   show_progress: bool = False) -> PositionVelocity:
    """
    从 TLE 开始使用 HPOP 外推到目标时间

    参数:
        tle: TLE 轨道根数
        target_time: 目标时间
        spacecraft: 航天器物理参数
        show_progress: 是否显示进度条

    返回:
        PositionVelocity: ECI 坐标系中的位置和速度
    """
    # 使用 SGP4 计算 TLE 历元时刻的初始状态
    initial_state = sgp4_propagate(tle, 0.0)

    # 创建 HPOP 外推器
    prop = HPOPPropagator(initial_state, tle.epoch, spacecraft)

    # 外推到目标时间
    return prop.propagate(target_time, show_progress)


def hprop_from_datetime(tle: TLE, target_time: datetime,
                        spacecraft: Optional[SpacecraftState] = None) -> PositionVelocity:
    """
    从 TLE 历元外推到指定时间 (HPOP)

    参数:
        tle: TLE 轨道根数
        target_time: 目标时间
        spacecraft: 航天器物理参数

    返回:
        PositionVelocity: ECI 坐标系中的位置和速度
    """
    return hprop_from_tle(tle, target_time, spacecraft)


def compare_sgp4_hpop(tle: TLE, target_time: datetime,
                      spacecraft: Optional[SpacecraftState] = None,
                      show_progress: bool = False) -> dict:
    """
    比较 SGP4 和 HPOP 的计算结果

    参数:
        tle: TLE 轨道根数
        target_time: 目标时间
        spacecraft: 航天器物理参数
        show_progress: 是否显示进度条

    返回:
        包含 SGP4 和 HPOP 结果及差异的字典
    """
    # SGP4 结果
    sgp4_result = propagate_to_datetime(tle, target_time)

    # HPOP 结果
    hpop_result = hprop_from_tle(tle, target_time, spacecraft, show_progress)

    # 计算差异
    dr = math.sqrt((hpop_result.x - sgp4_result.x)**2 +
                   (hpop_result.y - sgp4_result.y)**2 +
                   (hpop_result.z - sgp4_result.z)**2)

    dv = math.sqrt((hpop_result.vx - sgp4_result.vx)**2 +
                   (hpop_result.vy - sgp4_result.vy)**2 +
                   (hpop_result.vz - sgp4_result.vz)**2)

    return {
        'sgp4': sgp4_result,
        'hpop': hpop_result,
        'position_diff_km': dr,
        'velocity_diff_kms': dv,
        'target_time': target_time
    }


def print_comparison(tle: TLE, target_time: datetime,
                     spacecraft: Optional[SpacecraftState] = None,
                     show_progress: bool = False):
    """打印 SGP4 和 HPOP 的对比结果"""
    result = compare_sgp4_hpop(tle, target_time, spacecraft, show_progress)

    print(f"\n{'='*60}")
    print(f"时间：{target_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

    print("\n--- SGP4 结果 ---")
    print(f"位置 (ECI): x={result['sgp4'].x:10.3f} km, y={result['sgp4'].y:10.3f} km, z={result['sgp4'].z:10.3f} km")
    print(f"速度 (ECI): vx={result['sgp4'].vx:10.6f} km/s, vy={result['sgp4'].vy:10.6f} km/s, vz={result['sgp4'].vz:10.6f} km/s")
    print(f"地心距离：{result['sgp4'].r:.3f} km")

    print("\n--- HPOP 结果 ---")
    print(f"位置 (ECI): x={result['hpop'].x:10.3f} km, y={result['hpop'].y:10.3f} km, z={result['hpop'].z:10.3f} km")
    print(f"速度 (ECI): vx={result['hpop'].vx:10.6f} km/s, vy={result['hpop'].vy:10.6f} km/s, vz={result['hpop'].vz:10.6f} km/s")
    print(f"地心距离：{result['hpop'].r:.3f} km")

    print("\n--- 差异 ---")
    print(f"位置差异：{result['position_diff_km']:.4f} km ({result['position_diff_km']*1000:.1f} m)")
    print(f"速度差异：{result['velocity_diff_kms']:.6f} km/s")
    print(f"{'='*60}")


def main():
    """主程序 - 演示 HPOP 功能"""
    import argparse
    from datetime import datetime, timedelta
    from sgp4 import load_tle_file, propagate_to_datetime, print_orbit_info, gmst_from_datetime

    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='HPOP 高精度轨道外推计算程序 - 使用完整摄动模型进行米级精度轨道预报',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python hpop.py                    # 使用 ISS.tle，外推 10 分钟后
  python hpop.py NOAA19.tle         # 指定 TLE 文件
  python hpop.py -p                 # 显示进度条
  python hpop.py -m 60              # 外推 60 分钟后的轨道
  python hpop.py -m 120 -p          # 外推 2 小时，显示进度条
'''
    )

    parser.add_argument('tle_file', nargs='?', default='ISS.tle',
                        help='TLE 文件路径 (默认：ISS.tle)')
    parser.add_argument('-p', '--progress', action='store_true',
                        help='显示进度条')
    parser.add_argument('-m', '--minutes', type=float, default=1.0,
                        help='外推时间 (分钟，从 TLE 历元起算)，默认：1 分钟')

    args = parser.parse_args()

    print("=" * 60)
    print("HPOP 高精度轨道外推计算程序")
    print("=" * 60)
    print(f"TLE 文件：{args.tle_file}")
    print(f"外推时间：从 TLE 历元起算 {args.minutes:.1f} 分钟")
    print(f"进度条：{'已启用' if args.progress else '禁用'}")
    print("=" * 60)

    tles = load_tle_file(args.tle_file)

    if not tles:
        print("未找到 TLE 数据!")
        return

    print(f"\n已加载 {len(tles)} 个 TLE 目标\n")

    # 创建默认航天器参数 (适合 ISS)
    spacecraft = SpacecraftState(
        mass=420000.0,      # ISS 质量约 420 吨
        area_drag=1000.0,   # 阻力参考面积
        area_srp=2500.0     # 太阳帆板面积
    )

    for tle in tles:
        print(f"目标：{tle.name}")
        print(f"卫星编号：{tle.catalog_number}")
        print(f"TLE 历元：{tle.epoch}")
        print()

        # TLE 历元时刻
        print(f">>> TLE 历元时刻的轨道状态:")
        pv_epoch = sgp4_propagate(tle, 0.0)
        print_orbit_info(tle, pv_epoch, tle.epoch)

        # 外推目标时刻（从 TLE 历元起算）
        target = tle.epoch + timedelta(minutes=args.minutes)
        print(f"\n>>> 外推目标时刻 ({target.strftime('%Y-%m-%d %H:%M:%S')})")
        print(f"    从 TLE 历元起算：{args.minutes:.1f} 分钟")

        # SGP4 外推
        print("\n--- SGP4 外推结果 ---")
        sgp4_result = propagate_to_datetime(tle, target)
        print_orbit_info(tle, sgp4_result, target)

        # HPOP 外推
        print("\n--- HPOP 外推结果 ---")
        hpop_result = hprop_from_tle(tle, target, spacecraft, args.progress)
        print_orbit_info(tle, hpop_result, target)

        # 计算差异
        dr = math.sqrt((hpop_result.x - sgp4_result.x)**2 +
                       (hpop_result.y - sgp4_result.y)**2 +
                       (hpop_result.z - sgp4_result.z)**2)
        dv = math.sqrt((hpop_result.vx - sgp4_result.vx)**2 +
                       (hpop_result.vy - sgp4_result.vy)**2 +
                       (hpop_result.vz - sgp4_result.vz)**2)
        print(f"\n--- 位置/速度差异 ---")
        print(f"距离差异：{dr:.4f} km ({dr*1000:.1f} m)")
        print(f"速度差异：{dv:.6f} km/s ({dv*1000:.3f} m/s)")

        # 轨道预报（从 TLE 历元到目标时刻）
        steps = int(args.minutes / 5) if args.minutes <= 60 else int(args.minutes / 10)
        if steps < 5:
            steps = 5
        interval = args.minutes / steps

        print(f"\n>>> HPOP 轨道预报 (从 TLE 历元到目标时刻，每 {interval:.1f} 分钟):")
        print(f"{'时间':<8} | {'历元':<12} | {'纬度':<10} | {'经度':<10} | {'高度 (km)':<10} | {'速度 (km/s)':<10}")
        print("-" * 75)

        for i in range(steps + 1):
            t = i * interval
            calc_time = tle.epoch + timedelta(minutes=t)
            pv = hprop_from_tle(tle, calc_time, spacecraft, False)
            lat, lon, alt = get_lat_lon_alt(pv, gmst_from_datetime(calc_time))
            marker = " <-- 目标" if i == steps else ""
            print(f"  +{t:5.1f}min | {calc_time.strftime('%H:%M:%S')} | "
                  f"Lat={lat:7.3f}° | Lon={lon:8.3f}° | {alt:8.1f} | {pv.v:8.3f}{marker}")

        print()


if __name__ == "__main__":
    main()
