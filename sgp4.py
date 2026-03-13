"""
SGP4 轨道外推算法实现
Simplified General Perturbations 4 模型，用于 TLE 数据轨道计算

基于 Spacetrack Report #3 和 CelesTrak 参考实现
"""

import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Tuple, List, Optional

# ===== 常数定义 =====
SECONDS_PER_DAY = 86400.0
SECONDS_PER_MINUTE = 60.0
MINUTES_PER_DAY = 1440.0

# 地球常数 (WGS72 - SGP4 标准)
R_EARTH = 6378.135  # km - 地球半径
MU_EARTH = 398600.8  # km^3/s^2 - 地球重力常数
J2 = 1.082616e-3    # 二阶带谐系数
J3 = -2.53881e-6    # 三阶带谐系数
J4 = -1.65597e-6    # 四阶带谐系数

# 导出常数
QOMS2T = 1.88027916e-9
CK2 = 0.5 * J2 * R_EARTH * R_EARTH
CK4 = -0.375 * J4 * R_EARTH ** 4

# 角度转换
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi

# 地球自转角速度 (rad/min)
THETA_EARTH = 7.2921151467e-5 * 60.0


@dataclass
class TLE:
    """TLE 数据结构"""
    name: str
    catalog_number: int
    classification: str
    epoch_year: int
    epoch_day: float
    mean_motion_derivative: float
    mean_motion_second_derivative: float
    bstar: float
    inclination: float
    raan: float
    eccentricity: float
    arg_perigee: float
    mean_anomaly: float
    mean_motion: float
    orbit_number: int

    @property
    def epoch(self) -> datetime:
        """获取 TLE 历元时间"""
        year = self.epoch_year
        if year >= 57:
            year = 1900 + year
        else:
            year = 2000 + year

        day_fraction = self.epoch_day
        day = int(day_fraction)
        fraction = day_fraction - day

        epoch = datetime(year, 1, 1) + timedelta(days=day - 1, seconds=fraction * SECONDS_PER_DAY)
        return epoch


@dataclass
class PositionVelocity:
    """位置和速度 (ECI 坐标系)"""
    x: float  # km
    y: float  # km
    z: float  # km
    vx: float  # km/s
    vy: float  # km/s
    vz: float  # km/s

    @property
    def r(self) -> float:
        """位置矢量模 (km)"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def v(self) -> float:
        """速度矢量模 (km/s)"""
        return math.sqrt(self.vx**2 + self.vy**2 + self.vz**2)


def parse_tle_number(line: str, start: int, end: int) -> float:
    """解析 TLE 中的数字字段"""
    s = line[start:end].strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def parse_tle_bstar(line: str) -> float:
    """解析 TLE B* 阻力系数

    B* 格式：nnnnn±d 表示 0.nnnnn × 10^±d
    """
    if len(line) < 61:
        return 0.0

    bstar_raw = line[53:61].strip()
    if not bstar_raw or len(bstar_raw) < 6:
        return 0.0

    try:
        mantissa = float(bstar_raw[:5]) / 100000.0
        exp = int(bstar_raw[5:])
        return mantissa * (10 ** exp)
    except ValueError:
        return 0.0


def parse_tle_line1(line: str) -> dict:
    """解析 TLE 第一行"""
    if len(line) < 65:
        return {}

    epoch_str = line[18:32].strip()
    if len(epoch_str) >= 2:
        epoch_year = int(epoch_str[:2])
        epoch_day = float(epoch_str[2:])
    else:
        epoch_year = 0
        epoch_day = 0.0

    return {
        'catalog_number': int(line[2:7].strip()) if len(line) > 7 else 0,
        'classification': line[7].strip() if len(line) > 7 else 'U',
        'epoch_year': epoch_year,
        'epoch_day': epoch_day,
        'mean_motion_derivative': parse_tle_number(line, 33, 43),
        'mean_motion_second_derivative': parse_tle_number(line, 44, 52),
        'bstar': parse_tle_bstar(line),
        'checksum1': int(line[68]) if len(line) > 68 else 0,
    }


def parse_tle_line2(line: str) -> dict:
    """解析 TLE 第二行"""
    if len(line) < 68:
        return {}

    return {
        'catalog_number': int(line[2:7].strip()),
        'inclination': parse_tle_number(line, 8, 16),
        'raan': parse_tle_number(line, 17, 25),
        'eccentricity': int(line[26:33].strip() or '0') / 1e7,
        'arg_perigee': parse_tle_number(line, 34, 42),
        'mean_anomaly': parse_tle_number(line, 43, 51),
        'mean_motion': parse_tle_number(line, 52, 63),
        'orbit_number': int(line[63:68].strip()),
    }


def parse_tle(line1: str, line2: str, name: str = "") -> TLE:
    """解析 TLE 数据"""
    if name == "" and len(line1) > 0 and not line1.startswith('1'):
        name = line1.strip()
        line1, line2 = line2, line1

    data1 = parse_tle_line1(line1)
    data2 = parse_tle_line2(line2)

    return TLE(
        name=name,
        catalog_number=data1.get('catalog_number', 0),
        classification=data1.get('classification', 'U'),
        epoch_year=data1.get('epoch_year', 0),
        epoch_day=data1.get('epoch_day', 0.0),
        mean_motion_derivative=data1.get('mean_motion_derivative', 0.0),
        mean_motion_second_derivative=data1.get('mean_motion_second_derivative', 0.0),
        bstar=data1.get('bstar', 0.0),
        inclination=data2.get('inclination', 0.0),
        raan=data2.get('raan', 0.0),
        eccentricity=data2.get('eccentricity', 0.0),
        arg_perigee=data2.get('arg_perigee', 0.0),
        mean_anomaly=data2.get('mean_anomaly', 0.0),
        mean_motion=data2.get('mean_motion', 0.0),
        orbit_number=data2.get('orbit_number', 0),
    )


def load_tle_file(filepath: str) -> List[TLE]:
    """从文件加载 TLE 数据"""
    tles = []
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = [line.rstrip('\n\r') for line in f if line.strip()]

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('1'):
            if i + 1 < len(lines) and lines[i+1].startswith('2'):
                tles.append(parse_tle(line, lines[i+1], ""))
                i += 2
        elif not line.startswith('2'):
            if i + 2 < len(lines) and lines[i+1].startswith('1') and lines[i+2].startswith('2'):
                tles.append(parse_tle(lines[i+1], lines[i+2], name=line.strip()))
                i += 3
            else:
                i += 1
        else:
            i += 1
    return tles


def normalize_angle(angle: float) -> float:
    """将角度归一化到 [0, 2π)"""
    twopi = 2 * math.pi
    angle = angle % twopi
    if angle < 0:
        angle += twopi
    return angle


class SGP4Propagator:
    """
    SGP4 轨道外推器

    基于 Spacetrack Report #3 的简化实现
    """

    def __init__(self, tle: TLE):
        self.tle = tle

        # 角度转弧度
        self.i0 = tle.inclination * DEG2RAD
        self.raan0 = tle.raan * DEG2RAD
        self.argp0 = tle.arg_perigee * DEG2RAD
        self.M0 = tle.mean_anomaly * DEG2RAD

        # 轨道参数
        self.e0 = tle.eccentricity
        # n0 单位：rad/min
        self.n0 = tle.mean_motion * 2 * math.pi / MINUTES_PER_DAY
        self.bstar = tle.bstar

        # 平均运动导数
        self.ndot = tle.mean_motion_derivative * 2 * math.pi / (MINUTES_PER_DAY ** 2)
        self.nddot = tle.mean_motion_second_derivative * 2 * math.pi / (MINUTES_PER_DAY ** 3)

        self.epoch = tle.epoch
        self._init_constants()

    def _init_constants(self):
        """初始化 SGP4 常数"""
        # 半长轴计算
        # 注意单位：n0 是 rad/min，需要转换为 rad/s 或者使用 km^3/min^2 单位的 mu
        # MU_EARTH = 398600.8 km^3/s^2
        # 转换为 km^3/min^2: MU_min = MU_EARTH * 3600
        MU_min = MU_EARTH * 3600.0  # km^3/min^2
        a1 = (MU_min / self.n0**2)**(1/3)

        # 平均运动相关的常数
        self.cos_i0 = math.cos(self.i0)
        self.sin_i0 = math.sin(self.i0)

        # 偏心率辅助量
        beta0 = math.sqrt(1 - self.e0**2)

        # 半通径
        p = a1 * (1 - self.e0**2)
        self.p = p

        # 二阶摄动系数
        k2 = CK2
        a3ovk2 = -J3 / CK2 * R_EARTH**3

        # 辅助量
        theta2 = self.cos_i0**2
        x3thm1 = 3 * theta2 - 1
        x1mth2 = 1 - theta2
        x7thm1 = 7 * theta2 - 1

        # 长期摄动率
        # RAAN 变化率 (rad/min)
        dRAAN_dt = -1.5 * self.n0 * k2 / (p**2) * self.cos_i0

        # 近地点幅角变化率
        dargp_dt = 0.75 * self.n0 * k2 / (p**2) * (4 * theta2 - 1)

        # 平近点角变化率
        dM_dt = self.n0 + 0.75 * self.n0 * k2 / (p**2) * beta0 * x3thm1

        self.dRAAN_dt = dRAAN_dt
        self.dargp_dt = dargp_dt
        self.dM_dt = dM_dt

        # 存储半长轴
        self.a = a1

    def propagate(self, tsince: float) -> PositionVelocity:
        """
        SGP4 外推

        参数:
            tsince: 从 TLE 历元起算的时间 (分钟)

        返回:
            PositionVelocity: ECI 坐标系中的位置和速度
        """
        # 长期摄动下的轨道根数
        raan = normalize_angle(self.raan0 + self.dRAAN_dt * tsince)
        argp = normalize_angle(self.argp0 + self.dargp_dt * tsince)
        M = normalize_angle(self.M0 + self.dM_dt * tsince)

        # 平近点角长期项
        # 对于深空轨道需要更多摄动项，这里使用简化版本

        # 使用简化模型，不考虑阻力引起的半长轴和偏心率变化
        # 这对短时间预报（几天内）是足够的
        a = self.a
        e = self.e0

        # 求解开普勒方程
        E = self._solve_kepler(M, e)

        # 计算位置和速度
        return self._rv_from_kepler(a, e, E, self.i0, raan, argp)

    def _solve_kepler(self, M: float, e: float, tol: float = 1e-10,
                      max_iter: int = 50) -> float:
        """求解开普勒方程 M = E - e*sin(E)"""
        M = normalize_angle(M)
        if M > math.pi:
            M -= 2 * math.pi

        if e < 0.8:
            E = M
        else:
            E = math.pi

        for _ in range(max_iter):
            f = E - e * math.sin(E) - M
            f_prime = 1 - e * math.cos(E)

            if abs(f_prime) < 1e-12:
                break

            delta = -f / f_prime
            E += delta

            if abs(delta) < tol:
                break

        return E

    def _rv_from_kepler(self, a: float, e: float, E: float,
                        i: float, raan: float, argp: float) -> PositionVelocity:
        """从开普勒轨道根数计算位置和速度 (ECI 坐标系)

        使用标准的轨道力学公式
        """
        cos_E = math.cos(E)
        sin_E = math.sin(E)

        # 距离
        r = a * (1 - e * cos_E)

        # 轨道平面坐标 (perifocal frame)
        # x 指向近地点，y 在轨道平面内垂直于 x
        x_orb = a * (cos_E - e)
        y_orb = a * math.sqrt(1 - e**2) * sin_E

        # 速度分量 (perifocal)
        # v = sqrt(mu/a) / (1 - e*cos_E) * (-sin_E, sqrt(1-e^2)*cos_E)
        sqrt_mu_a = math.sqrt(MU_EARTH / a)
        xdot_orb = -sqrt_mu_a * sin_E / (1 - e * cos_E)
        ydot_orb = sqrt_mu_a * math.sqrt(1 - e**2) * cos_E / (1 - e * cos_E)

        # 旋转矩阵从 perifocal 到 ECI
        # R = R3(-raan) * R1(-i) * R3(-argp)
        cos_raan = math.cos(raan)
        sin_raan = math.sin(raan)
        cos_argp = math.cos(argp)
        sin_argp = math.sin(argp)
        cos_i = math.cos(i)
        sin_i = math.sin(i)

        # 组合旋转矩阵元素
        # r_ECI = R * r_perifocal
        r11 = cos_raan * cos_argp - sin_raan * sin_argp * cos_i
        r12 = -cos_raan * sin_argp - sin_raan * cos_argp * cos_i
        r13 = sin_raan * sin_i

        r21 = sin_raan * cos_argp + cos_raan * sin_argp * cos_i
        r22 = -sin_raan * sin_argp + cos_raan * cos_argp * cos_i
        r23 = -cos_raan * sin_i

        r31 = sin_argp * sin_i
        r32 = cos_argp * sin_i
        r33 = cos_i

        # 位置 ECI
        x = r11 * x_orb + r12 * y_orb
        y = r21 * x_orb + r22 * y_orb
        z = r31 * x_orb + r32 * y_orb

        # 速度 ECI
        vx = r11 * xdot_orb + r12 * ydot_orb
        vy = r21 * xdot_orb + r22 * ydot_orb
        vz = r31 * xdot_orb + r32 * ydot_orb

        return PositionVelocity(x, y, z, vx, vy, vz)


def gmst_from_datetime(dt: datetime) -> float:
    """计算格林尼治恒星时 (rad)"""
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

    T = (JD - 2451545.0) / 36525.0
    gmst = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + 0.000387933 * T**2
    gmst = gmst % 360.0

    return gmst * DEG2RAD


def get_lat_lon_alt(pv: PositionVelocity, gmst: float) -> Tuple[float, float, float]:
    """计算星下点纬度和经度"""
    r = pv.r
    if r < 1e-10:
        return 0.0, 0.0, 0.0

    dec = math.asin(max(-1, min(1, pv.z / r)))
    ra = math.atan2(pv.y, pv.x)
    lon = normalize_angle(ra - gmst)
    if lon > math.pi:
        lon -= 2 * math.pi

    lat = dec
    alt = r - R_EARTH

    return lat * RAD2DEG, lon * RAD2DEG, alt


def sgp4_propagate(tle: TLE, tsince: float) -> PositionVelocity:
    """SGP4 外推函数接口"""
    prop = SGP4Propagator(tle)
    return prop.propagate(tsince)


def propagate_to_datetime(tle: TLE, target_time: datetime) -> PositionVelocity:
    """外推到指定时间"""
    delta = target_time - tle.epoch
    tsince = delta.total_seconds() / 60.0
    return sgp4_propagate(tle, tsince)


def print_orbit_info(tle: TLE, pv: PositionVelocity, timestamp: datetime):
    """打印轨道信息"""
    lat, lon, alt = get_lat_lon_alt(pv, gmst_from_datetime(timestamp))

    print(f"时间：{timestamp.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"位置 (ECI): x={pv.x:10.3f} km, y={pv.y:10.3f} km, z={pv.z:10.3f} km")
    print(f"速度 (ECI): vx={pv.vx:10.6f} km/s, vy={pv.vy:10.6f} km/s, vz={pv.vz:10.6f} km/s")
    print(f"地心距离：{pv.r:.3f} km")
    print(f"星下点：纬度={lat:8.4f}°, 经度={lon:8.4f}°")
    print(f"高度：{alt:.3f} km")
    print("-" * 60)


def main():
    """主程序"""
    import sys
    import argparse

    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='SGP4/HPOP 轨道外推计算程序 - 根据 TLE 数据计算卫星轨道位置和速度',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python sgp4.py                    # 使用 ISS.tle，外推 10 分钟后
  python sgp4.py NOAA19.tle         # 指定 TLE 文件
  python sgp4.py --hpop             # 使用 HPOP 高精度模式
  python sgp4.py --hpop -p          # HPOP 模式 + 进度条
  python sgp4.py -m 60              # 外推 60 分钟后的轨道
  python sgp4.py --hpop -m 120 -p   # HPOP 外推 2 小时，显示进度条
'''
    )

    parser.add_argument('tle_file', nargs='?', default='ISS.tle',
                        help='TLE 文件路径 (默认：ISS.tle)')
    parser.add_argument('--hpop', action='store_true',
                        help='使用 HPOP 高精度外推模式 (默认：SGP4)')
    parser.add_argument('-p', '--progress', action='store_true',
                        help='显示进度条 (仅 HPOP 模式)')
    parser.add_argument('-m', '--minutes', type=float, default=1.0,
                        help='外推时间 (分钟，从当前时刻起算)，默认：1 分钟')

    args = parser.parse_args()

    print("=" * 60)
    print("SGP4/HPOP 轨道外推计算程序")
    print("=" * 60)
    print(f"TLE 文件：{args.tle_file}")
    print(f"外推时间：当前时刻 + {args.minutes:.1f} 分钟")
    print(f"模式：{'HPOP 高精度' if args.hpop else 'SGP4 快速'}")
    if args.hpop and args.progress:
        print(f"进度条：已启用")
    print("=" * 60)

    tles = load_tle_file(args.tle_file)

    if not tles:
        print("未找到 TLE 数据!")
        return

    print(f"\n已加载 {len(tles)} 个 TLE 目标\n")

    # 如果启用 HPOP，导入 hpop 模块
    if args.hpop:
        try:
            from hpop import SpacecraftState, compare_sgp4_hpop
            spacecraft = SpacecraftState(mass=420000.0, area_drag=1000.0, area_srp=2500.0)
        except ImportError:
            print("警告：hpop 模块不可用，回退到 SGP4 模式\n")
            args.hpop = False

    for tle in tles:
        print(f"目标：{tle.name}")
        print(f"卫星编号：{tle.catalog_number}")
        print(f"TLE 历元：{tle.epoch}")
        print(f"偏心率：{tle.eccentricity:.7f}")
        print(f"轨道倾角：{tle.inclination:.4f}°")
        print(f"平均运动：{tle.mean_motion:.8f} rev/day")
        print(f"B* 阻力系数：{tle.bstar:.2e} 1/ER")
        print()

        # TLE 历元时刻
        print(f">>> TLE 历元时刻的轨道状态:")
        pv = sgp4_propagate(tle, 0.0)
        print_orbit_info(tle, pv, tle.epoch)

        # 当前时刻
        now = datetime.now()
        delta_from_epoch = (now - tle.epoch).total_seconds() / 60.0
        print(f">>> 当前时刻 ({now.strftime('%Y-%m-%d %H:%M:%S')})")
        print(f"    距离 TLE 历元：{delta_from_epoch:.1f} 分钟 ({delta_from_epoch/1440:.1f} 天)")

        # 目标时刻
        target = now + timedelta(minutes=args.minutes)
        print(f"\n>>> 外推目标时刻 ({target.strftime('%Y-%m-%d %H:%M:%S')})")
        print(f"    距离当前时刻：{args.minutes:.1f} 分钟")

        if args.hpop:
            from hpop import hprop_from_tle
            pv = hprop_from_tle(tle, target, spacecraft, args.progress)
        else:
            pv = propagate_to_datetime(tle, target)
        print_orbit_info(tle, pv, target)

        # 轨道预报（从当前时刻到目标时刻）
        steps = int(args.minutes / 5) if args.minutes <= 60 else int(args.minutes / 10)
        if steps < 5:
            steps = 5
        interval = args.minutes / steps

        print(f">>> 轨道预报 (从当前时刻到目标时刻，每 {interval:.1f} 分钟):")
        for i in range(steps + 1):
            t = i * interval
            calc_time = now + timedelta(minutes=t)
            if args.hpop:
                from hpop import hprop_from_tle
                pv = hprop_from_tle(tle, calc_time, spacecraft, False)
            else:
                pv = propagate_to_datetime(tle, calc_time)
            lat, lon, alt = get_lat_lon_alt(pv, gmst_from_datetime(calc_time))
            marker = " <-- 目标" if i == steps else ""
            print(f"  +{t:5.1f}min | {calc_time.strftime('%H:%M:%S')} | "
                  f"Lat={lat:7.3f}° Lon={lon:8.3f}° Alt={alt:7.1f}km | V={pv.v:.3f} km/s{marker}")

        print()


if __name__ == "__main__":
    main()
