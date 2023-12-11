<a href="https://pkg.go.dev/github.com/hebl/gofa"><img src="https://pkg.go.dev/badge/github.com/hebl/gofa.svg" alt="Go Reference"></a>

[英文](README.md)

# GoFA (基础天文标准库，Golang standards Of Fundamental Astronomy)

GoFA 是一个纯Go实现的基础天文计算库，基于国际天文联合会官方的`基础天文标准库`【`Standards Of Fundamental Astronomy (SOFA)`, [http://iausofa.org](http://iausofa.org)】
实现。

GoFA 开发文档： <https://pkg.go.dev/github.com/hebl/gofa>。

GoFA 当前版本为`v1.19`，与SOFA `v19`([2023-10-11](http://iausofa.org/2023_1011_C/))
保持同步。

官方文档中的示例程序在目录[examples](examples)下。

所有的函数已经通过测试，[测试](test)程序参考自`t_sofa_c.c`。

## 安装

```shell
go get -u github.com/hebl/gofa
```

或者

```shell
go install github.com/hebl/gofa@latest
```

## 函数

### 向量/矩阵函数库 [vml.go](vml.go)

- 初始化 (4)
- 拷贝/扩展等 (5)
- 构建旋转函数 (3)
- 向量操作 (17)
- 矩阵操作 (2)
- 矩阵-向量乘法 (4)
- 旋转向量 (2)

### 角度计算 [angle.go](angle.go)

- 球面/笛卡儿坐标转换 (6)
- 球面角距离 (4)
- 角度运算 (8)

### 天体测量 (38) [astrometry.go](astrometry.go)

### 日历 (7) [jd.go](jd.go)

### 时间尺度 (20) [ts.go](ts.go)

### 坐标系 [coord.go](coord.go)

- 黄道坐标系 (6)
- 银道坐标系 (2)
- 地平/赤道坐标系 (3)
- 地心坐标与大地坐标转换 (5)

### 地球转动与恒星时 (15) [erast.go](erast.go)

### 星历 (3) [ephem.go](ephem.go)

### 基础参数 (14) [fundargs.go](fundargs.go)

### 球心投影 (6) [projection.go](projection.go)

### 岁差/章动/极移 (64) [pn.go](pn.go)

### 星表转换 (9) [catalog.go](catalog.go)

## 版本

### v1.19

Version 1.19 包含有192个天文函数（包括时间、地球转动、恒星时、星历、岁差、岁差/章动、天体测量和坐标转换等），55个向量/矩阵函数。

## 许可证

- [SOFA License](sofa_copyr.txt)
- [MIT License](LICENSE)

2023-10-23

