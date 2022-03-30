// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestEpv00(t *testing.T) {
	var pvh, pvb [2][3]float64
	var j int

	j = gofa.Epv00(2400000.5, 53411.52501161, &pvh, &pvb)

	vvd(t, pvh[0][0], -0.7757238809297706813, 1e-14, "Epv00", "ph(x)")
	vvd(t, pvh[0][1], 0.5598052241363340596, 1e-14, "Epv00", "ph(y)")
	vvd(t, pvh[0][2], 0.2426998466481686993, 1e-14, "Epv00", "ph(z)")

	vvd(t, pvh[1][0], -0.1091891824147313846e-1, 1e-15, "Epv00", "vh(x)")
	vvd(t, pvh[1][1], -0.1247187268440845008e-1, 1e-15, "Epv00", "vh(y)")
	vvd(t, pvh[1][2], -0.5407569418065039061e-2, 1e-15, "Epv00", "vh(z)")

	vvd(t, pvb[0][0], -0.7714104440491111971, 1e-14, "Epv00", "pb(x)")
	vvd(t, pvb[0][1], 0.5598412061824171323, 1e-14, "Epv00", "pb(y)")
	vvd(t, pvb[0][2], 0.2425996277722452400, 1e-14, "Epv00", "pb(z)")

	vvd(t, pvb[1][0], -0.1091874268116823295e-1, 1e-15, "Epv00", "vb(x)")
	vvd(t, pvb[1][1], -0.1246525461732861538e-1, 1e-15, "Epv00", "vb(y)")
	vvd(t, pvb[1][2], -0.5404773180966231279e-2, 1e-15, "Epv00", "vb(z)")

	viv(t, j, 0, "Epv00", "j")
}

func TestMoon98(t *testing.T) {
	var pv [2][3]float64

	gofa.Moon98(2400000.5, 43999.9, &pv)

	wpv := [2][3]float64{{-0.2601295959971044180e-2, 0.6139750944302742189e-3, 0.2640794528229828909e-3}, {-0.1244321506649895021e-3, -0.5219076942678119398e-3, -0.1716132214378462047e-3}}
	vpv(t, pv, wpv, 1e-11, "Moon98")
}

func TestPlan94(t *testing.T) {
	var pv [2][3]float64
	var j int

	j = gofa.Plan94(2400000.5, 1e6, 0, &pv)

	vvd(t, pv[0][0], 0.0, 0.0, "Plan94", "x 1")
	vvd(t, pv[0][1], 0.0, 0.0, "Plan94", "y 1")
	vvd(t, pv[0][2], 0.0, 0.0, "Plan94", "z 1")

	vvd(t, pv[1][0], 0.0, 0.0, "Plan94", "xd 1")
	vvd(t, pv[1][1], 0.0, 0.0, "Plan94", "yd 1")
	vvd(t, pv[1][2], 0.0, 0.0, "Plan94", "zd 1")

	viv(t, j, -1, "Plan94", "j 1")

	j = gofa.Plan94(2400000.5, 1e6, 10, &pv)

	viv(t, j, -1, "Plan94", "j 2")

	j = gofa.Plan94(2400000.5, -320000, 3, &pv)

	vvd(t, pv[0][0], 0.9308038666832975759, 1e-11, "Plan94", "x 3")
	vvd(t, pv[0][1], 0.3258319040261346000, 1e-11, "Plan94", "y 3")
	vvd(t, pv[0][2], 0.1422794544481140560, 1e-11, "Plan94", "z 3")

	vvd(t, pv[1][0], -0.6429458958255170006e-2, 1e-11, "Plan94", "xd 3")
	vvd(t, pv[1][1], 0.1468570657704237764e-1, 1e-11, "Plan94", "yd 3")
	vvd(t, pv[1][2], 0.6406996426270981189e-2, 1e-11, "Plan94", "zd 3")

	viv(t, j, 1, "Plan94", "j 3")

	j = gofa.Plan94(2400000.5, 43999.9, 1, &pv)

	vvd(t, pv[0][0], 0.2945293959257430832, 1e-11, "Plan94", "x 4")
	vvd(t, pv[0][1], -0.2452204176601049596, 1e-11, "Plan94", "y 4")
	vvd(t, pv[0][2], -0.1615427700571978153, 1e-11, "Plan94", "z 4")

	vvd(t, pv[1][0], 0.1413867871404614441e-1, 1e-11, "Plan94", "xd 4")
	vvd(t, pv[1][1], 0.1946548301104706582e-1, 1e-11, "Plan94", "yd 4")
	vvd(t, pv[1][2], 0.8929809783898904786e-2, 1e-11, "Plan94", "zd 4")

	viv(t, j, 0, "Plan94", "j 4")
}