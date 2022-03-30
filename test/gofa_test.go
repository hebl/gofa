// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"fmt"
	"math"
	"testing"
)

// Test from t_sofa.c

const (
	EPS7  = 1e-7
	EPS12 = 1e-12
)

func vvd(t *testing.T, val, want, eps float64, fun, param string) {
	var a, f float64
	a = val - want
	if a != 0.0 && math.Abs(a) > math.Abs(eps) {
		t.Fail()
		f = math.Abs(want / a)
		t.Errorf("%s failed: %s want %v got %v ((1/%.3f)) \n", fun, param, want, val, f)
	} else {
		f = math.Abs(want / a)
		t.Logf("%s passed: %s want %v got %v ((1/%.3f)) \n", fun, param, want, val, f)
	}
}

func viv(t *testing.T, val, want int, fun, param string) {
	if val != want {
		t.Errorf("%s failed: %s want %v got %v\n", fun, param, want, val)
	} else {
		t.Logf("%s passed: %s want %v got %v\n", fun, param, want, val)
	}
}

func vbv(t *testing.T, val, want byte, fun, param string) {
	if val != want {
		t.Errorf("%s failed: %s want %s got %s\n", fun, param, string(want), string(val))
	} else {
		t.Logf("%s passed: %s want %s got %s\n", fun, param, string(want), string(val))
	}
}

func vrv(t *testing.T, val, want [3][3]float64, eps float64, fun string) {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			vvd(t, val[i][j], want[i][j], eps, fun, fmt.Sprintf("[%d][%d]", i+1, j+1))
		}
	}
}

func vpv(t *testing.T, val, want [2][3]float64, eps float64, fun string) {
	for i := 0; i < 2; i++ {
		for j := 0; j < 3; j++ {
			vvd(t, val[i][j], want[i][j], eps, fun, fmt.Sprintf("[%d][%d]", i+1, j+1))
		}
	}
}
