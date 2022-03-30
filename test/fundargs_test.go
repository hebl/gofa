// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestFad03(t *testing.T) {
	vvd(t, gofa.Fad03(0.80), 1.946709205396925672, 1e-12, "Fad03", "")
}
func TestFae03(t *testing.T) {
	vvd(t, gofa.Fae03(0.80), 1.744713738913081846, 1e-12, "Fae03", "")
}
func TestFaf03(t *testing.T) {
	vvd(t, gofa.Faf03(0.80), 0.2597711366745499518, 1e-12, "Faf03", "")
}
func TestFaju03(t *testing.T) {
	vvd(t, gofa.Faju03(0.80), 5.275711665202481138, 1e-12, "Faju03", "")
}
func TestFal03(t *testing.T) {
	vvd(t, gofa.Fal03(0.80), 5.132369751108684150, 1e-12, "Fal03", "")
}
func TestFalp03(t *testing.T) {
	vvd(t, gofa.Falp03(0.80), 6.226797973505507345, 1e-12, "Falp03", "")
}
func TestFama03(t *testing.T) {
	vvd(t, gofa.Fama03(0.80), 3.275506840277781492, 1e-12, "Fama03", "")
}
func TestFame03(t *testing.T) {
	vvd(t, gofa.Fame03(0.80), 5.417338184297289661, 1e-12, "Fame03", "")
}
func TestFane03(t *testing.T) {
	vvd(t, gofa.Fane03(0.80), 2.079343830860413523, 1e-12, "Fane03", "")
}
func TestFaom03(t *testing.T) {
	vvd(t, gofa.Faom03(0.80), -5.973618440951302183, 1e-12, "Faom03", "")
}
func TestFapa03(t *testing.T) {
	vvd(t, gofa.Fapa03(0.80), 0.1950884762240000000e-1, 1e-12, "Fapa03", "")

}
func TestFasa03(t *testing.T) {
	vvd(t, gofa.Fasa03(0.80), 5.371574539440827046, 1e-12, "Fasa03", "")
}
func TestFaur03(t *testing.T) {
	vvd(t, gofa.Faur03(0.80), 5.180636450180413523, 1e-12, "Faur03", "")
}
func TestFave03(t *testing.T) {
	vvd(t, gofa.Fave03(0.80), 3.424900460533758000, 1e-12, "Fave03", "")
}
