//go:build !noadx
// +build !noadx

// Copyright 2020 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

package fptower

import "golang.org/x/sys/cpu"

// supportAdx will be set only on amd64 that has MULX and ADDX instructions
var (
	supportAdx = cpu.X86.HasADX && cpu.X86.HasBMI2
	_          = supportAdx // used in asm
)
