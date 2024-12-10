// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package mimc

import (
	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
)

// Option defines option for altering the behavior of the MiMC hasher.
// See the descriptions of functions returning instances of this type for
// particular options.
type Option func(*mimcConfig)

type mimcConfig struct {
	byteOrder fr.ByteOrder
}

// default options
func mimcOptions(opts ...Option) mimcConfig {
	// apply options
	opt := mimcConfig{
		byteOrder: fr.BigEndian,
	}
	for _, option := range opts {
		option(&opt)
	}
	return opt
}

// WithByteOrder sets the byte order used to decode the input
// in the Write method. Default is BigEndian.
func WithByteOrder(byteOrder fr.ByteOrder) Option {
	return func(opt *mimcConfig) {
		opt.byteOrder = byteOrder
	}
}
