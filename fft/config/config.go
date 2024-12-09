package config

import "strings"

type Config struct {

	// TODO this should be in the finite field package API
	GeneratorFullMultiplicativeGroup string // generator of \mathbb{F}_r^{*}

	// TODO should be generated by goff
	GeneratorMaxTwoAdicSubgroup string // generator of the maximum subgroup of size 2^<something>

	// TODO should be generated by goff
	LogTwoOrderMaxTwoAdicSubgroup string // log_2 of the max order of the max two adic subgroup

	// Directory path of the fodler where the files are generated
	OutputDir string

	// Package name of the generated package
	Package string

	// ImportPathFiniteField path to the finite field package
	FieldPackagePath string

	// FF the name of the package corresponding to the finite field
	FF string
}

// NewFFTConfig returns a data structure with needed information to generate apis for the FFT
func NewConfig(genFullMultiplicativeGroup,
	generatorMaxTwoAdicSubgroup,
	logTwoOrderMaxTwoAdicSubgroup,
	outputDir,
	fieldPackagePath string) Config {
	splittedPath := strings.Split(fieldPackagePath, "/")
	finiteFieldPackage := splittedPath[len(splittedPath)-1]
	return Config{
		GeneratorFullMultiplicativeGroup: genFullMultiplicativeGroup,
		GeneratorMaxTwoAdicSubgroup:      generatorMaxTwoAdicSubgroup,
		OutputDir:                        outputDir,
		LogTwoOrderMaxTwoAdicSubgroup:    logTwoOrderMaxTwoAdicSubgroup,
		FieldPackagePath:                 fieldPackagePath,
		FF:                               finiteFieldPackage,
	}
}
