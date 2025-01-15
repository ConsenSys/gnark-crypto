package gkr

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

type Config struct {
	config.FieldDependency
	GenerateTests           bool
	RetainTestCaseRawInfo   bool
	OutsideGkrPackage       bool
	TestVectorsRelativePath string
}

func Generate(config Config, baseDir string, bgen *bavard.BatchGenerator) error {
	if !git.HasChanges("./gkr/template/") {
		return nil
	}
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "gkr.go"), Templates: []string{"gkr.go.tmpl"}},
	}

	if config.GenerateTests {
		entries = append(entries,
			bavard.Entry{File: filepath.Join(baseDir, "gkr_test.go"), Templates: []string{"gkr.test.go.tmpl", "gkr.test.vectors.go.tmpl"}})
	}

	return bgen.Generate(config, "gkr", "./gkr/template/", entries...)
}
