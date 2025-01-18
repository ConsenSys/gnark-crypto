package mimc

import (
	"fmt"
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	conf.Package = "mimc"
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "mimc.go"), Templates: []string{"mimc.go.tmpl"}},
		{File: filepath.Join(baseDir, "options.go"), Templates: []string{"options.go.tmpl"}},
	}
	entriesTest := []bavard.Entry{
		{File: filepath.Join(baseDir, "mimc_test.go"), Templates: []string{"tests/mimc_test.go.tmpl"}},
	}

	if err := bgen.Generate(conf, conf.Package, "./crypto/hash/mimc/template", entries...); err != nil {
		return fmt.Errorf("generate package: %w", err)
	}
	if err := bgen.Generate(conf, "mimc_test", "./crypto/hash/mimc/template", entriesTest...); err != nil {
		return fmt.Errorf("generate tests: %w", err)
	}
	return nil
}
