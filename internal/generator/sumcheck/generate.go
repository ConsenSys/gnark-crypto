package sumcheck

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

func Generate(conf config.FieldDependency, baseDir string, bgen *bavard.BatchGenerator) error {
	if !git.HasChanges("./sumcheck/template/") {
		return nil
	}
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "sumcheck.go"), Templates: []string{"sumcheck.go.tmpl"}},
		{File: filepath.Join(baseDir, "sumcheck_test.go"), Templates: []string{"sumcheck.test.go.tmpl"}},
	}
	return bgen.Generate(conf, "sumcheck", "./sumcheck/template/", entries...)
}
